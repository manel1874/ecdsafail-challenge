use std::collections::{BTreeMap, BTreeSet};

use quantum_ecc::circuit::{Op, OperationType, NO_BIT, NO_QUBIT};
use quantum_ecc::point_add;

#[derive(Clone, Debug, Eq, Ord, PartialEq, PartialOrd)]
struct CondSig(Vec<u64>);

#[derive(Clone, Debug, Eq, Ord, PartialEq, PartialOrd)]
struct FaceKey {
    cond: CondSig,
    a: u64,
    b: u64,
}

#[derive(Clone, Debug)]
struct FaceBucket {
    face: FaceKey,
    ops: Vec<usize>,
    third_parity: BTreeSet<u64>,
}

fn main() {
    let target_face = requested_face();
    if target_face.is_some() && std::env::var_os("TRACE_OP_SITES").is_none() {
        std::env::set_var("TRACE_OP_SITES", "1");
    }

    let ops = point_add::build();
    let sites = point_add::take_last_op_sites();
    let ccx = ops
        .iter()
        .filter(|op| op.kind == OperationType::CCX)
        .count();
    let ccz = ops
        .iter()
        .filter(|op| op.kind == OperationType::CCZ)
        .count();
    println!("ops={} ccx={} ccz={}", ops.len(), ccx, ccz);

    let (exact_save, exact_buckets) = diagonal_run_face_savings(&ops);
    let (global_save, global_buckets) = face_savings_for_slice(&ops, 0);
    let conditioned = conditioned_parity_candidates(&ops);
    print_diagonal_run_summary(&ops, &sites, 24);

    println!(
        "literal_ccz_shared_face_exact_diagonal_run_save={}",
        exact_save
    );
    print_buckets("exact", &exact_buckets, 20);

    println!(
        "literal_ccz_shared_face_strict_condition_global_save={}",
        global_save
    );
    print_buckets("strict_global", &global_buckets, 20);
    print_conditioned_candidates(&conditioned, 20);
    if std::env::var("TODD_SITE_COUNTS").ok().as_deref() == Some("1") {
        print_ccz_site_counts(&ops, &sites, 40);
    }
    if let Some(site_filter) = std::env::var("TODD_GROUP_SITE").ok() {
        print_ordered_site_groups(&ops, &sites, &site_filter, 40);
    }
    if let Some(site_filter) = std::env::var("TODD_FACE_FACTOR_SITE").ok() {
        print_face_factor_segments(&ops, &sites, &site_filter, 40);
    }

    if let Some((a, b)) = target_face {
        inspect_face(&ops, &sites, a, b);
    }
}

fn diagonal_run_face_savings(ops: &[Op]) -> (usize, Vec<FaceBucket>) {
    let mut total = 0usize;
    let mut all_buckets = Vec::new();
    let mut start = 0usize;
    while start < ops.len() {
        while start < ops.len() && !is_diagonal_phase_op(ops[start].kind) {
            start += 1;
        }
        let mut end = start;
        while end < ops.len() && is_diagonal_phase_op(ops[end].kind) {
            end += 1;
        }
        if start < end {
            let (save, mut buckets) = face_savings_for_slice(&ops[start..end], start);
            total += save;
            all_buckets.append(&mut buckets);
        }
        start = end.saturating_add(1);
    }
    all_buckets.sort_by(|a, b| bucket_save(b).cmp(&bucket_save(a)));
    (total, all_buckets)
}

fn face_savings_for_slice(ops: &[Op], base: usize) -> (usize, Vec<FaceBucket>) {
    let mut buckets: BTreeMap<FaceKey, FaceBucket> = BTreeMap::new();
    let mut stack = Vec::new();
    for (local_i, op) in ops.iter().enumerate() {
        match op.kind {
            OperationType::PushCondition => {
                stack.push(op.c_condition.0);
                continue;
            }
            OperationType::PopCondition => {
                stack.pop();
                continue;
            }
            _ => {}
        }

        if op.kind != OperationType::CCZ {
            continue;
        }
        let idx = base + local_i;
        let cond = cond_sig(&stack, op.c_condition);
        let q = [op.q_control1.0, op.q_control2.0, op.q_target.0];
        for &(i, j, k) in &[(0usize, 1usize, 2usize), (0, 2, 1), (1, 2, 0)] {
            let (a, b) = if q[i] <= q[j] {
                (q[i], q[j])
            } else {
                (q[j], q[i])
            };
            let face = FaceKey {
                cond: cond.clone(),
                a,
                b,
            };
            let bucket = buckets.entry(face.clone()).or_insert_with(|| FaceBucket {
                face,
                ops: Vec::new(),
                third_parity: BTreeSet::new(),
            });
            bucket.ops.push(idx);
            if !bucket.third_parity.insert(q[k]) {
                bucket.third_parity.remove(&q[k]);
            }
        }
    }

    let mut used = BTreeSet::new();
    let mut candidates: Vec<FaceBucket> = buckets
        .into_values()
        .filter(|bucket| bucket_save(bucket) > 0)
        .collect();
    candidates.sort_by(|a, b| {
        bucket_save(b)
            .cmp(&bucket_save(a))
            .then_with(|| a.ops[0].cmp(&b.ops[0]))
    });

    let mut chosen = Vec::new();
    let mut total = 0usize;
    for bucket in candidates {
        if bucket.ops.iter().any(|idx| used.contains(idx)) {
            continue;
        }
        total += bucket_save(&bucket);
        for idx in &bucket.ops {
            used.insert(*idx);
        }
        chosen.push(bucket);
    }
    (total, chosen)
}

fn bucket_save(bucket: &FaceBucket) -> usize {
    let replacement = usize::from(!bucket.third_parity.is_empty());
    bucket.ops.len().saturating_sub(replacement)
}

fn print_buckets(label: &str, buckets: &[FaceBucket], limit: usize) {
    println!(
        "{}_bucket_count={} showing_top={}",
        label,
        buckets.len(),
        limit.min(buckets.len())
    );
    for bucket in buckets.iter().take(limit) {
        println!(
            "{} save={} ops={} cond={} face=({}, {}) xor_terms={} first_op={}",
            label,
            bucket_save(bucket),
            bucket.ops.len(),
            format_cond(&bucket.face.cond),
            bucket.face.a,
            bucket.face.b,
            bucket.third_parity.len(),
            bucket.ops[0],
        );
    }
}

fn is_diagonal_phase_op(kind: OperationType) -> bool {
    matches!(
        kind,
        OperationType::Neg | OperationType::Z | OperationType::CZ | OperationType::CCZ
    )
}

fn print_diagonal_run_summary(ops: &[Op], sites: &[point_add::OpSite], limit: usize) {
    let mut rows = Vec::<(usize, usize, usize, usize, BTreeSet<String>)>::new();
    let mut start = 0usize;
    while start < ops.len() {
        while start < ops.len() && !is_diagonal_phase_op(ops[start].kind) {
            start += 1;
        }
        let mut end = start;
        let mut ccz = 0usize;
        let mut phase_ops = 0usize;
        let mut site_set = BTreeSet::new();
        while end < ops.len() && is_diagonal_phase_op(ops[end].kind) {
            phase_ops += 1;
            if ops[end].kind == OperationType::CCZ {
                ccz += 1;
                site_set.insert(format_site(sites, end));
            }
            end += 1;
        }
        if ccz > 1 {
            rows.push((ccz, phase_ops, start, end.saturating_sub(1), site_set));
        }
        start = end.saturating_add(1);
    }
    rows.sort_by(|left, right| right.0.cmp(&left.0).then_with(|| left.2.cmp(&right.2)));
    let total_multi_ccz: usize = rows.iter().map(|row| row.0).sum();
    println!(
        "diagonal_ccz_runs multi_run_count={} multi_run_ccz={} showing_top={}",
        rows.len(),
        total_multi_ccz,
        limit.min(rows.len())
    );
    for (ccz, phase_ops, start, end, site_set) in rows.iter().take(limit) {
        let sites_preview = site_set
            .iter()
            .take(4)
            .cloned()
            .collect::<Vec<_>>()
            .join(",");
        println!(
            "diagonal_ccz_run ccz={} phase_ops={} span={}..{} sites={}{}",
            ccz,
            phase_ops,
            start,
            end,
            sites_preview,
            if site_set.len() > 4 { ",..." } else { "" },
        );
    }
}

fn cond_sig(stack: &[u64], direct: quantum_ecc::circuit::BitId) -> CondSig {
    let mut bits = stack.to_vec();
    if direct != NO_BIT {
        bits.push(direct.0);
    }
    bits.sort_unstable();
    bits.dedup();
    CondSig(bits)
}

fn format_cond(cond: &CondSig) -> String {
    if cond.0.is_empty() {
        "none".to_owned()
    } else {
        cond.0
            .iter()
            .map(u64::to_string)
            .collect::<Vec<_>>()
            .join("&")
    }
}

fn requested_face() -> Option<(u64, u64)> {
    let a = std::env::var("TODD_FACE_A").ok()?.parse::<u64>().ok()?;
    let b = std::env::var("TODD_FACE_B").ok()?.parse::<u64>().ok()?;
    Some(if a <= b { (a, b) } else { (b, a) })
}

#[derive(Clone, Debug)]
struct FaceMember {
    op_index: usize,
    third: u64,
    cond: CondSig,
}

fn inspect_face(ops: &[Op], sites: &[point_add::OpSite], a: u64, b: u64) {
    let requested_cond = std::env::var("TODD_FACE_COND").ok();
    let first_filter = std::env::var("TODD_FACE_FIRST")
        .ok()
        .and_then(|s| s.parse::<usize>().ok());
    let print_limit = std::env::var("TODD_FACE_LIMIT")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(48);

    let mut members = face_members(ops, a, b);
    if let Some(cond) = requested_cond.as_deref() {
        members.retain(|m| format_cond(&m.cond) == cond);
    }
    if let Some(first) = first_filter {
        let Some(cond) = members
            .iter()
            .find(|m| m.op_index == first)
            .map(|m| m.cond.clone())
        else {
            println!(
                "inspect_face face=({}, {}) first_op={} not present after cond filter",
                a, b, first
            );
            return;
        };
        members.retain(|m| m.cond == cond);
    }

    if members.is_empty() {
        println!("inspect_face face=({}, {}) no literal CCZ members", a, b);
        return;
    }
    members.sort_by_key(|m| m.op_index);

    let mut odd = BTreeSet::new();
    let mut thirds = BTreeMap::<u64, usize>::new();
    let mut conds = BTreeMap::<CondSig, usize>::new();
    for member in &members {
        if !odd.insert(member.third) {
            odd.remove(&member.third);
        }
        *thirds.entry(member.third).or_insert(0) += 1;
        *conds.entry(member.cond.clone()).or_insert(0) += 1;
    }

    let first = members.first().unwrap().op_index;
    let last = members.last().unwrap().op_index;
    println!(
        "inspect_face face=({}, {}) members={} span={}..{} width={} unique_thirds={} xor_terms={} cond_contexts={}",
        a,
        b,
        members.len(),
        first,
        last,
        last.saturating_sub(first),
        thirds.len(),
        odd.len(),
        conds.len(),
    );
    for (cond, count) in conds.iter().take(print_limit) {
        println!("inspect_cond cond={} members={}", format_cond(cond), count);
    }

    print_member_sample("member_head", &members, ops, sites, print_limit);
    let tail_start = members.len().saturating_sub(print_limit);
    print_member_sample(
        "member_tail",
        &members[tail_start..],
        ops,
        sites,
        print_limit,
    );

    let mut third_counts: Vec<_> = thirds.into_iter().collect();
    third_counts.sort_by(|left, right| right.1.cmp(&left.1).then_with(|| left.0.cmp(&right.0)));
    for (third, count) in third_counts.iter().take(print_limit) {
        println!("third_count q{} count={}", third, count);
    }
    for third in odd.iter().take(print_limit) {
        println!("xor_term q{}", third);
    }

    inspect_barriers(ops, sites, a, b, &members, &odd, print_limit);
}

fn face_members(ops: &[Op], a: u64, b: u64) -> Vec<FaceMember> {
    let mut out = Vec::new();
    let mut stack = Vec::new();
    for (idx, op) in ops.iter().enumerate() {
        match op.kind {
            OperationType::PushCondition => {
                stack.push(op.c_condition.0);
                continue;
            }
            OperationType::PopCondition => {
                stack.pop();
                continue;
            }
            _ => {}
        }
        if op.kind != OperationType::CCZ {
            continue;
        }
        let q = [op.q_control1.0, op.q_control2.0, op.q_target.0];
        for &(i, j, k) in &[(0usize, 1usize, 2usize), (0, 2, 1), (1, 2, 0)] {
            let (fa, fb) = if q[i] <= q[j] {
                (q[i], q[j])
            } else {
                (q[j], q[i])
            };
            if fa == a && fb == b {
                out.push(FaceMember {
                    op_index: idx,
                    third: q[k],
                    cond: cond_sig(&stack, op.c_condition),
                });
            }
        }
    }
    out
}

fn print_member_sample(
    label: &str,
    members: &[FaceMember],
    ops: &[Op],
    sites: &[point_add::OpSite],
    limit: usize,
) {
    for member in members.iter().take(limit) {
        let op = ops[member.op_index];
        println!(
            "{} op={} third=q{} cond={} gate={} site={}",
            label,
            member.op_index,
            member.third,
            format_cond(&member.cond),
            format_op(op),
            format_site(sites, member.op_index),
        );
    }
}

fn inspect_barriers(
    ops: &[Op],
    sites: &[point_add::OpSite],
    a: u64,
    b: u64,
    members: &[FaceMember],
    odd_terms: &BTreeSet<u64>,
    print_limit: usize,
) {
    let first = members.first().unwrap().op_index;
    let last = members.last().unwrap().op_index;
    let mut support = odd_terms.clone();
    support.insert(a);
    support.insert(b);

    let mut member_set = BTreeSet::new();
    for member in members {
        member_set.insert(member.op_index);
    }

    let mut face_write_barriers = Vec::new();
    let mut support_write_barriers = Vec::new();
    let mut condition_transitions = Vec::new();
    let mut diagonal_touches = 0usize;
    let mut non_diag_reads = 0usize;
    let mut segments: Vec<Vec<usize>> = vec![Vec::new()];

    for idx in first..=last {
        let op = ops[idx];
        if member_set.contains(&idx) {
            segments.last_mut().unwrap().push(idx);
            continue;
        }
        if matches!(
            op.kind,
            OperationType::PushCondition | OperationType::PopCondition
        ) {
            condition_transitions.push(idx);
            if !segments.last().unwrap().is_empty() {
                segments.push(Vec::new());
            }
            continue;
        }
        let touches_face = op_touches_any(op, &[a, b]);
        let writes_face = op_writes_any(op, &[a, b]);
        let touches_support = op_touches_set(op, &support);
        let writes_support = op_writes_set(op, &support);
        if is_diagonal_phase_op(op.kind) && touches_support {
            diagonal_touches += 1;
        } else if touches_support && !writes_support {
            non_diag_reads += 1;
        }
        if writes_face || (touches_face && !is_diagonal_phase_op(op.kind) && may_not_commute(op)) {
            face_write_barriers.push(idx);
        }
        if writes_support {
            support_write_barriers.push(idx);
            if !segments.last().unwrap().is_empty() {
                segments.push(Vec::new());
            }
        }
    }
    segments.retain(|segment| !segment.is_empty());
    let mut segment_saves = Vec::new();
    for segment in &segments {
        let mut odd = BTreeSet::new();
        for &idx in segment {
            if let Some(member) = members.iter().find(|m| m.op_index == idx) {
                if !odd.insert(member.third) {
                    odd.remove(&member.third);
                }
            }
        }
        let save = segment.len().saturating_sub(usize::from(!odd.is_empty()));
        if save > 0 {
            segment_saves.push((
                save,
                segment.len(),
                segment[0],
                *segment.last().unwrap(),
                odd.len(),
            ));
        }
    }
    segment_saves.sort_by(|left, right| {
        right
            .0
            .cmp(&left.0)
            .then_with(|| right.1.cmp(&left.1))
            .then_with(|| left.2.cmp(&right.2))
    });

    println!(
        "barrier_summary face_write_or_noncommuting={} support_write={} condition_transitions={} diagonal_support_touches={} non_diag_support_reads={} conservative_segment_count={} conservative_segment_save_sum={}",
        face_write_barriers.len(),
        support_write_barriers.len(),
        condition_transitions.len(),
        diagonal_touches,
        non_diag_reads,
        segments.len(),
        segment_saves.iter().map(|s| s.0).sum::<usize>(),
    );
    print_barriers(
        "face_barrier",
        &face_write_barriers,
        ops,
        sites,
        print_limit,
    );
    print_barriers(
        "support_write_barrier",
        &support_write_barriers,
        ops,
        sites,
        print_limit,
    );
    print_barriers(
        "condition_transition",
        &condition_transitions,
        ops,
        sites,
        print_limit,
    );
    for (save, len, start, end, odd_count) in segment_saves.iter().take(print_limit) {
        println!(
            "conservative_segment save={} members={} span={}..{} xor_terms={}",
            save, len, start, end, odd_count
        );
    }
}

fn print_barriers(
    label: &str,
    indices: &[usize],
    ops: &[Op],
    sites: &[point_add::OpSite],
    limit: usize,
) {
    for &idx in indices.iter().take(limit) {
        println!(
            "{} op={} gate={} site={}",
            label,
            idx,
            format_op(ops[idx]),
            format_site(sites, idx),
        );
    }
}

fn op_touches_set(op: Op, support: &BTreeSet<u64>) -> bool {
    op_qubits(op).iter().any(|q| support.contains(q))
}

fn op_writes_set(op: Op, support: &BTreeSet<u64>) -> bool {
    op_write_qubits(op).iter().any(|q| support.contains(q))
}

fn op_touches_any(op: Op, qs: &[u64]) -> bool {
    let touched = op_qubits(op);
    qs.iter().any(|q| touched.contains(q))
}

fn op_writes_any(op: Op, qs: &[u64]) -> bool {
    let written = op_write_qubits(op);
    qs.iter().any(|q| written.contains(q))
}

fn op_qubits(op: Op) -> Vec<u64> {
    [op.q_control2, op.q_control1, op.q_target]
        .into_iter()
        .filter(|q| *q != NO_QUBIT)
        .map(|q| q.0)
        .collect()
}

fn op_write_qubits(op: Op) -> Vec<u64> {
    match op.kind {
        OperationType::X
        | OperationType::CX
        | OperationType::CCX
        | OperationType::R
        | OperationType::Hmr => vec![op.q_target.0],
        OperationType::Swap => vec![op.q_control1.0, op.q_target.0],
        _ => Vec::new(),
    }
}

fn may_not_commute(op: Op) -> bool {
    !matches!(
        op.kind,
        OperationType::Neg
            | OperationType::Z
            | OperationType::CZ
            | OperationType::CCZ
            | OperationType::Register
            | OperationType::AppendToRegister
            | OperationType::BitInvert
            | OperationType::BitStore0
            | OperationType::BitStore1
            | OperationType::DebugPrint
    )
}

fn format_op(op: Op) -> String {
    let mut parts = vec![format!("{:?}", op.kind)];
    for q in [op.q_control2, op.q_control1, op.q_target] {
        if q != NO_QUBIT {
            parts.push(format!("q{}", q.0));
        }
    }
    if op.c_target != NO_BIT {
        parts.push(format!("b{}", op.c_target.0));
    }
    if op.c_condition != NO_BIT {
        parts.push(format!("if=b{}", op.c_condition.0));
    }
    parts.join("/")
}

fn format_site(sites: &[point_add::OpSite], idx: usize) -> String {
    if sites.len() <= idx {
        return "unavailable".to_owned();
    }
    let (file, line, context) = sites[idx];
    if context == 0 {
        format!("{file}:{line}")
    } else {
        format!("{file}:{line}#0x{context:08x}")
    }
}

#[derive(Clone, Debug)]
struct ConditionedSegment {
    face_a: u64,
    face_b: u64,
    save: usize,
    members: usize,
    first: usize,
    last: usize,
    terms: usize,
}

fn conditioned_parity_candidates(ops: &[Op]) -> Vec<ConditionedSegment> {
    let writes_by_qubit = writes_by_qubit(ops);
    let writes_by_bit = writes_by_bit(ops);
    let mut buckets: BTreeMap<(u64, u64), Vec<FaceMember>> = BTreeMap::new();
    let mut stack = Vec::new();

    for (idx, op) in ops.iter().enumerate() {
        match op.kind {
            OperationType::PushCondition => {
                stack.push(op.c_condition.0);
                continue;
            }
            OperationType::PopCondition => {
                stack.pop();
                continue;
            }
            _ => {}
        }
        if op.kind != OperationType::CCZ {
            continue;
        }
        let cond = cond_sig(&stack, op.c_condition);
        let q = [op.q_control1.0, op.q_control2.0, op.q_target.0];
        for &(i, j, k) in &[(0usize, 1usize, 2usize), (0, 2, 1), (1, 2, 0)] {
            let (a, b) = if q[i] <= q[j] {
                (q[i], q[j])
            } else {
                (q[j], q[i])
            };
            buckets.entry((a, b)).or_default().push(FaceMember {
                op_index: idx,
                third: q[k],
                cond: cond.clone(),
            });
        }
    }

    let mut out = Vec::new();
    for ((a, b), mut members) in buckets {
        if members.len() <= 1 {
            continue;
        }
        members.sort_by_key(|m| m.op_index);
        out.extend(conditioned_segments_for_bucket(
            a,
            b,
            &members,
            &writes_by_qubit,
            &writes_by_bit,
        ));
    }
    out.sort_by(|left, right| {
        right
            .save
            .cmp(&left.save)
            .then_with(|| right.members.cmp(&left.members))
            .then_with(|| left.first.cmp(&right.first))
    });
    out
}

fn conditioned_segments_for_bucket(
    a: u64,
    b: u64,
    members: &[FaceMember],
    writes_by_qubit: &[Vec<usize>],
    writes_by_bit: &[Vec<usize>],
) -> Vec<ConditionedSegment> {
    let mut out = Vec::new();
    let mut segment: Vec<&FaceMember> = Vec::new();
    let mut q_support = BTreeSet::new();
    let mut bit_support = BTreeSet::new();
    q_support.insert(a);
    q_support.insert(b);

    for member in members {
        if let Some(prev) = segment.last() {
            if has_any_write_between(&q_support, writes_by_qubit, prev.op_index, member.op_index)
                || has_any_write_between(
                    &bit_support,
                    writes_by_bit,
                    prev.op_index,
                    member.op_index,
                )
            {
                push_conditioned_segment(a, b, &segment, &mut out);
                segment.clear();
                q_support.clear();
                bit_support.clear();
                q_support.insert(a);
                q_support.insert(b);
            }
        }
        q_support.insert(member.third);
        for &bit in &member.cond.0 {
            bit_support.insert(bit);
        }
        segment.push(member);
    }
    push_conditioned_segment(a, b, &segment, &mut out);
    out
}

fn push_conditioned_segment(
    a: u64,
    b: u64,
    segment: &[&FaceMember],
    out: &mut Vec<ConditionedSegment>,
) {
    if segment.is_empty() {
        return;
    }
    let mut terms = BTreeSet::new();
    for member in segment {
        let mut term = member.cond.0.clone();
        term.push(member.third);
        terms.insert(term);
    }
    let replacement = usize::from(!terms.is_empty());
    let save = segment.len().saturating_sub(replacement);
    if save == 0 {
        return;
    }
    out.push(ConditionedSegment {
        face_a: a,
        face_b: b,
        save,
        members: segment.len(),
        first: segment.first().unwrap().op_index,
        last: segment.last().unwrap().op_index,
        terms: terms.len(),
    });
}

fn writes_by_qubit(ops: &[Op]) -> Vec<Vec<usize>> {
    let max_q = ops.iter().flat_map(|op| op_qubits(*op)).max().unwrap_or(0) as usize;
    let mut out = vec![Vec::new(); max_q + 1];
    for (idx, op) in ops.iter().enumerate() {
        for q in op_write_qubits(*op) {
            if let Some(slot) = out.get_mut(q as usize) {
                slot.push(idx);
            }
        }
    }
    out
}

fn writes_by_bit(ops: &[Op]) -> Vec<Vec<usize>> {
    let max_b = ops
        .iter()
        .filter_map(|op| (op.c_target != NO_BIT).then_some(op.c_target.0))
        .max()
        .unwrap_or(0) as usize;
    let mut out = vec![Vec::new(); max_b + 1];
    for (idx, op) in ops.iter().enumerate() {
        if op.c_target != NO_BIT {
            out[op.c_target.0 as usize].push(idx);
        }
    }
    out
}

fn has_any_write_between(
    support: &BTreeSet<u64>,
    writes_by_id: &[Vec<usize>],
    lo: usize,
    hi: usize,
) -> bool {
    support.iter().any(|&id| {
        writes_by_id
            .get(id as usize)
            .is_some_and(|writes| has_write_between(writes, lo, hi))
    })
}

fn has_write_between(writes: &[usize], lo: usize, hi: usize) -> bool {
    let pos = writes.partition_point(|&idx| idx <= lo);
    writes.get(pos).is_some_and(|&idx| idx < hi)
}

fn print_conditioned_candidates(candidates: &[ConditionedSegment], limit: usize) {
    let total: usize = candidates.iter().map(|c| c.save).sum();
    println!(
        "conditioned_parity_write_safe_save={} segment_count={} showing_top={}",
        total,
        candidates.len(),
        limit.min(candidates.len())
    );
    for c in candidates.iter().take(limit) {
        println!(
            "conditioned_parity save={} members={} terms={} face=({}, {}) span={}..{}",
            c.save, c.members, c.terms, c.face_a, c.face_b, c.first, c.last
        );
    }
}

fn print_ccz_site_counts(ops: &[Op], sites: &[point_add::OpSite], limit: usize) {
    if sites.len() != ops.len() {
        println!(
            "ccz_site_counts unavailable sites_len={} ops_len={}",
            sites.len(),
            ops.len()
        );
        return;
    }
    let mut counts = BTreeMap::<String, usize>::new();
    for (idx, op) in ops.iter().enumerate() {
        if op.kind != OperationType::CCZ {
            continue;
        }
        *counts.entry(format_site(sites, idx)).or_insert(0) += 1;
    }
    let mut rows: Vec<_> = counts.into_iter().collect();
    rows.sort_by(|left, right| right.1.cmp(&left.1).then_with(|| left.0.cmp(&right.0)));
    println!(
        "ccz_site_count_sites={} showing_top={}",
        rows.len(),
        limit.min(rows.len())
    );
    for (site, count) in rows.into_iter().take(limit) {
        println!("ccz_site count={} site={}", count, site);
    }
}

fn print_ordered_site_groups(
    ops: &[Op],
    sites: &[point_add::OpSite],
    site_filter: &str,
    limit: usize,
) {
    if sites.len() != ops.len() {
        println!(
            "ordered_site_groups unavailable sites_len={} ops_len={}",
            sites.len(),
            ops.len()
        );
        return;
    }
    let mut groups =
        BTreeMap::<(u64, u64), (usize, BTreeSet<u64>, BTreeSet<String>, usize, usize)>::new();
    for (idx, op) in ops.iter().enumerate() {
        if op.kind != OperationType::CCZ {
            continue;
        }
        let site = format_site(sites, idx);
        if !site.contains(site_filter) {
            continue;
        }
        let entry = groups
            .entry((op.q_control2.0, op.q_control1.0))
            .or_insert_with(|| (0, BTreeSet::new(), BTreeSet::new(), idx, idx));
        entry.0 += 1;
        entry.1.insert(op.q_target.0);
        if op.c_condition != NO_BIT {
            entry.2.insert(op.c_condition.0.to_string());
        }
        entry.3 = entry.3.min(idx);
        entry.4 = entry.4.max(idx);
    }
    let mut rows: Vec<_> = groups.into_iter().collect();
    rows.sort_by(|left, right| {
        right
            .1
             .0
            .cmp(&left.1 .0)
            .then_with(|| left.0.cmp(&right.0))
    });
    let total: usize = rows.iter().map(|(_, row)| row.0).sum();
    let theoretical_save: usize = rows.iter().map(|(_, row)| row.0.saturating_sub(1)).sum();
    println!(
        "ordered_site_groups filter={} total_ccz={} groups={} theoretical_same_ordered_face_save={} showing_top={}",
        site_filter,
        total,
        rows.len(),
        theoretical_save,
        limit.min(rows.len())
    );
    for ((ctrl, other), (count, targets, conds, first, last)) in rows.into_iter().take(limit) {
        println!(
            "ordered_site_group count={} save_ceiling={} ordered_face=(q{}, q{}) unique_targets={} direct_conditions={} span={}..{}",
            count,
            count.saturating_sub(1),
            ctrl,
            other,
            targets.len(),
            conds.len(),
            first,
            last
        );
    }
}

#[derive(Clone, Debug)]
struct FaceFactorSegment {
    face_a: u64,
    face_b: u64,
    save: usize,
    members: usize,
    first: usize,
    last: usize,
    start_insert: Option<usize>,
    end_insert: Option<usize>,
    conditions: usize,
    targets: usize,
}

fn print_face_factor_segments(
    ops: &[Op],
    sites: &[point_add::OpSite],
    site_filter: &str,
    limit: usize,
) {
    if sites.len() != ops.len() {
        println!(
            "face_factor_segments unavailable sites_len={} ops_len={}",
            sites.len(),
            ops.len()
        );
        return;
    }

    let writes_by_qubit = writes_by_qubit(ops);
    let depth_before = condition_depth_before(ops);
    let mut stack = Vec::new();
    let mut groups = BTreeMap::<(u64, u64), Vec<FaceMember>>::new();
    let mut total_site_ccz = 0usize;

    for (idx, op) in ops.iter().enumerate() {
        match op.kind {
            OperationType::PushCondition => {
                stack.push(op.c_condition.0);
                continue;
            }
            OperationType::PopCondition => {
                stack.pop();
                continue;
            }
            _ => {}
        }
        if op.kind != OperationType::CCZ {
            continue;
        }
        let site = format_site(sites, idx);
        if !site.contains(site_filter) {
            continue;
        }
        total_site_ccz += 1;
        groups
            .entry((op.q_control2.0, op.q_control1.0))
            .or_default()
            .push(FaceMember {
                op_index: idx,
                third: op.q_target.0,
                cond: cond_sig(&stack, op.c_condition),
            });
    }

    let mut segments = Vec::new();
    for ((a, b), mut members) in groups {
        if members.len() <= 1 {
            continue;
        }
        members.sort_by_key(|m| m.op_index);
        collect_face_factor_segments(
            a,
            b,
            &members,
            &writes_by_qubit,
            &depth_before,
            ops.len(),
            &mut segments,
        );
    }
    segments.sort_by(|left, right| {
        right
            .save
            .cmp(&left.save)
            .then_with(|| right.members.cmp(&left.members))
            .then_with(|| left.first.cmp(&right.first))
    });

    let raw_save: usize = segments.iter().map(|s| s.save).sum();
    let neutral_save: usize = segments
        .iter()
        .filter(|s| s.start_insert.is_some() && s.end_insert.is_some())
        .map(|s| s.save)
        .sum();
    let best_extra_qubit_score_delta = segments
        .first()
        .map(|s| (1_364_229.770_f64 - s.save as f64) * 1153.0 - 1_364_229.770_f64 * 1152.0)
        .unwrap_or(0.0);

    println!(
        "face_factor_segments filter={} total_site_ccz={} segment_count={} raw_save={} stack_neutral_save={} best_extra_qubit_score_delta={:.3} showing_top={}",
        site_filter,
        total_site_ccz,
        segments.len(),
        raw_save,
        neutral_save,
        best_extra_qubit_score_delta,
        limit.min(segments.len())
    );
    for segment in segments.iter().take(limit) {
        println!(
            "face_factor save={} members={} face=(q{}, q{}) targets={} conds={} span={}..{} start_insert={} end_insert={}",
            segment.save,
            segment.members,
            segment.face_a,
            segment.face_b,
            segment.targets,
            segment.conditions,
            segment.first,
            segment.last,
            format_insert(segment.start_insert),
            format_insert(segment.end_insert),
        );
    }
}

fn collect_face_factor_segments(
    a: u64,
    b: u64,
    members: &[FaceMember],
    writes_by_qubit: &[Vec<usize>],
    depth_before: &[usize],
    ops_len: usize,
    out: &mut Vec<FaceFactorSegment>,
) {
    let mut segment: Vec<&FaceMember> = Vec::new();
    let mut face = BTreeSet::new();
    face.insert(a);
    face.insert(b);

    for member in members {
        if let Some(prev) = segment.last() {
            if has_any_write_between(&face, writes_by_qubit, prev.op_index, member.op_index) {
                push_face_factor_segment(
                    a,
                    b,
                    &segment,
                    writes_by_qubit,
                    depth_before,
                    ops_len,
                    out,
                );
                segment.clear();
            }
        }
        segment.push(member);
    }
    push_face_factor_segment(a, b, &segment, writes_by_qubit, depth_before, ops_len, out);
}

fn push_face_factor_segment(
    a: u64,
    b: u64,
    segment: &[&FaceMember],
    writes_by_qubit: &[Vec<usize>],
    depth_before: &[usize],
    ops_len: usize,
    out: &mut Vec<FaceFactorSegment>,
) {
    if segment.len() <= 1 {
        return;
    }
    let first = segment.first().unwrap().op_index;
    let last = segment.last().unwrap().op_index;
    let prev_write = previous_face_write(writes_by_qubit, a, b, first);
    let next_write = next_face_write(writes_by_qubit, a, b, last);
    let start_lo = prev_write.map_or(0, |idx| idx.saturating_add(1));
    let end_hi = next_write.unwrap_or(ops_len);
    let start_insert = latest_zero_depth_boundary(depth_before, start_lo, first);
    let end_insert = earliest_zero_depth_boundary(depth_before, last.saturating_add(1), end_hi);

    let mut conditions = BTreeSet::new();
    let mut targets = BTreeSet::new();
    for member in segment {
        conditions.insert(member.cond.clone());
        targets.insert(member.third);
    }
    out.push(FaceFactorSegment {
        face_a: a,
        face_b: b,
        save: segment.len() - 1,
        members: segment.len(),
        first,
        last,
        start_insert,
        end_insert,
        conditions: conditions.len(),
        targets: targets.len(),
    });
}

fn condition_depth_before(ops: &[Op]) -> Vec<usize> {
    let mut out = Vec::with_capacity(ops.len());
    let mut depth = 0usize;
    for op in ops {
        out.push(depth);
        match op.kind {
            OperationType::PushCondition => depth += 1,
            OperationType::PopCondition => depth = depth.saturating_sub(1),
            _ => {}
        }
    }
    out
}

fn previous_face_write(
    writes_by_qubit: &[Vec<usize>],
    a: u64,
    b: u64,
    before: usize,
) -> Option<usize> {
    [a, b]
        .into_iter()
        .filter_map(|q| {
            let writes = writes_by_qubit.get(q as usize)?;
            let pos = writes.partition_point(|&idx| idx < before);
            pos.checked_sub(1).map(|i| writes[i])
        })
        .max()
}

fn next_face_write(writes_by_qubit: &[Vec<usize>], a: u64, b: u64, after: usize) -> Option<usize> {
    [a, b]
        .into_iter()
        .filter_map(|q| {
            let writes = writes_by_qubit.get(q as usize)?;
            let pos = writes.partition_point(|&idx| idx <= after);
            writes.get(pos).copied()
        })
        .min()
}

fn latest_zero_depth_boundary(depth_before: &[usize], lo: usize, hi: usize) -> Option<usize> {
    if lo > hi || lo >= depth_before.len() {
        return None;
    }
    let hi = hi.min(depth_before.len().saturating_sub(1));
    (lo..=hi).rev().find(|&idx| depth_before[idx] == 0)
}

fn earliest_zero_depth_boundary(depth_before: &[usize], lo: usize, hi: usize) -> Option<usize> {
    if lo > hi {
        return None;
    }
    if lo == depth_before.len() {
        return Some(lo);
    }
    let hi = hi.min(depth_before.len().saturating_sub(1));
    (lo..=hi).find(|&idx| depth_before[idx] == 0)
}

fn format_insert(insert: Option<usize>) -> String {
    insert.map_or_else(|| "none".to_owned(), |idx| idx.to_string())
}
