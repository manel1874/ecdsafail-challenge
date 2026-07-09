use std::collections::BTreeMap;

use quantum_ecc::circuit::{Op, OperationType, NO_QUBIT};
use quantum_ecc::point_add;

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
struct CcxKey {
    a: u64,
    b: u64,
    target: u64,
}

#[derive(Clone, Debug)]
struct MbuPair {
    compute: usize,
    clear: usize,
    key: CcxKey,
}

fn main() {
    std::env::set_var("TRACE_OP_SITES", "1");
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
    println!(
        "and_island_report ops={} sites={} ccx={} ccz={}",
        ops.len(),
        sites.len(),
        ccx,
        ccz
    );

    let limit = std::env::var("AND_REPORT_LIMIT")
        .ok()
        .and_then(|value| value.parse::<usize>().ok())
        .unwrap_or(40);

    print_ccx_site_counts(&ops, &sites, limit);
    print_ccx_context_prefix_counts(&ops, &sites);

    let candidates = conservative_mbu_pairs(&ops);
    println!(
        "conservative_mbu_pair_candidates={} toffoli_save_ceiling={}",
        candidates.len(),
        candidates.len()
    );
    print_mbu_pair_sites(&candidates, &ops, &sites, limit);
    print_mbu_pair_samples(&candidates, &ops, &sites, limit.min(24));
}

fn print_ccx_site_counts(ops: &[Op], sites: &[point_add::OpSite], limit: usize) {
    if ops.len() != sites.len() {
        println!(
            "ccx_site_counts unavailable sites_len={} ops_len={}",
            sites.len(),
            ops.len()
        );
        return;
    }
    let mut full = BTreeMap::<String, usize>::new();
    let mut base = BTreeMap::<String, usize>::new();
    for (idx, op) in ops.iter().enumerate() {
        if op.kind != OperationType::CCX {
            continue;
        }
        *full.entry(format_site(sites, idx, true)).or_insert(0) += 1;
        *base.entry(format_site(sites, idx, false)).or_insert(0) += 1;
    }
    print_counts("ccx_site", full, limit);
    print_counts("ccx_base_site", base, limit);
}

fn print_ccx_context_prefix_counts(ops: &[Op], sites: &[point_add::OpSite]) {
    if ops.len() != sites.len() {
        return;
    }
    let mut counts = BTreeMap::<u32, usize>::new();
    for (idx, op) in ops.iter().enumerate() {
        if op.kind != OperationType::CCX {
            continue;
        }
        let (_, _, context) = sites[idx];
        *counts.entry(context >> 24).or_insert(0) += 1;
    }
    let mut rows: Vec<_> = counts.into_iter().collect();
    rows.sort_by(|left, right| right.1.cmp(&left.1).then_with(|| left.0.cmp(&right.0)));
    println!("ccx_context_prefixes={}", rows.len());
    for (prefix, count) in rows {
        println!("ccx_context_prefix prefix=0x{prefix:02x} count={count}");
    }
}

fn print_counts(label: &str, counts: BTreeMap<String, usize>, limit: usize) {
    let mut rows: Vec<_> = counts.into_iter().collect();
    rows.sort_by(|left, right| right.1.cmp(&left.1).then_with(|| left.0.cmp(&right.0)));
    println!(
        "{}_count={} showing_top={}",
        label,
        rows.len(),
        limit.min(rows.len())
    );
    for (site, count) in rows.into_iter().take(limit) {
        println!("{label} count={count} site={site}");
    }
}

fn conservative_mbu_pairs(ops: &[Op]) -> Vec<MbuPair> {
    let writes = writes_by_qubit(ops);
    let touches = touches_by_qubit(ops);
    let mut out = Vec::new();
    for (idx, op) in ops.iter().enumerate() {
        if op.kind != OperationType::CCX {
            continue;
        }
        let Some(key) = ccx_key(*op) else {
            continue;
        };
        if key.target < 512 {
            continue;
        }
        if !previous_write_is_zeroing(ops, &writes, key.target, idx) {
            continue;
        }
        let Some(next) = next_write(&writes, key.target, idx) else {
            continue;
        };
        if ops[next].kind != OperationType::CCX || ccx_key(ops[next]) != Some(key) {
            continue;
        }
        if has_touch_between(&touches, key.target, idx, next) {
            continue;
        }
        if has_write_between(&writes, key.a, idx, next)
            || has_write_between(&writes, key.b, idx, next)
        {
            continue;
        }
        out.push(MbuPair {
            compute: idx,
            clear: next,
            key,
        });
    }
    out
}

fn previous_write_is_zeroing(ops: &[Op], writes: &[Vec<usize>], q: u64, before: usize) -> bool {
    let Some(row) = writes.get(q as usize) else {
        return true;
    };
    let pos = row.partition_point(|&idx| idx < before);
    let Some(prev_pos) = pos.checked_sub(1) else {
        return true;
    };
    matches!(
        ops[row[prev_pos]].kind,
        OperationType::R | OperationType::Hmr
    )
}

fn next_write(writes: &[Vec<usize>], q: u64, after: usize) -> Option<usize> {
    let row = writes.get(q as usize)?;
    let pos = row.partition_point(|&idx| idx <= after);
    row.get(pos).copied()
}

fn has_write_between(writes: &[Vec<usize>], q: u64, lo: usize, hi: usize) -> bool {
    writes.get(q as usize).is_some_and(|row| {
        let pos = row.partition_point(|&idx| idx <= lo);
        row.get(pos).is_some_and(|&idx| idx < hi)
    })
}

fn has_touch_between(touches: &[Vec<usize>], q: u64, lo: usize, hi: usize) -> bool {
    touches.get(q as usize).is_some_and(|row| {
        let pos = row.partition_point(|&idx| idx <= lo);
        row.get(pos).is_some_and(|&idx| idx < hi)
    })
}

fn writes_by_qubit(ops: &[Op]) -> Vec<Vec<usize>> {
    let max_q = ops.iter().flat_map(op_qubits).max().unwrap_or(0) as usize;
    let mut out = vec![Vec::new(); max_q + 1];
    for (idx, op) in ops.iter().enumerate() {
        for q in op_write_qubits(*op) {
            if let Some(row) = out.get_mut(q as usize) {
                row.push(idx);
            }
        }
    }
    out
}

fn touches_by_qubit(ops: &[Op]) -> Vec<Vec<usize>> {
    let max_q = ops.iter().flat_map(op_qubits).max().unwrap_or(0) as usize;
    let mut out = vec![Vec::new(); max_q + 1];
    for (idx, op) in ops.iter().enumerate() {
        for q in op_qubits(op) {
            if let Some(row) = out.get_mut(q as usize) {
                row.push(idx);
            }
        }
    }
    out
}

fn op_qubits(op: &Op) -> Vec<u64> {
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

fn ccx_key(op: Op) -> Option<CcxKey> {
    if op.kind != OperationType::CCX {
        return None;
    }
    let a = op.q_control1.0.min(op.q_control2.0);
    let b = op.q_control1.0.max(op.q_control2.0);
    Some(CcxKey {
        a,
        b,
        target: op.q_target.0,
    })
}

fn print_mbu_pair_sites(pairs: &[MbuPair], _ops: &[Op], sites: &[point_add::OpSite], limit: usize) {
    let mut clear_sites = BTreeMap::<String, usize>::new();
    let mut compute_sites = BTreeMap::<String, usize>::new();
    for pair in pairs {
        *clear_sites
            .entry(format_site(sites, pair.clear, false))
            .or_insert(0) += 1;
        *compute_sites
            .entry(format_site(sites, pair.compute, false))
            .or_insert(0) += 1;
    }
    print_counts("mbu_clear_base_site", clear_sites, limit);
    print_counts("mbu_compute_base_site", compute_sites, limit);
}

fn print_mbu_pair_samples(
    pairs: &[MbuPair],
    _ops: &[Op],
    sites: &[point_add::OpSite],
    limit: usize,
) {
    for pair in pairs.iter().take(limit) {
        println!(
            "mbu_pair compute={} clear={} span={} key=(q{},q{})->q{} compute_site={} clear_site={}",
            pair.compute,
            pair.clear,
            pair.clear.saturating_sub(pair.compute),
            pair.key.a,
            pair.key.b,
            pair.key.target,
            format_site(sites, pair.compute, true),
            format_site(sites, pair.clear, true),
        );
    }
}

fn format_site(sites: &[point_add::OpSite], idx: usize, include_context: bool) -> String {
    if sites.len() <= idx {
        return "unavailable".to_owned();
    }
    let (file, line, context) = sites[idx];
    if include_context && context != 0 {
        format!("{file}:{line}#0x{context:08x}")
    } else {
        format!("{file}:{line}")
    }
}
