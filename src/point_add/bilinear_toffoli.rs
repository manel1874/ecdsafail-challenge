//! Opt-in bilinear Toffoli duplicate-and-destroy tooling.
//!
//! This is the Toffoli-count analogue of the parity-table duplicate-and-destroy
//! move used by T-count optimizers: profile expensive nonlinear columns on the
//! verifier-shaped support, find GF(2) relations among their bilinear fire masks,
//! and remove relations that cannot be observed before they cancel.

use std::collections::{BTreeSet, HashMap};

use alloy_primitives::U256;
use sha3::{
    digest::{ExtendableOutput, Update, XofReader},
    Shake256,
};

use crate::circuit::{analyze_ops, Op, OperationType, QubitId, NO_BIT};
use crate::sim::Simulator;
use crate::weierstrass_elliptic_curve::WeierstrassEllipticCurve;

const DEFAULT_PROFILE_SHOTS: usize = 9024;
const DEFAULT_PROFILE_TOP: usize = 80;
const NO_SLOT: u32 = u32::MAX;

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
struct SpanKey {
    target: u64,
    epoch: u32,
}

#[derive(Clone, Debug)]
struct CcxSig {
    op_index: usize,
    target: QubitId,
    control1: QubitId,
    control2: QubitId,
    span: SpanKey,
    executed: u64,
    fired: u64,
    h1: u64,
    h2: u64,
    l1: u64,
    l2: u64,
}

#[derive(Clone, Debug)]
struct PairCandidate {
    i: usize,
    j: usize,
    span: SpanKey,
    exec_sum: u64,
    fired: u64,
    distance: usize,
}

#[derive(Clone, Debug)]
struct TripleCandidate {
    i: usize,
    j: usize,
    k: usize,
    span: SpanKey,
    exec_sum: u64,
    distance: usize,
}

#[derive(Default)]
struct DropSpec {
    singles: BTreeSet<usize>,
    pairs: Vec<(usize, usize)>,
    groups: Vec<Vec<usize>>,
}

pub fn run(ops: Vec<Op>) -> Vec<Op> {
    if std::env::var("BILINEAR_TOFFOLI_PROFILE").ok().as_deref() == Some("1") {
        profile_and_print(&ops);
    }

    let drops = parse_drop_spec();
    if drops.singles.is_empty() && drops.pairs.is_empty() && drops.groups.is_empty() {
        return ops;
    }

    apply_drops(ops, drops)
}

fn apply_drops(ops: Vec<Op>, drops: DropSpec) -> Vec<Op> {
    let mut drop_indices = drops.singles;
    for (a, b) in drops.pairs {
        drop_indices.insert(a);
        drop_indices.insert(b);
    }
    for group in drops.groups {
        for idx in group {
            drop_indices.insert(idx);
        }
    }

    let before_tof = toffoli_ops(&ops);
    let before_ops = ops.len();
    let mut dropped = 0usize;
    let mut out = Vec::with_capacity(ops.len().saturating_sub(drop_indices.len()));
    for (i, op) in ops.into_iter().enumerate() {
        if drop_indices.contains(&i) {
            if op.kind == OperationType::CCX {
                dropped += 1;
                continue;
            }
            eprintln!(
                "BILINEAR_TOFFOLI: requested drop index {i}, but op is {:?}; keeping it",
                op.kind
            );
        }
        out.push(op);
    }
    eprintln!(
        "BILINEAR_TOFFOLI: dropped_ccx={} ops {}->{} emitted_tof {}->{}",
        dropped,
        before_ops,
        out.len(),
        before_tof,
        toffoli_ops(&out),
    );
    out
}

fn profile_and_print(ops: &[Op]) {
    let shots = env_usize("BILINEAR_TOFFOLI_PROFILE_SHOTS", DEFAULT_PROFILE_SHOTS);
    let top = env_usize("BILINEAR_TOFFOLI_PROFILE_TOP", DEFAULT_PROFILE_TOP);
    let seeds = env_usize("BILINEAR_TOFFOLI_PROFILE_SEEDS", 1).max(1);

    let (total_qubits, num_bits, _num_regs, regs) = analyze_ops(ops.iter());
    if regs.len() != 4 || regs.iter().any(|r| r.len() != 256) {
        eprintln!(
            "BILINEAR_PROFILE: expected four 256-wide registers, got {}; skipping",
            regs.len()
        );
        return;
    }

    let spans = ccx_span_keys(ops, total_qubits as usize);
    let mut gate_slots = vec![NO_SLOT; ops.len()];
    let mut sigs = Vec::new();
    for (i, op) in ops.iter().enumerate() {
        if op.kind == OperationType::CCX {
            let Some(span) = spans[i] else {
                continue;
            };
            gate_slots[i] = sigs.len() as u32;
            sigs.push(CcxSig {
                op_index: i,
                target: op.q_target,
                control1: op.q_control1,
                control2: op.q_control2,
                span,
                executed: 0,
                fired: 0,
                h1: 0,
                h2: 0,
                l1: 0,
                l2: 0,
            });
        }
    }

    eprintln!(
        "BILINEAR_PROFILE: seeds={} shots_per_seed={} ops={} ccx={} qubits={} bits={}",
        seeds,
        shots,
        ops.len(),
        sigs.len(),
        total_qubits,
        num_bits,
    );

    let curve = secp256k1();
    let mut total_profiled_shots = 0usize;
    let mut block_id = 0u64;
    for seed in 0..seeds {
        let mut xof = if seed == 0 {
            official_xof(ops)
        } else {
            support_profile_xof(ops, seed as u64)
        };
        let mut targets = Vec::with_capacity(shots);
        let mut offsets = Vec::with_capacity(shots);
        for _ in 0..shots {
            let mut rb = [[0u8; 32]; 2];
            xof.read(&mut rb[0]);
            xof.read(&mut rb[1]);
            let k1 = U256::from_le_bytes(rb[0]);
            let k2 = U256::from_le_bytes(rb[1]);
            let t = curve.mul(curve.gx, curve.gy, k1);
            let o = curve.mul(curve.gx, curve.gy, k2);
            if t.0 == o.0 {
                continue;
            }
            if t.0.is_zero() && t.1.is_zero() {
                continue;
            }
            if o.0.is_zero() && o.1.is_zero() {
                continue;
            }
            targets.push(t);
            offsets.push(o);
        }

        total_profiled_shots += targets.len();
        let mut sim = Simulator::new(total_qubits as usize, num_bits as usize, &mut xof);
        const BATCH: usize = 64;
        let num_batches = (targets.len() + BATCH - 1) / BATCH;
        for batch in 0..num_batches {
            let bs = BATCH.min(targets.len() - batch * BATCH);
            let live_mask = if bs == 64 { u64::MAX } else { (1u64 << bs) - 1 };
            sim.clear_for_shot();
            for shot in 0..bs {
                let i = batch * BATCH + shot;
                sim.set_register(&regs[0], targets[i].0, shot);
                sim.set_register(&regs[1], targets[i].1, shot);
                sim.set_register(&regs[2], offsets[i].0, shot);
                sim.set_register(&regs[3], offsets[i].1, shot);
            }
            replay_with_signatures(&mut sim, ops, &gate_slots, &mut sigs, live_mask, block_id);
            block_id += 1;
        }
    }

    print_duplicate_candidates(&sigs, total_profiled_shots.max(1) as f64, top);
}

fn replay_with_signatures<R: XofReader>(
    sim: &mut Simulator<'_, R>,
    ops: &[Op],
    gate_slots: &[u32],
    sigs: &mut [CcxSig],
    live_mask: u64,
    block_id: u64,
) {
    let mut condition_stack = Vec::new();
    let mut current_base_condition = live_mask;

    for (i, op) in ops.iter().enumerate() {
        let mut cond = current_base_condition;
        if op.c_condition != NO_BIT {
            cond &= sim.bit(op.c_condition);
        }
        cond &= live_mask;

        let slot = gate_slots[i];
        if slot != NO_SLOT {
            let sig = &mut sigs[slot as usize];
            let fire = cond & sim.qubit(op.q_control1) & sim.qubit(op.q_control2);
            sig.executed += cond.count_ones() as u64;
            sig.fired += fire.count_ones() as u64;
            let salt = block_id.wrapping_mul(0x9e37_79b9_7f4a_7c15);
            sig.h1 ^= mix64(fire ^ salt);
            sig.h2 = sig
                .h2
                .wrapping_add(mix64(fire.wrapping_add(salt.rotate_left(17))));
            let (l1, l2) = linear_mix_mask(fire, block_id);
            sig.l1 ^= l1;
            sig.l2 ^= l2;
        }

        match op.kind {
            OperationType::CCX => {
                let v = cond & sim.qubit(op.q_control1) & sim.qubit(op.q_control2);
                *sim.qubit_mut(op.q_target) ^= v;
            }
            OperationType::CX => {
                let v = cond & sim.qubit(op.q_control1);
                *sim.qubit_mut(op.q_target) ^= v;
            }
            OperationType::Swap => {
                let mut q_c1 = sim.qubit(op.q_control1);
                let mut q_t = sim.qubit(op.q_target);
                q_c1 ^= q_t;
                q_t ^= cond & q_c1;
                q_c1 ^= q_t;
                *sim.qubit_mut(op.q_control1) = q_c1;
                *sim.qubit_mut(op.q_target) = q_t;
            }
            OperationType::X => {
                *sim.qubit_mut(op.q_target) ^= cond;
            }
            OperationType::CCZ => {
                let v = cond
                    & sim.qubit(op.q_target)
                    & sim.qubit(op.q_control1)
                    & sim.qubit(op.q_control2);
                sim.phase ^= v;
            }
            OperationType::CZ => {
                let v = cond & sim.qubit(op.q_target) & sim.qubit(op.q_control1);
                sim.phase ^= v;
            }
            OperationType::Z => {
                let v = cond & sim.qubit(op.q_target);
                sim.phase ^= v;
            }
            OperationType::Neg => {
                sim.phase ^= cond;
            }
            OperationType::Hmr => {
                let mut buf = [0u8; 8];
                sim.xof.read(&mut buf);
                let rng_val = u64::from_le_bytes(buf);
                *sim.bit_mut(op.c_target) &= !cond;
                *sim.bit_mut(op.c_target) ^= rng_val & cond;
                sim.phase ^= sim.qubit(op.q_target) & rng_val & cond;
                *sim.qubit_mut(op.q_target) &= !cond;
            }
            OperationType::R => {
                let mut buf = [0u8; 8];
                sim.xof.read(&mut buf);
                let rng_val = u64::from_le_bytes(buf);
                sim.phase ^= sim.qubit(op.q_target) & rng_val & cond;
                *sim.qubit_mut(op.q_target) &= !cond;
            }
            OperationType::BitInvert => {
                *sim.bit_mut(op.c_target) ^= cond;
            }
            OperationType::BitStore0 => {
                *sim.bit_mut(op.c_target) &= !cond;
            }
            OperationType::BitStore1 => {
                *sim.bit_mut(op.c_target) |= cond;
            }
            OperationType::AppendToRegister
            | OperationType::Register
            | OperationType::DebugPrint => {}
            OperationType::PushCondition => {
                condition_stack.push(current_base_condition);
                current_base_condition &= sim.bit(op.c_condition);
            }
            OperationType::PopCondition => {
                if let Some(val) = condition_stack.pop() {
                    current_base_condition = val;
                }
            }
        }
    }
}

fn print_duplicate_candidates(sigs: &[CcxSig], denom: f64, top: usize) {
    let mut singles: Vec<&CcxSig> = sigs
        .iter()
        .filter(|sig| sig.executed > 0 && sig.fired == 0)
        .collect();
    singles.sort_by(|a, b| {
        b.executed
            .cmp(&a.executed)
            .then_with(|| a.op_index.cmp(&b.op_index))
    });

    let mut groups: HashMap<(SpanKey, u64, u64, u64), Vec<&CcxSig>> = HashMap::new();
    let mut span_groups: HashMap<SpanKey, Vec<&CcxSig>> = HashMap::new();
    for sig in sigs {
        if sig.executed == 0 {
            continue;
        }
        span_groups.entry(sig.span).or_default().push(sig);
        groups
            .entry((sig.span, sig.fired, sig.h1, sig.h2))
            .or_default()
            .push(sig);
    }

    let mut pairs = Vec::new();
    for ((span, fired, _h1, _h2), group) in groups {
        if group.len() < 2 {
            continue;
        }
        for a in 0..group.len() {
            for b in a + 1..group.len() {
                let ga = group[a];
                let gb = group[b];
                pairs.push(PairCandidate {
                    i: ga.op_index,
                    j: gb.op_index,
                    span,
                    exec_sum: ga.executed + gb.executed,
                    fired,
                    distance: gb.op_index.abs_diff(ga.op_index),
                });
            }
        }
    }
    pairs.sort_by(|a, b| {
        b.exec_sum
            .cmp(&a.exec_sum)
            .then_with(|| a.distance.cmp(&b.distance))
            .then_with(|| a.i.cmp(&b.i))
            .then_with(|| a.j.cmp(&b.j))
    });

    let triple_group_max = env_usize("BILINEAR_TOFFOLI_TRIPLE_GROUP_MAX", 192);
    let mut triples = Vec::new();
    for (span, mut group) in span_groups {
        if group.len() < 3 || group.len() > triple_group_max {
            continue;
        }
        group.sort_by_key(|sig| sig.op_index);
        let mut by_linear: HashMap<(u64, u64), Vec<usize>> = HashMap::new();
        for (pos, sig) in group.iter().enumerate() {
            by_linear.entry((sig.l1, sig.l2)).or_default().push(pos);
        }
        for a in 0..group.len() {
            for b in a + 1..group.len() {
                let need = (group[a].l1 ^ group[b].l1, group[a].l2 ^ group[b].l2);
                let Some(cs) = by_linear.get(&need) else {
                    continue;
                };
                for &c in cs {
                    if c <= b {
                        continue;
                    }
                    let ga = group[a];
                    let gb = group[b];
                    let gc = group[c];
                    triples.push(TripleCandidate {
                        i: ga.op_index,
                        j: gb.op_index,
                        k: gc.op_index,
                        span,
                        exec_sum: ga.executed + gb.executed + gc.executed,
                        distance: gc.op_index - ga.op_index,
                    });
                }
            }
        }
    }
    triples.sort_by(|a, b| {
        b.exec_sum
            .cmp(&a.exec_sum)
            .then_with(|| a.distance.cmp(&b.distance))
            .then_with(|| a.i.cmp(&b.i))
            .then_with(|| a.j.cmp(&b.j))
            .then_with(|| a.k.cmp(&b.k))
    });

    eprintln!("=== bilinear zero-product single-CCX candidates ===");
    eprintln!(
        "{:<12} {:>12} {:>10} {:>8} {:>8} {:>8}",
        "index", "exec_sum", "exec/shot", "target", "c1", "c2"
    );
    for sig in singles.into_iter().take(top) {
        eprintln!(
            "{:<12} {:>12} {:>10.3} {:>8} {:>8} {:>8}",
            sig.op_index,
            sig.executed,
            sig.executed as f64 / denom,
            sig.target.0,
            sig.control1.0,
            sig.control2.0,
        );
    }

    eprintln!("=== bilinear duplicate-pair drop candidates ===");
    eprintln!(
        "{:<25} {:>12} {:>12} {:>10} {:>8} {:>8} {:>8}",
        "pair", "exec_sum", "fire_sum", "exec/shot", "target", "epoch", "gap"
    );
    for cand in pairs.into_iter().take(top) {
        eprintln!(
            "{:<25} {:>12} {:>12} {:>10.3} {:>8} {:>8} {:>8}",
            format!("{}:{}", cand.i, cand.j),
            cand.exec_sum,
            cand.fired,
            cand.exec_sum as f64 / denom,
            cand.span.target,
            cand.span.epoch,
            cand.distance,
        );
    }

    eprintln!("=== bilinear xor-zero triple drop candidates ===");
    eprintln!(
        "{:<37} {:>12} {:>10} {:>8} {:>8} {:>8}",
        "triple", "exec_sum", "exec/shot", "target", "epoch", "span"
    );
    for cand in triples.into_iter().take(top) {
        eprintln!(
            "{:<37} {:>12} {:>10.3} {:>8} {:>8} {:>8}",
            format!("{}:{}:{}", cand.i, cand.j, cand.k),
            cand.exec_sum,
            cand.exec_sum as f64 / denom,
            cand.span.target,
            cand.span.epoch,
            cand.distance,
        );
    }
}

fn ccx_span_keys(ops: &[Op], num_qubits: usize) -> Vec<Option<SpanKey>> {
    let mut epochs = vec![0u32; num_qubits];
    let mut out = vec![None; ops.len()];
    for (i, op) in ops.iter().enumerate() {
        match op.kind {
            OperationType::CCX => {
                mark_read(&mut epochs, op.q_control1);
                mark_read(&mut epochs, op.q_control2);
                out[i] = Some(SpanKey {
                    target: op.q_target.0,
                    epoch: epochs[op.q_target.0 as usize],
                });
            }
            OperationType::CX => {
                mark_read(&mut epochs, op.q_control1);
            }
            OperationType::Swap | OperationType::CZ => {
                mark_read(&mut epochs, op.q_control1);
                mark_read(&mut epochs, op.q_target);
            }
            OperationType::CCZ => {
                mark_read(&mut epochs, op.q_control1);
                mark_read(&mut epochs, op.q_control2);
                mark_read(&mut epochs, op.q_target);
            }
            OperationType::Z | OperationType::R | OperationType::Hmr => {
                mark_read(&mut epochs, op.q_target);
            }
            OperationType::X
            | OperationType::Neg
            | OperationType::Register
            | OperationType::AppendToRegister
            | OperationType::BitInvert
            | OperationType::BitStore0
            | OperationType::BitStore1
            | OperationType::PushCondition
            | OperationType::PopCondition
            | OperationType::DebugPrint => {}
        }
    }
    out
}

fn mark_read(epochs: &mut [u32], q: QubitId) {
    let idx = q.0 as usize;
    if idx < epochs.len() {
        epochs[idx] = epochs[idx].wrapping_add(1);
    }
}

fn parse_drop_spec() -> DropSpec {
    let mut spec = DropSpec::default();
    for token in env_tokens("BILINEAR_TOFFOLI_DROP") {
        match token.parse::<usize>() {
            Ok(idx) => {
                spec.singles.insert(idx);
            }
            Err(_) => eprintln!("BILINEAR_TOFFOLI_DROP: bad index '{token}'"),
        }
    }
    for token in env_tokens("BILINEAR_TOFFOLI_DROP_PAIRS") {
        let Some((a, b)) = token.split_once(':') else {
            eprintln!("BILINEAR_TOFFOLI_DROP_PAIRS: expected i:j, got '{token}'");
            continue;
        };
        let Ok(a) = a.parse::<usize>() else {
            eprintln!("BILINEAR_TOFFOLI_DROP_PAIRS: bad index '{a}'");
            continue;
        };
        let Ok(b) = b.parse::<usize>() else {
            eprintln!("BILINEAR_TOFFOLI_DROP_PAIRS: bad index '{b}'");
            continue;
        };
        spec.pairs.push((a, b));
    }
    for token in env_tokens("BILINEAR_TOFFOLI_DROP_GROUPS") {
        let mut group = Vec::new();
        let mut ok = true;
        for field in token.split(':') {
            match field.parse::<usize>() {
                Ok(idx) => group.push(idx),
                Err(_) => {
                    eprintln!("BILINEAR_TOFFOLI_DROP_GROUPS: bad group '{token}'");
                    ok = false;
                    break;
                }
            }
        }
        if ok && group.len() >= 2 {
            spec.groups.push(group);
        }
    }
    spec
}

fn env_tokens(env_name: &str) -> Vec<String> {
    let mut text = String::new();
    if let Ok(value) = std::env::var(env_name) {
        text.push_str(&value);
        text.push('\n');
    }
    let file_env = format!("{env_name}_FILE");
    if let Ok(path) = std::env::var(file_env) {
        match std::fs::read_to_string(&path) {
            Ok(value) => {
                text.push_str(&value);
                text.push('\n');
            }
            Err(error) => eprintln!("{env_name}: failed to read {path}: {error}"),
        }
    }
    text.lines()
        .flat_map(|line| {
            line.split('#')
                .next()
                .unwrap_or("")
                .split(|ch: char| ch == ',' || ch.is_whitespace())
                .filter(|s| !s.is_empty())
                .map(str::to_owned)
                .collect::<Vec<_>>()
        })
        .collect()
}

fn env_usize(key: &str, default: usize) -> usize {
    std::env::var(key)
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(default)
}

fn toffoli_ops(ops: &[Op]) -> usize {
    ops.iter()
        .filter(|op| matches!(op.kind, OperationType::CCX | OperationType::CCZ))
        .count()
}

fn mix64(mut x: u64) -> u64 {
    x ^= x >> 30;
    x = x.wrapping_mul(0xbf58_476d_1ce4_e5b9);
    x ^= x >> 27;
    x = x.wrapping_mul(0x94d0_49bb_1331_11eb);
    x ^ (x >> 31)
}

fn linear_mix_mask(mask: u64, block_id: u64) -> (u64, u64) {
    let r1 = ((block_id.wrapping_mul(13).wrapping_add(7)) & 63) as u32;
    let r2 = ((block_id.wrapping_mul(29).wrapping_add(19)) & 63) as u32;
    let r3 = ((block_id.wrapping_mul(43).wrapping_add(31)) & 63) as u32;
    let r4 = ((block_id.wrapping_mul(53).wrapping_add(47)) & 63) as u32;
    (
        mask.rotate_left(r1) ^ mask.rotate_left(r2),
        mask.rotate_left(r3) ^ mask.rotate_left(r4),
    )
}

fn official_xof(ops: &[Op]) -> sha3::Shake256Reader {
    let mut hasher = Shake256::default();
    hasher.update(b"quantum_ecc-fiat-shamir-v2");
    hasher.update(&(ops.len() as u64).to_le_bytes());
    for op in ops {
        absorb_op(&mut hasher, op);
    }
    hasher.finalize_xof()
}

fn support_profile_xof(ops: &[Op], tag: u64) -> sha3::Shake256Reader {
    let mut hasher = Shake256::default();
    hasher.update(b"quantum_ecc-bilinear-toffoli-profile-v1");
    hasher.update(&tag.to_le_bytes());
    hasher.update(&(ops.len() as u64).to_le_bytes());
    for op in ops {
        absorb_op(&mut hasher, op);
    }
    hasher.finalize_xof()
}

fn absorb_op(hasher: &mut Shake256, op: &Op) {
    hasher.update(&[op.kind as u8]);
    hasher.update(&op.q_control2.0.to_le_bytes());
    hasher.update(&op.q_control1.0.to_le_bytes());
    hasher.update(&op.q_target.0.to_le_bytes());
    hasher.update(&op.c_target.0.to_le_bytes());
    hasher.update(&op.c_condition.0.to_le_bytes());
    hasher.update(&op.r_target.0.to_le_bytes());
}

fn secp256k1() -> WeierstrassEllipticCurve {
    WeierstrassEllipticCurve {
        modulus: U256::from_str_radix(
            "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F",
            16,
        )
        .unwrap(),
        a: U256::from(0),
        b: U256::from(7),
        gx: U256::from_str_radix(
            "79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798",
            16,
        )
        .unwrap(),
        gy: U256::from_str_radix(
            "483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8",
            16,
        )
        .unwrap(),
        order: U256::from_str_radix(
            "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141",
            16,
        )
        .unwrap(),
    }
}
