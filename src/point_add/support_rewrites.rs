//! Experimental support-restricted nonlinear rewrites.
//!
//! This module is deliberately opt-in. It exists to test two ideas:
//!
//! 1. Replace a `CCX`/`CCZ` with a Clifford expression that is wrong on one
//!    specified control pattern, after profiling shows that pattern is absent
//!    on the verifier-shaped support.
//! 2. Replace one quantum control by a known classical `BitId`, turning
//!    `CCX(q, bit, t)` into a bit-conditioned `CX` and `CCZ(a, b, bit)` into a
//!    bit-conditioned `CZ`.
//!
//! The pass is not enabled by default because these are support claims, not
//! purely syntactic identities. Use the profiler first, then validate every
//! applied candidate with the trusted evaluator.

use alloy_primitives::U256;
use sha3::{
    digest::{ExtendableOutput, Update, XofReader},
    Shake256,
};

use crate::circuit::{analyze_ops, BitId, Op, OperationType, QubitId, NO_BIT, NO_QUBIT, NO_REG};
use crate::sim::Simulator;
use crate::weierstrass_elliptic_curve::WeierstrassEllipticCurve;

const DEFAULT_PROFILE_SHOTS: usize = 9024;
const DEFAULT_PROFILE_TOP: usize = 40;
const NO_SLOT: u32 = u32::MAX;

#[derive(Clone, Debug)]
struct Profile {
    op_index: usize,
    kind: OperationType,
    executed: u64,
    effect: u64,
    ccx_patterns: [u64; 4],
    ccz_patterns: [u64; 8],
}

#[derive(Clone, Copy, Debug)]
struct BitCtrlRewrite {
    slot: u8,
    bit: BitId,
}

/// Run the opt-in profiler and support-rewrite pass.
pub fn run(ops: Vec<Op>) -> Vec<Op> {
    if std::env::var("SUPPORT_REWRITE_PROFILE").ok().as_deref() == Some("1") {
        profile_and_print(&ops);
    }

    let ccx_missing = parse_missing_map("SUPPORT_REWRITE_CCX_MISSING");
    let ccz_missing = parse_missing_map("SUPPORT_REWRITE_CCZ_MISSING");
    let ccx_bitctrl = parse_bitctrl_map("SUPPORT_REWRITE_CCX_BITCTRL");
    let ccz_bitctrl = parse_bitctrl_map("SUPPORT_REWRITE_CCZ_BITCTRL");

    if ccx_missing.is_empty()
        && ccz_missing.is_empty()
        && ccx_bitctrl.is_empty()
        && ccz_bitctrl.is_empty()
    {
        return ops;
    }

    let before_tof = toffoli_ops(&ops);
    let before_ops = ops.len();
    let mut out = Vec::with_capacity(ops.len());
    let mut stats = RewriteStats::default();

    for (i, op) in ops.iter().copied().enumerate() {
        if let Some(rewrite) = ccx_bitctrl.get(&i) {
            if op.kind == OperationType::CCX {
                emit_ccx_bitctrl(&mut out, op, *rewrite);
                stats.ccx_bitctrl += 1;
                continue;
            }
            eprintln!(
                "SUPPORT_REWRITE_CCX_BITCTRL: index {i} is {:?}, not CCX; keeping op",
                op.kind
            );
        }

        if let Some(rewrite) = ccz_bitctrl.get(&i) {
            if op.kind == OperationType::CCZ {
                emit_ccz_bitctrl(&mut out, op, *rewrite);
                stats.ccz_bitctrl += 1;
                continue;
            }
            eprintln!(
                "SUPPORT_REWRITE_CCZ_BITCTRL: index {i} is {:?}, not CCZ; keeping op",
                op.kind
            );
        }

        if let Some(&missing) = ccx_missing.get(&i) {
            if op.kind == OperationType::CCX {
                if missing < 4 {
                    emit_ccx_missing_pattern(&mut out, op, missing);
                    stats.ccx_missing += 1;
                    continue;
                }
                eprintln!(
                    "SUPPORT_REWRITE_CCX_MISSING: index {i} has invalid missing pattern {missing}; keeping op",
                );
            } else {
                eprintln!(
                    "SUPPORT_REWRITE_CCX_MISSING: index {i} is {:?}, not CCX; keeping op",
                    op.kind
                );
            }
        }

        if let Some(&missing) = ccz_missing.get(&i) {
            if op.kind == OperationType::CCZ {
                if missing < 8 {
                    emit_ccz_missing_pattern(&mut out, op, missing);
                    stats.ccz_missing += 1;
                    continue;
                }
                eprintln!(
                    "SUPPORT_REWRITE_CCZ_MISSING: index {i} has invalid missing pattern {missing}; keeping op",
                );
            } else {
                eprintln!(
                    "SUPPORT_REWRITE_CCZ_MISSING: index {i} is {:?}, not CCZ; keeping op",
                    op.kind
                );
            }
        }

        out.push(op);
    }

    let after_tof = toffoli_ops(&out);
    eprintln!(
        "SUPPORT_REWRITE: ccx_missing={} ccz_missing={} ccx_bitctrl={} ccz_bitctrl={} ops {}->{} emitted_tof {}->{}",
        stats.ccx_missing,
        stats.ccz_missing,
        stats.ccx_bitctrl,
        stats.ccz_bitctrl,
        before_ops,
        out.len(),
        before_tof,
        after_tof,
    );

    if std::env::var("SUPPORT_REWRITE_PREFILTER").ok().as_deref() == Some("1") {
        prefilter_missing_patterns(&ops, &out, &ccx_missing, &ccz_missing);
    }

    out
}

#[derive(Default)]
struct RewriteStats {
    ccx_missing: usize,
    ccz_missing: usize,
    ccx_bitctrl: usize,
    ccz_bitctrl: usize,
}

fn toffoli_ops(ops: &[Op]) -> usize {
    ops.iter()
        .filter(|op| matches!(op.kind, OperationType::CCX | OperationType::CCZ))
        .count()
}

fn profile_and_print(ops: &[Op]) {
    let shots = std::env::var("SUPPORT_REWRITE_PROFILE_SHOTS")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(DEFAULT_PROFILE_SHOTS);
    let top = std::env::var("SUPPORT_REWRITE_PROFILE_TOP")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(DEFAULT_PROFILE_TOP);
    let seeds = std::env::var("SUPPORT_REWRITE_PROFILE_SEEDS")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(1)
        .max(1);

    let (total_qubits, num_bits, _num_regs, regs) = analyze_ops(ops.iter());
    if regs.len() != 4 || regs.iter().any(|r| r.len() != 256) {
        eprintln!(
            "SUPPORT_PROFILE: expected four 256-wide registers, got {}; skipping",
            regs.len()
        );
        return;
    }

    let mut gate_slots = vec![NO_SLOT; ops.len()];
    let mut profiles = Vec::new();
    for (i, op) in ops.iter().enumerate() {
        if matches!(op.kind, OperationType::CCX | OperationType::CCZ) {
            gate_slots[i] = profiles.len() as u32;
            profiles.push(Profile {
                op_index: i,
                kind: op.kind,
                executed: 0,
                effect: 0,
                ccx_patterns: [0; 4],
                ccz_patterns: [0; 8],
            });
        }
    }

    eprintln!(
        "SUPPORT_PROFILE: seeds={} shots_per_seed={} ops={} nonlinear_gates={} qubits={} bits={}",
        seeds,
        shots,
        ops.len(),
        profiles.len(),
        total_qubits,
        num_bits,
    );

    let curve = secp256k1();
    let mut total_profiled_shots = 0usize;

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
            replay_with_profile(&mut sim, ops, &gate_slots, &mut profiles, live_mask);
        }
    }

    let denom = total_profiled_shots.max(1) as f64;
    print_profile_candidates(&profiles, denom, top);
}

fn replay_with_profile<R: XofReader>(
    sim: &mut Simulator<'_, R>,
    ops: &[Op],
    gate_slots: &[u32],
    profiles: &mut [Profile],
    live_mask: u64,
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
            let profile = &mut profiles[slot as usize];
            profile.executed += cond.count_ones() as u64;
            match op.kind {
                OperationType::CCX => {
                    let a = sim.qubit(op.q_control1);
                    let b = sim.qubit(op.q_control2);
                    profile.ccx_patterns[0] += (cond & !a & !b).count_ones() as u64;
                    profile.ccx_patterns[1] += (cond & !a & b).count_ones() as u64;
                    profile.ccx_patterns[2] += (cond & a & !b).count_ones() as u64;
                    profile.ccx_patterns[3] += (cond & a & b).count_ones() as u64;
                    profile.effect += (cond & a & b).count_ones() as u64;
                }
                OperationType::CCZ => {
                    let a = sim.qubit(op.q_control1);
                    let b = sim.qubit(op.q_control2);
                    let c = sim.qubit(op.q_target);
                    for pattern in 0..8 {
                        let mut mask = cond;
                        mask &= if (pattern & 0b100) == 0 { !a } else { a };
                        mask &= if (pattern & 0b010) == 0 { !b } else { b };
                        mask &= if (pattern & 0b001) == 0 { !c } else { c };
                        profile.ccz_patterns[pattern] += mask.count_ones() as u64;
                    }
                    profile.effect += (cond & a & b & c).count_ones() as u64;
                }
                _ => {}
            }
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

fn prefilter_missing_patterns(
    original_ops: &[Op],
    rewritten_ops: &[Op],
    ccx_missing: &std::collections::BTreeMap<usize, u8>,
    ccz_missing: &std::collections::BTreeMap<usize, u8>,
) {
    if ccx_missing.is_empty() && ccz_missing.is_empty() {
        return;
    }
    let shots = std::env::var("SUPPORT_REWRITE_PREFILTER_SHOTS")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(DEFAULT_PROFILE_SHOTS);
    let max_idx = ccx_missing
        .keys()
        .chain(ccz_missing.keys())
        .copied()
        .max()
        .unwrap_or(0);
    let (total_qubits, num_bits, _num_regs, regs) = analyze_ops(rewritten_ops.iter());
    if regs.len() != 4 || regs.iter().any(|r| r.len() != 256) {
        eprintln!("SUPPORT_PREFILTER: expected four 256-wide registers; skipping");
        return;
    }

    let curve = secp256k1();
    let mut xof = official_xof(rewritten_ops);
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

    let mut hits: std::collections::BTreeMap<usize, u64> = ccx_missing
        .keys()
        .chain(ccz_missing.keys())
        .map(|idx| (*idx, 0))
        .collect();
    let mut execs: std::collections::BTreeMap<usize, u64> =
        hits.keys().map(|idx| (*idx, 0)).collect();

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
        replay_prefilter_until(
            &mut sim,
            original_ops,
            max_idx,
            ccx_missing,
            ccz_missing,
            &mut hits,
            &mut execs,
            live_mask,
        );
    }

    eprintln!(
        "SUPPORT_PREFILTER: shots={} max_idx={} candidates={}",
        targets.len(),
        max_idx,
        hits.len()
    );
    for (idx, hit) in hits {
        let exec = execs.get(&idx).copied().unwrap_or(0);
        eprintln!(
            "SUPPORT_PREFILTER idx={} exec={} missing_hits={}",
            idx, exec, hit
        );
    }
}

#[allow(clippy::too_many_arguments)]
fn replay_prefilter_until<R: XofReader>(
    sim: &mut Simulator<'_, R>,
    ops: &[Op],
    max_idx: usize,
    ccx_missing: &std::collections::BTreeMap<usize, u8>,
    ccz_missing: &std::collections::BTreeMap<usize, u8>,
    hits: &mut std::collections::BTreeMap<usize, u64>,
    execs: &mut std::collections::BTreeMap<usize, u64>,
    live_mask: u64,
) {
    let mut condition_stack = Vec::new();
    let mut current_base_condition = live_mask;

    for (i, op) in ops.iter().enumerate() {
        if i > max_idx {
            break;
        }
        let mut cond = current_base_condition;
        if op.c_condition != NO_BIT {
            cond &= sim.bit(op.c_condition);
        }
        cond &= live_mask;

        if let Some(&missing) = ccx_missing.get(&i) {
            let a = sim.qubit(op.q_control1);
            let b = sim.qubit(op.q_control2);
            let mut mask = cond;
            mask &= if (missing & 0b10) == 0 { !a } else { a };
            mask &= if (missing & 0b01) == 0 { !b } else { b };
            *hits.entry(i).or_default() += mask.count_ones() as u64;
            *execs.entry(i).or_default() += cond.count_ones() as u64;
        }
        if let Some(&missing) = ccz_missing.get(&i) {
            let a = sim.qubit(op.q_control1);
            let b = sim.qubit(op.q_control2);
            let c = sim.qubit(op.q_target);
            let mut mask = cond;
            mask &= if (missing & 0b100) == 0 { !a } else { a };
            mask &= if (missing & 0b010) == 0 { !b } else { b };
            mask &= if (missing & 0b001) == 0 { !c } else { c };
            *hits.entry(i).or_default() += mask.count_ones() as u64;
            *execs.entry(i).or_default() += cond.count_ones() as u64;
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

fn print_profile_candidates(profiles: &[Profile], denom: f64, top: usize) {
    let mut candidates = Vec::new();
    let mut rare_effect = Vec::new();

    for p in profiles {
        if p.executed == 0 {
            continue;
        }
        match p.kind {
            OperationType::CCX => {
                for missing in 0..4 {
                    if p.ccx_patterns[missing] == 0 {
                        candidates.push((p.executed, p.effect, p.op_index, p.kind, missing as u8));
                    }
                }
            }
            OperationType::CCZ => {
                for missing in 0..8 {
                    if p.ccz_patterns[missing] == 0 {
                        candidates.push((p.executed, p.effect, p.op_index, p.kind, missing as u8));
                    }
                }
            }
            _ => {}
        }
        rare_effect.push((p.effect, p.executed, p.op_index, p.kind));
    }

    candidates.sort_by(|a, b| b.0.cmp(&a.0).then_with(|| a.2.cmp(&b.2)));
    rare_effect.sort_by(|a, b| a.0.cmp(&b.0).then_with(|| b.1.cmp(&a.1)));

    eprintln!("=== support-rewrite one-missing-pattern candidates ===");
    eprintln!(
        "{:<12} {:<4} {:<7} {:>12} {:>12} {:>10} {:>10}",
        "index", "kind", "missing", "exec_sum", "effect_sum", "exec/shot", "effect%"
    );
    for (executed, effect, index, kind, missing) in candidates.into_iter().take(top) {
        let effect_pct = if executed == 0 {
            0.0
        } else {
            (effect as f64) * 100.0 / (executed as f64)
        };
        eprintln!(
            "{:<12} {:<4} {:<7} {:>12} {:>12} {:>10.3} {:>9.5}%",
            index,
            kind_name(kind),
            missing,
            executed,
            effect,
            (executed as f64) / denom,
            effect_pct,
        );
    }

    eprintln!("=== rare actual-effect nonlinear gates ===");
    eprintln!(
        "{:<12} {:<4} {:>12} {:>12} {:>10} {:>10}",
        "index", "kind", "effect_sum", "exec_sum", "exec/shot", "effect%"
    );
    for (effect, executed, index, kind) in rare_effect.into_iter().take(top) {
        let effect_pct = if executed == 0 {
            0.0
        } else {
            (effect as f64) * 100.0 / (executed as f64)
        };
        eprintln!(
            "{:<12} {:<4} {:>12} {:>12} {:>10.3} {:>9.5}%",
            index,
            kind_name(kind),
            effect,
            executed,
            (executed as f64) / denom,
            effect_pct,
        );
    }
}

fn kind_name(kind: OperationType) -> &'static str {
    match kind {
        OperationType::CCX => "CCX",
        OperationType::CCZ => "CCZ",
        _ => "op",
    }
}

fn emit_ccx_missing_pattern(out: &mut Vec<Op>, op: Op, missing: u8) {
    let a = op.q_control1;
    let b = op.q_control2;
    match missing {
        // Wrong only at 00: 1 ^ a ^ b.
        0 => {
            push_x(out, op.q_target, op.c_condition);
            push_cx(out, a, op.q_target, op.c_condition);
            push_cx(out, b, op.q_target, op.c_condition);
        }
        // Wrong only at 01: b.
        1 => push_cx(out, b, op.q_target, op.c_condition),
        // Wrong only at 10: a.
        2 => push_cx(out, a, op.q_target, op.c_condition),
        // Wrong only at 11: 0.
        3 => {}
        _ => unreachable!(),
    }
}

fn emit_ccz_missing_pattern(out: &mut Vec<Op>, op: Op, missing: u8) {
    let a = op.q_control1;
    let b = op.q_control2;
    let c = op.q_target;
    let cond = op.c_condition;
    let vars = [a, b, c];

    for bit in 0..3 {
        if (missing & (1 << (2 - bit))) != 0 {
            push_x(out, vars[bit], cond);
        }
    }

    push_neg(out, cond);
    push_z(out, a, cond);
    push_z(out, b, cond);
    push_z(out, c, cond);
    push_cz(out, a, b, cond);
    push_cz(out, a, c, cond);
    push_cz(out, b, c, cond);

    for bit in (0..3).rev() {
        if (missing & (1 << (2 - bit))) != 0 {
            push_x(out, vars[bit], cond);
        }
    }
}

fn emit_ccx_bitctrl(out: &mut Vec<Op>, op: Op, rewrite: BitCtrlRewrite) {
    let keep = match rewrite.slot {
        1 => op.q_control2,
        2 => op.q_control1,
        _ => {
            eprintln!(
                "SUPPORT_REWRITE_CCX_BITCTRL: invalid control slot {}; keeping original CCX",
                rewrite.slot
            );
            out.push(op);
            return;
        }
    };
    push_condition(out, rewrite.bit);
    push_cx(out, keep, op.q_target, op.c_condition);
    push_pop(out);
}

fn emit_ccz_bitctrl(out: &mut Vec<Op>, op: Op, rewrite: BitCtrlRewrite) {
    let (left, right) = match rewrite.slot {
        1 => (op.q_control2, op.q_target),
        2 => (op.q_control1, op.q_target),
        3 => (op.q_control1, op.q_control2),
        _ => {
            eprintln!(
                "SUPPORT_REWRITE_CCZ_BITCTRL: invalid control slot {}; keeping original CCZ",
                rewrite.slot
            );
            out.push(op);
            return;
        }
    };
    push_condition(out, rewrite.bit);
    push_cz(out, left, right, op.c_condition);
    push_pop(out);
}

fn push_x(out: &mut Vec<Op>, target: QubitId, cond: BitId) {
    let mut op = empty_op(OperationType::X);
    op.q_target = target;
    op.c_condition = cond;
    out.push(op);
}

fn push_z(out: &mut Vec<Op>, target: QubitId, cond: BitId) {
    let mut op = empty_op(OperationType::Z);
    op.q_target = target;
    op.c_condition = cond;
    out.push(op);
}

fn push_cx(out: &mut Vec<Op>, ctrl: QubitId, target: QubitId, cond: BitId) {
    if ctrl == target {
        eprintln!(
            "SUPPORT_REWRITE: invalid CX alias q{}; dropping replacement op",
            ctrl.0
        );
        return;
    }
    let mut op = empty_op(OperationType::CX);
    op.q_control1 = ctrl;
    op.q_target = target;
    op.c_condition = cond;
    out.push(op);
}

fn push_cz(out: &mut Vec<Op>, a: QubitId, b: QubitId, cond: BitId) {
    if a == b {
        push_z(out, a, cond);
        return;
    }
    let mut op = empty_op(OperationType::CZ);
    op.q_control1 = a;
    op.q_target = b;
    op.c_condition = cond;
    out.push(op);
}

fn push_neg(out: &mut Vec<Op>, cond: BitId) {
    let mut op = empty_op(OperationType::Neg);
    op.c_condition = cond;
    out.push(op);
}

fn push_condition(out: &mut Vec<Op>, cond: BitId) {
    let mut op = empty_op(OperationType::PushCondition);
    op.c_condition = cond;
    out.push(op);
}

fn push_pop(out: &mut Vec<Op>) {
    out.push(empty_op(OperationType::PopCondition));
}

fn empty_op(kind: OperationType) -> Op {
    Op {
        kind,
        q_control2: NO_QUBIT,
        q_control1: NO_QUBIT,
        q_target: NO_QUBIT,
        c_target: NO_BIT,
        c_condition: NO_BIT,
        r_target: NO_REG,
    }
}

fn parse_missing_map(env_name: &str) -> std::collections::BTreeMap<usize, u8> {
    let mut out = std::collections::BTreeMap::new();
    for token in env_tokens(env_name) {
        let Some((idx, missing)) = token.split_once(':') else {
            eprintln!("{env_name}: expected idx:missing, got '{token}'");
            continue;
        };
        let Ok(idx) = idx.parse::<usize>() else {
            eprintln!("{env_name}: bad index '{idx}'");
            continue;
        };
        let Ok(missing) = missing.parse::<u8>() else {
            eprintln!("{env_name}: bad missing pattern '{missing}'");
            continue;
        };
        out.insert(idx, missing);
    }
    out
}

fn parse_bitctrl_map(env_name: &str) -> std::collections::BTreeMap<usize, BitCtrlRewrite> {
    let mut out = std::collections::BTreeMap::new();
    for token in env_tokens(env_name) {
        let fields: Vec<&str> = token.split(':').collect();
        if fields.len() != 3 {
            eprintln!("{env_name}: expected idx:slot:bit, got '{token}'");
            continue;
        }
        let Ok(idx) = fields[0].parse::<usize>() else {
            eprintln!("{env_name}: bad index '{}'", fields[0]);
            continue;
        };
        let Ok(slot) = fields[1].parse::<u8>() else {
            eprintln!("{env_name}: bad slot '{}'", fields[1]);
            continue;
        };
        let Ok(bit) = fields[2].parse::<u64>() else {
            eprintln!("{env_name}: bad bit '{}'", fields[2]);
            continue;
        };
        out.insert(
            idx,
            BitCtrlRewrite {
                slot,
                bit: BitId(bit),
            },
        );
    }
    out
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
    let mut tokens = Vec::new();
    for line in text.lines() {
        let line = line.split('#').next().unwrap_or("").trim();
        tokens.extend(
            line.split(|ch: char| ch == ',' || ch.is_whitespace())
                .filter(|s| !s.is_empty())
                .map(str::to_owned),
        );
    }
    tokens
}

fn official_xof(ops: &[Op]) -> sha3::Shake256Reader {
    let mut hasher = Shake256::default();
    hasher.update(b"quantum_ecc-fiat-shamir-v2");
    hasher.update(&(ops.len() as u64).to_le_bytes());
    for op in ops {
        hasher.update(&[op.kind as u8]);
        hasher.update(&op.q_control2.0.to_le_bytes());
        hasher.update(&op.q_control1.0.to_le_bytes());
        hasher.update(&op.q_target.0.to_le_bytes());
        hasher.update(&op.c_target.0.to_le_bytes());
        hasher.update(&op.c_condition.0.to_le_bytes());
        hasher.update(&op.r_target.0.to_le_bytes());
    }
    hasher.finalize_xof()
}

fn support_profile_xof(ops: &[Op], tag: u64) -> sha3::Shake256Reader {
    let mut hasher = Shake256::default();
    hasher.update(b"quantum_ecc-support-profile-v1");
    hasher.update(&tag.to_le_bytes());
    hasher.update(&(ops.len() as u64).to_le_bytes());
    for op in ops {
        hasher.update(&[op.kind as u8]);
        hasher.update(&op.q_control2.0.to_le_bytes());
        hasher.update(&op.q_control1.0.to_le_bytes());
        hasher.update(&op.q_target.0.to_le_bytes());
        hasher.update(&op.c_target.0.to_le_bytes());
        hasher.update(&op.c_condition.0.to_le_bytes());
        hasher.update(&op.r_target.0.to_le_bytes());
    }
    hasher.finalize_xof()
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
