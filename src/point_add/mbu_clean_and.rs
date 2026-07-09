//! Exact measurement-based uncompute for clean AND pairs.
//!
//! Rewrites the clearing half of
//!
//!   CCX(a,b,t); ...; CCX(a,b,t)
//!
//! to
//!
//!   HMR(t -> m); CZ(a,b) if m
//!
//! when `t` is a clean zero ancilla at the compute, `t` is not written until
//! the clear, and `a,b` are not written between the pair. This is Gidney's
//! temporary logical-AND uncompute: it preserves the quantum action exactly
//! while removing one counted Toffoli.

use std::collections::{BTreeMap, BTreeSet};

use crate::circuit::{BitId, Op, OperationType, NO_BIT, NO_QUBIT};

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
struct CcxKey {
    a: u64,
    b: u64,
    target: u64,
}

pub(crate) fn run(ops: Vec<Op>) -> Vec<Op> {
    if std::env::var("MBU_CLEAN_AND_ENABLE").ok().as_deref() != Some("1") {
        return ops;
    }
    if std::env::var("MBU_CLEAN_AND_DISABLE").ok().as_deref() == Some("1") {
        return ops;
    }

    let mut candidates = conservative_immediate_candidates(&ops);
    if let Some(clear) = env_usize("MBU_CLEAN_AND_ONLY_CLEAR") {
        candidates.retain(|candidate| candidate.clear == clear);
    }
    if let Some(limit) = env_usize("MBU_CLEAN_AND_LIMIT") {
        candidates.truncate(limit);
    }
    if std::env::var("MBU_CLEAN_AND_TRACE").ok().as_deref() == Some("1") {
        for candidate in candidates
            .iter()
            .take(env_usize("MBU_CLEAN_AND_TRACE_LIMIT").unwrap_or(16))
        {
            eprintln!(
                "MBU_CLEAN_AND_CAND clear={} reset={} key=(q{},q{})->q{}",
                candidate.clear,
                candidate.reset,
                candidate.key.a,
                candidate.key.b,
                candidate.key.target
            );
        }
    }
    if candidates.is_empty() {
        return ops;
    }

    let before_tof = toffoli_ops(&ops);
    let before_ops = ops.len();
    let mut next_bit = next_fresh_bit(&ops);
    let mut rewritten = 0usize;
    let clear_indices: BTreeMap<usize, CcxKey> = candidates
        .iter()
        .map(|candidate| (candidate.clear, candidate.key))
        .collect();
    let reset_indices: BTreeSet<usize> =
        candidates.iter().map(|candidate| candidate.reset).collect();
    let mut out = Vec::with_capacity(ops.len() + candidates.len());

    for (idx, op) in ops.into_iter().enumerate() {
        if let Some(key) = clear_indices.get(&idx).copied() {
            let bit = BitId(next_bit);
            next_bit += 1;
            emit_immediate_mbu_clear(&mut out, op, key, bit);
            rewritten += 1;
            continue;
        }
        if reset_indices.contains(&idx) {
            continue;
        }
        out.push(op);
    }

    eprintln!(
        "MBU_CLEAN_AND: rewritten={} ops {}->{} emitted_tof {}->{}",
        rewritten,
        before_ops,
        out.len(),
        before_tof,
        toffoli_ops(&out),
    );
    out
}

#[derive(Clone, Copy, Debug)]
struct DeferredCandidate {
    clear: usize,
    reset: usize,
    key: CcxKey,
}

fn conservative_immediate_candidates(ops: &[Op]) -> Vec<DeferredCandidate> {
    let writes = writes_by_qubit(ops);
    let touches = touches_by_qubit(ops);
    let depth_before = condition_depth_before(ops);
    let mut out = Vec::new();
    for (idx, op) in ops.iter().enumerate() {
        if op.kind != OperationType::CCX || op.c_condition != NO_BIT || depth_before[idx] != 0 {
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
        if ops[next].kind != OperationType::CCX
            || ops[next].c_condition != NO_BIT
            || depth_before[next] != 0
            || ccx_key(ops[next]) != Some(key)
        {
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
        let Some(reset) = next_touch(&touches, key.target, next) else {
            continue;
        };
        if ops[reset].kind != OperationType::R
            || ops[reset].c_condition != NO_BIT
            || depth_before[reset] != 0
        {
            continue;
        }
        if has_measure_between(ops, next, reset) {
            continue;
        }
        out.push(DeferredCandidate {
            clear: next,
            reset,
            key,
        });
    }
    out
}

fn emit_immediate_mbu_clear(out: &mut Vec<Op>, clear: Op, key: CcxKey, bit: BitId) {
    debug_assert_eq!(clear.kind, OperationType::CCX);

    let mut hmr = Op::empty();
    hmr.kind = OperationType::Hmr;
    hmr.q_target = clear.q_target;
    hmr.c_target = bit;
    out.push(hmr);

    let mut cz = Op::empty();
    cz.kind = OperationType::CZ;
    cz.q_control1 = crate::circuit::QubitId(key.a);
    cz.q_target = crate::circuit::QubitId(key.b);
    cz.c_condition = bit;
    out.push(cz);
}

fn next_fresh_bit(ops: &[Op]) -> u64 {
    ops.iter()
        .flat_map(|op| [op.c_target, op.c_condition])
        .filter(|bit| *bit != NO_BIT)
        .map(|bit| bit.0)
        .max()
        .map_or(0, |bit| bit + 1)
}

fn env_usize(name: &str) -> Option<usize> {
    std::env::var(name)
        .ok()
        .and_then(|value| value.parse::<usize>().ok())
}

fn toffoli_ops(ops: &[Op]) -> usize {
    ops.iter()
        .filter(|op| matches!(op.kind, OperationType::CCX | OperationType::CCZ))
        .count()
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

fn next_touch(touches: &[Vec<usize>], q: u64, after: usize) -> Option<usize> {
    let row = touches.get(q as usize)?;
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

fn has_measure_between(ops: &[Op], lo: usize, hi: usize) -> bool {
    ops[lo.saturating_add(1)..hi]
        .iter()
        .any(|op| matches!(op.kind, OperationType::R | OperationType::Hmr))
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
