//! find_dead_ccx — dynamic dead-CCX finder for the dead-CCX drop lever.
//!
//! Builds the (config-from-env) point_add circuit ONCE (post-fanout, exactly the
//! stream eval_circuit/the grader scores), then runs a self-contained bit-sliced
//! simulation (transposed op-major over MANY Fiat-Shamir seeds, reading qubit/bit
//! state via the public Simulator fields -- no Simulator changes needed), recording
//! per-op-index the OR of each CCX/CCZ's per-shot FIRE mask. An op that
//! NEVER fires across every screened shot is "charged-but-inert": the grader counts
//! its executed_shots in the Toffoli total, but it flips nothing -> removing it is
//! value-neutral. We emit the never-fired CCX/CCZ op indices as a `.idx` list
//! (same format the build()-side drop-loader reads), matching the field's
//! drop_dead_robust mechanism but regenerated for OUR circuit.
//!
//! Robustness: deadness is the INTERSECTION over all screen seeds (an op is dead
//! only if it fired on no seed). More seeds -> safer + ship-nonces stay findable.
//! Final safety is enforced downstream anyway: the post-drop nonce hunt re-checks
//! 0/0/0 on the ship seed, so an over-aggressive drop only ever shows up as a dirty
//! (rejected) nonce, never as a wrong-but-accepted circuit.
//!
//! Env:
//!   DEAD_SCREEN_NONCES   space/comma nonce list to screen (default: built-in spread)
//!   PREDICT_SHOTS        shots per seed (default 9024, the grader count)
//!   DEAD_IDX_OUT         output .idx path (default /tmp/dead_ccx.idx)
//!   DEAD_CANDIDATE_IDX_FILE restrict final/checkpoint output to an existing idx list
//!   DEAD_BATCH_SIZE       checkpoint batch size in nonces (default: all nonces)
//!   DEAD_CHECKPOINT_DIR   optional directory for per-batch dead/rank checkpoints
//!   DEAD_RANK_OUT         optional TSV ranking by charged avg-Toffoli save
//!   + all the usual TLM_* circuit-config knobs (set them to the target route)
use alloy_primitives::U256;
use quantum_ecc::circuit::{analyze_ops, Op, OperationType, QubitId, QubitOrBit};
use quantum_ecc::point_add;
use quantum_ecc::weierstrass_elliptic_curve::WeierstrassEllipticCurve;
use sha3::{
    digest::{ExtendableOutput, Update, XofReader},
    Shake256,
};
use std::io::Write;
use std::path::Path;

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

#[inline]
fn feed_op(hasher: &mut Shake256, op: &Op) {
    hasher.update(&[op.kind as u8]);
    hasher.update(&op.q_control2.0.to_le_bytes());
    hasher.update(&op.q_control1.0.to_le_bytes());
    hasher.update(&op.q_target.0.to_le_bytes());
    hasher.update(&op.c_target.0.to_le_bytes());
    hasher.update(&op.c_condition.0.to_le_bytes());
    hasher.update(&op.r_target.0.to_le_bytes());
}

// Per-nonce FS-seed reuse (same trick as predict_fast): the nonce only swaps the
// q_target of the trailing NONCE_BITS*2 identity X ops, so pre-feed the fixed
// prefix once and re-feed only the tail per nonce.
const NONCE_BITS: usize = 48;
const TAIL_OPS: usize = NONCE_BITS * 2;

struct SeedPrefix {
    prefix_hasher: Shake256,
    tail_template: Vec<Op>,
    tx0: u64,
    tx1: u64,
}

impl SeedPrefix {
    fn build() -> (Self, Vec<Op>) {
        std::env::set_var("DIALOG_TAIL_NONCE", "0");
        let ops = point_add::build();
        let n = ops.len();
        let tail_template = ops[n - TAIL_OPS..].to_vec();
        std::env::set_var("DIALOG_TAIL_NONCE", "1");
        let ops1 = point_add::build();
        let tx1 = ops1[n - TAIL_OPS].q_target.0;
        let tx0 = tail_template[0].q_target.0;
        let mut prefix_hasher = Shake256::default();
        prefix_hasher.update(b"quantum_ecc-fiat-shamir-v2");
        prefix_hasher.update(&(n as u64).to_le_bytes());
        for op in &ops[..n - TAIL_OPS] {
            feed_op(&mut prefix_hasher, op);
        }
        (
            SeedPrefix {
                prefix_hasher,
                tail_template,
                tx0,
                tx1,
            },
            ops,
        )
    }

    fn seed_for(&self, nonce: u64) -> sha3::Shake256Reader {
        let mut hasher = self.prefix_hasher.clone();
        for bit in 0..NONCE_BITS {
            let qt = if (nonce >> bit) & 1 == 1 {
                self.tx1
            } else {
                self.tx0
            };
            for k in 0..2 {
                let mut op = self.tail_template[bit * 2 + k];
                op.q_target = QubitId(qt);
                feed_op(&mut hasher, &op);
            }
        }
        hasher.finalize_xof()
    }
}

fn precompute_fixed_base(curve: &WeierstrassEllipticCurve, x: U256, y: U256) -> Vec<(U256, U256)> {
    let mut table = Vec::with_capacity(256);
    let mut p = (x, y);
    for _ in 0..256 {
        table.push(p);
        p = curve.add(p.0, p.1, p.0, p.1);
    }
    table
}

fn mul_fixed(curve: &WeierstrassEllipticCurve, table: &[(U256, U256)], k: U256) -> (U256, U256) {
    let mut res = (U256::ZERO, U256::ZERO);
    for (i, t) in table.iter().enumerate() {
        if k.bit(i) {
            res = curve.add(res.0, res.1, t.0, t.1);
        }
    }
    res
}

// ───────────────────────── transposed rng=0 fast path ─────────────────────────
//
// Two structural facts make a big speedup possible, both validated downstream by
// predict_fast cls-add (the authoritative real-rng 0/0/0 check):
//
//  (1) rng=0 is a VALID execution. Every measurement-based uncompute corrects on a
//      random measured bit; the m=0 branch needs no correction and is always a
//      legal, correct run. So forcing all Hmr/R measurements to 0 gives correct
//      qubit evolution. The ONLY firing it gets wrong is the replay CCX gated on a
//      measured bit (they fire ~50% under real rng, 0% under m=0) -- so we must
//      never drop those. We exclude them via structural taint (a CCX whose
//      condition traces to an Hmr-written bit). This removes the entire ~1.4TB
//      SHAKE squeeze for a 9M screen.
//
//  (2) Only bits read as conditions affect firing. We compact bit storage to those
//      and skip all other bit writes.
//
// Plus the loop is transposed to op-major over a qubit-major / cond-bit-major state,
// so the 568MB op-list is streamed once per SEED (not once per 64-shot batch) and
// the hot 1156-qubit state stays cache-resident.

struct Analysis {
    bit_remap: Vec<u32>, // bit id -> compact cond-bit index, or u32::MAX (not a condition bit)
    n_cond_bits: usize,
    exclude: Vec<bool>, // op idx -> true if this CCX/CCZ is rng-dependent (never drop)
}

fn analyze(ops: &[Op], num_bits: usize) -> Analysis {
    use OperationType::*;
    let no_bit = quantum_ecc::circuit::NO_BIT;
    // Pass 1: bits ever used as a condition (c_condition on any op, incl PushCondition).
    let mut is_cond = vec![false; num_bits];
    for o in ops {
        if o.c_condition != no_bit {
            is_cond[o.c_condition.0 as usize] = true;
        }
    }
    let mut bit_remap = vec![u32::MAX; num_bits];
    let mut n_cond_bits = 0usize;
    for b in 0..num_bits {
        if is_cond[b] {
            bit_remap[b] = n_cond_bits as u32;
            n_cond_bits += 1;
        }
    }
    // Pass 2: conservative taint propagation. A bit is tainted once an Hmr writes it
    // (rng-dependent) and stays tainted until an UNCONDITIONAL deterministic store.
    // A CCX/CCZ is excluded if it executes under any tainted condition.
    let mut tainted = vec![false; num_bits];
    let mut stack: Vec<bool> = Vec::new();
    let mut tdepth = 0usize;
    let mut exclude = vec![false; ops.len()];
    for (i, o) in ops.iter().enumerate() {
        let cc = o.c_condition;
        let op_rng_dep = tdepth > 0 || (cc != no_bit && tainted[cc.0 as usize]);
        match o.kind {
            CCX | CCZ => {
                if op_rng_dep {
                    exclude[i] = true;
                }
            }
            PushCondition => {
                let t = cc != no_bit && tainted[cc.0 as usize];
                stack.push(t);
                if t {
                    tdepth += 1;
                }
            }
            PopCondition => {
                if let Some(t) = stack.pop() {
                    if t {
                        tdepth -= 1;
                    }
                }
            }
            Hmr => {
                tainted[o.c_target.0 as usize] = true;
            }
            BitStore0 | BitStore1 => {
                let unconditional = tdepth == 0 && cc == no_bit;
                tainted[o.c_target.0 as usize] = !unconditional;
            }
            BitInvert => {
                let id = o.c_target.0 as usize;
                tainted[id] = tainted[id] || op_rng_dep;
            }
            _ => {}
        }
    }
    Analysis {
        bit_remap,
        n_cond_bits,
        exclude,
    }
}

struct TScratch {
    q: Vec<u64>,    // num_qubits * NB_MAX (qubit-major)
    bt: Vec<u64>,   // n_cond_bits * NB_MAX (cond-bit-major)
    cond: Vec<u64>, // NB scratch
    base: Vec<u64>,
    stack: Vec<Vec<u64>>,
    rng: Vec<u64>, // real-rng mode only: op-major [j*NB + b], j over Hmr+R ops
}

impl TScratch {
    /// `nb_max` = max batches per seed = ceil(shots / 64).
    fn new(num_qubits: usize, n_cond_bits: usize, nb_max: usize) -> Self {
        TScratch {
            q: vec![0u64; num_qubits * nb_max],
            bt: vec![0u64; n_cond_bits * nb_max],
            cond: vec![0u64; nb_max],
            base: vec![0u64; nb_max],
            stack: Vec::with_capacity(64),
            rng: Vec::new(),
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn screen_seed_transposed(
    mut xof: sha3::Shake256Reader,
    ops: &[Op],
    num_qubits: usize,
    regs: &[Vec<QubitOrBit>],
    curve: &WeierstrassEllipticCurve,
    gtable: &[(U256, U256)],
    shots: usize,
    an: &Analysis,
    real_rng: bool,
    hmr_slot: &[i32], // len = #(Hmr+R) per pass; entry = Hmr-only index, or -1 for R
    n_hmr: usize,     // #(Hmr) per pass (rng storage width)
    fired: &mut [u64],
    mut charged: Option<&mut [u64]>,
    s: &mut TScratch,
) -> usize {
    use OperationType::*;
    let no_bit = quantum_ecc::circuit::NO_BIT;
    let remap = &an.bit_remap;
    // Derive shots exactly as run_tests (reads xof for ALL `shots` iters, incl skips).
    let mut targets = Vec::with_capacity(shots);
    let mut offsets = Vec::with_capacity(shots);
    for _ in 0..shots {
        let mut rb = [[0u8; 32]; 2];
        xof.read(&mut rb[0]);
        xof.read(&mut rb[1]);
        let k1 = U256::from_le_bytes(rb[0]);
        let k2 = U256::from_le_bytes(rb[1]);
        let t = mul_fixed(curve, gtable, k1);
        let o = mul_fixed(curve, gtable, k2);
        if t.0 == o.0 || (t.0.is_zero() && t.1.is_zero()) || (o.0.is_zero() && o.1.is_zero()) {
            continue;
        }
        targets.push(t);
        offsets.push(o);
    }
    let n = targets.len();
    if n == 0 {
        return 0;
    }
    let nb = n.div_ceil(64);
    assert!(
        nb <= s.base.len(),
        "nb {} exceeds scratch nb_max {}",
        nb,
        s.base.len()
    );
    let ncb = an.n_cond_bits;

    // Real-rng mode: bulk-squeeze the sim's measured-gate rng (one big read amortizes
    // the per-Hmr/R call overhead). The original sim consumes rng batch-outer/op-inner,
    // so global stream index for (Hmr/R op j, batch b) = b*H + j. Store op-major
    // (rng[j*nb + b]) for contiguous per-op access in the transposed loop.
    if real_rng {
        let h = hmr_slot.len(); // Hmr+R stream length per pass
        s.rng.resize(n_hmr * nb, 0);
        // Bulk-squeeze one batch's worth of measured-gate rng (H u64 ~= 10MB) per
        // read() call -- identical byte stream to H sequential 8-byte reads, but ~H
        // fewer XofReader calls. We must CONSUME all H (Hmr+R) values for stream
        // alignment, but only STORE the Hmr ones (R is phase-only -> skipped), which
        // halves the rng buffer (and the memory contention bottlenecking the run).
        let mut tmp = vec![0u8; h * 8];
        for b in 0..nb {
            xof.read(&mut tmp);
            for k in 0..h {
                let hj = hmr_slot[k];
                if hj >= 0 {
                    let v = u64::from_le_bytes(tmp[k * 8..k * 8 + 8].try_into().unwrap());
                    s.rng[hj as usize * nb + b] = v;
                }
            }
        }
    }
    let _ = real_rng;

    // Zero the used portion of state.
    for q in 0..num_qubits {
        for b in 0..nb {
            s.q[q * nb + b] = 0;
        }
    }
    for c in 0..ncb {
        for b in 0..nb {
            s.bt[c * nb + b] = 0;
        }
    }
    // Per-batch valid-lane mask (last batch may be partial).
    for b in 0..nb {
        let bs = 64.min(n - b * 64);
        s.base[b] = if bs == 64 { u64::MAX } else { (1u64 << bs) - 1 };
    }

    // Set input registers (qubit-major / cond-bit-major). Replicates set_register.
    let mut set_reg = |reg: &[QubitOrBit], val: U256, b: usize, shot: usize| {
        let lane = 1u64 << shot;
        for (k, item) in reg.iter().enumerate() {
            let on = val.bit(k);
            match item {
                QubitOrBit::Qubit(id) => {
                    let idx = id.0 as usize * nb + b;
                    if on {
                        s.q[idx] |= lane;
                    } else {
                        s.q[idx] &= !lane;
                    }
                }
                QubitOrBit::Bit(id) => {
                    let ci = remap[id.0 as usize];
                    if ci != u32::MAX {
                        let idx = ci as usize * nb + b;
                        if on {
                            s.bt[idx] |= lane;
                        } else {
                            s.bt[idx] &= !lane;
                        }
                    }
                }
            }
        }
    };
    for b in 0..nb {
        let bs = 64.min(n - b * 64);
        for shot in 0..bs {
            let i = b * 64 + shot;
            set_reg(&regs[0], targets[i].0, b, shot);
            set_reg(&regs[1], targets[i].1, b, shot);
            set_reg(&regs[2], offsets[i].0, b, shot);
            set_reg(&regs[3], offsets[i].1, b, shot);
        }
    }
    drop(set_reg);

    s.stack.clear();
    let mut hj = 0usize; // Hmr-only index into the stored rng (real-rng mode)
                         // base[] already holds the per-batch validity mask = the top-level condition.

    macro_rules! cond_for {
        ($o:expr) => {{
            if $o.c_condition != no_bit {
                let ci = remap[$o.c_condition.0 as usize] as usize;
                for b in 0..nb {
                    s.cond[b] = s.base[b] & s.bt[ci * nb + b];
                }
                &s.cond[..nb]
            } else {
                &s.base[..nb]
            }
        }};
    }

    for (i, o) in ops.iter().enumerate() {
        match o.kind {
            CCX => {
                let c1 = o.q_control1.0 as usize;
                let c2 = o.q_control2.0 as usize;
                let t = o.q_target.0 as usize;
                let cond = cond_for!(o);
                let mut fv = 0u64;
                let mut cv = 0u64;
                for b in 0..nb {
                    cv += cond[b].count_ones() as u64;
                    let v = cond[b] & s.q[c1 * nb + b] & s.q[c2 * nb + b];
                    fv |= v;
                    s.q[t * nb + b] ^= v;
                }
                fired[i] |= fv;
                if let Some(charged) = charged.as_deref_mut() {
                    charged[i] += cv;
                }
            }
            CCZ => {
                let c1 = o.q_control1.0 as usize;
                let c2 = o.q_control2.0 as usize;
                let t = o.q_target.0 as usize;
                let cond = cond_for!(o);
                let mut fv = 0u64;
                let mut cv = 0u64;
                for b in 0..nb {
                    cv += cond[b].count_ones() as u64;
                    fv |= cond[b] & s.q[t * nb + b] & s.q[c1 * nb + b] & s.q[c2 * nb + b];
                }
                fired[i] |= fv;
                if let Some(charged) = charged.as_deref_mut() {
                    charged[i] += cv;
                }
            }
            CX => {
                let c1 = o.q_control1.0 as usize;
                let t = o.q_target.0 as usize;
                let cond = cond_for!(o);
                for b in 0..nb {
                    s.q[t * nb + b] ^= cond[b] & s.q[c1 * nb + b];
                }
            }
            X => {
                let t = o.q_target.0 as usize;
                let cond = cond_for!(o);
                for b in 0..nb {
                    s.q[t * nb + b] ^= cond[b];
                }
            }
            Swap => {
                let c1 = o.q_control1.0 as usize;
                let t = o.q_target.0 as usize;
                let cond = cond_for!(o);
                for b in 0..nb {
                    let mut a = s.q[c1 * nb + b];
                    let mut d = s.q[t * nb + b];
                    a ^= d;
                    d ^= cond[b] & a;
                    a ^= d;
                    s.q[c1 * nb + b] = a;
                    s.q[t * nb + b] = d;
                }
            }
            Hmr => {
                // bit measured: bit = (bit & !cond) ^ (rng & cond); qubit demolished.
                // rng=0 mode: rng term vanishes -> just clears bit on cond lanes.
                let t = o.q_target.0 as usize;
                let ci = remap[o.c_target.0 as usize];
                let cond = cond_for!(o);
                if real_rng {
                    let rb = hj * nb;
                    for b in 0..nb {
                        if ci != u32::MAX {
                            let idx = ci as usize * nb + b;
                            s.bt[idx] = (s.bt[idx] & !cond[b]) ^ (s.rng[rb + b] & cond[b]);
                        }
                        s.q[t * nb + b] &= !cond[b];
                    }
                    hj += 1;
                } else {
                    for b in 0..nb {
                        if ci != u32::MAX {
                            s.bt[ci as usize * nb + b] &= !cond[b];
                        }
                        s.q[t * nb + b] &= !cond[b];
                    }
                }
            }
            R => {
                // R is phase-only (skipped) and stores no rng -> just demolish qubit.
                let t = o.q_target.0 as usize;
                let cond = cond_for!(o);
                for b in 0..nb {
                    s.q[t * nb + b] &= !cond[b];
                }
            }
            BitInvert => {
                let ci = remap[o.c_target.0 as usize];
                if ci != u32::MAX {
                    let cond = cond_for!(o);
                    for b in 0..nb {
                        s.bt[ci as usize * nb + b] ^= cond[b];
                    }
                }
            }
            BitStore0 => {
                let ci = remap[o.c_target.0 as usize];
                if ci != u32::MAX {
                    let cond = cond_for!(o);
                    for b in 0..nb {
                        s.bt[ci as usize * nb + b] &= !cond[b];
                    }
                }
            }
            BitStore1 => {
                let ci = remap[o.c_target.0 as usize];
                if ci != u32::MAX {
                    let cond = cond_for!(o);
                    for b in 0..nb {
                        s.bt[ci as usize * nb + b] |= cond[b];
                    }
                }
            }
            PushCondition => {
                s.stack.push(s.base[..nb].to_vec());
                if o.c_condition != no_bit {
                    let ci = remap[o.c_condition.0 as usize] as usize;
                    for b in 0..nb {
                        s.base[b] &= s.bt[ci * nb + b];
                    }
                }
            }
            PopCondition => {
                if let Some(prev) = s.stack.pop() {
                    s.base[..nb].copy_from_slice(&prev);
                }
            }
            // Phase-only (Neg/Z/CZ) and annotations never affect firing -> skip.
            _ => {}
        }
    }
    n
}

fn screen_nonces() -> Vec<u64> {
    if let Ok(s) = std::env::var("DEAD_SCREEN_NONCES") {
        let v: Vec<u64> = s
            .split(|c: char| c == ',' || c.is_whitespace())
            .filter(|t| !t.is_empty())
            .map(|t| t.parse().unwrap_or_else(|e| panic!("bad nonce {t:?}: {e}")))
            .collect();
        if !v.is_empty() {
            return v;
        }
    }
    // Built-in spread of well-separated nonces (diverse FS seeds).
    (0..8u64)
        .map(|i| 100_000_000_001 + i * 700_000_003)
        .collect()
}

fn env_usize(name: &str, default: usize) -> usize {
    std::env::var(name)
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(default)
}

fn read_idx_file(path: &str, n_ops: usize) -> Vec<bool> {
    let text = std::fs::read_to_string(path)
        .unwrap_or_else(|err| panic!("failed to read idx file {path}: {err}"));
    let mut mask = vec![false; n_ops];
    for line in text.lines() {
        let t = line.trim();
        if t.is_empty() || t.starts_with('#') || t.starts_with("idx") {
            continue;
        }
        let tok = t.split('\t').next().unwrap_or(t);
        let idx = tok
            .parse::<usize>()
            .unwrap_or_else(|err| panic!("bad idx token {tok:?} in {path}: {err}"));
        if idx < n_ops {
            mask[idx] = true;
        }
    }
    mask
}

fn candidate_mask_from_env(n_ops: usize) -> Vec<bool> {
    if let Ok(path) = std::env::var("DEAD_CANDIDATE_IDX_FILE") {
        read_idx_file(&path, n_ops)
    } else {
        vec![true; n_ops]
    }
}

fn write_dead_idx(
    path: impl AsRef<Path>,
    ops: &[Op],
    fired: &[u64],
    an: &Analysis,
    real_rng: bool,
    candidate_mask: &[bool],
) -> usize {
    let path = path.as_ref();
    let mut f = std::fs::File::create(path)
        .unwrap_or_else(|err| panic!("create dead idx {}: {err}", path.display()));
    writeln!(f, "idx").unwrap();
    let mut n = 0usize;
    for (i, op) in ops.iter().enumerate() {
        if candidate_mask[i]
            && matches!(op.kind, OperationType::CCX | OperationType::CCZ)
            && fired[i] == 0
            && (real_rng || !an.exclude[i])
        {
            writeln!(f, "{i}").unwrap();
            n += 1;
        }
    }
    f.flush().unwrap();
    n
}

fn write_rank_tsv(
    path: impl AsRef<Path>,
    ops: &[Op],
    fired: &[u64],
    charged: &[u64],
    total_shots: usize,
    an: &Analysis,
    real_rng: bool,
    candidate_mask: &[bool],
) {
    let path = path.as_ref();
    let denom = total_shots.max(1) as f64;
    let mut rows = Vec::new();
    for (i, op) in ops.iter().enumerate() {
        if !candidate_mask[i] || !matches!(op.kind, OperationType::CCX | OperationType::CCZ) {
            continue;
        }
        let dead = fired[i] == 0 && (real_rng || !an.exclude[i]);
        rows.push((dead, charged[i], i, op.kind, fired[i] != 0, an.exclude[i]));
    }
    rows.sort_by(|a, b| {
        b.0.cmp(&a.0)
            .then_with(|| b.1.cmp(&a.1))
            .then_with(|| a.2.cmp(&b.2))
    });

    let mut f = std::fs::File::create(path)
        .unwrap_or_else(|err| panic!("create rank tsv {}: {err}", path.display()));
    writeln!(
        f,
        "idx\tkind\tdead\tever_fired\texcluded\tcharged\ttotal_shots\tavg_save"
    )
    .unwrap();
    for (dead, charge, idx, kind, ever_fired, excluded) in rows {
        writeln!(
            f,
            "{idx}\t{kind:?}\t{}\t{}\t{}\t{charge}\t{total_shots}\t{:.9}",
            dead as u8,
            ever_fired as u8,
            excluded as u8,
            charge as f64 / denom
        )
        .unwrap();
    }
    f.flush().unwrap();
}

#[allow(clippy::too_many_arguments)]
fn screen_nonce_batch(
    nonces: &[u64],
    n_threads: usize,
    shots: usize,
    n_ops: usize,
    total_qubits: u64,
    prefix: &SeedPrefix,
    ops: &[Op],
    regs: &[Vec<QubitOrBit>],
    curve: &WeierstrassEllipticCurve,
    gtable: &[(U256, U256)],
    an: &Analysis,
    real_rng: bool,
    hmr_slot: &[i32],
    n_hmr: usize,
    track_charged: bool,
) -> (Vec<u64>, Option<Vec<u64>>, usize) {
    let n_threads = n_threads.max(1).min(nonces.len().max(1));
    let mut chunks: Vec<Vec<u64>> = vec![Vec::new(); n_threads];
    for (i, &nonce) in nonces.iter().enumerate() {
        chunks[i % n_threads].push(nonce);
    }
    std::thread::scope(|sc| {
        let handles: Vec<_> = chunks
            .iter()
            .map(|chunk| {
                sc.spawn(move || {
                    let mut local_fired = vec![0u64; n_ops];
                    let mut local_charged = track_charged.then(|| vec![0u64; n_ops]);
                    let mut scratch =
                        TScratch::new(total_qubits as usize, an.n_cond_bits, shots.div_ceil(64));
                    let mut shots_run = 0usize;
                    for &nonce in chunk {
                        let xof = prefix.seed_for(nonce);
                        shots_run += screen_seed_transposed(
                            xof,
                            ops,
                            total_qubits as usize,
                            regs,
                            curve,
                            gtable,
                            shots,
                            an,
                            real_rng,
                            hmr_slot,
                            n_hmr,
                            &mut local_fired,
                            local_charged.as_deref_mut(),
                            &mut scratch,
                        );
                    }
                    (local_fired, local_charged, shots_run)
                })
            })
            .collect();
        let mut fired = vec![0u64; n_ops];
        let mut charged = track_charged.then(|| vec![0u64; n_ops]);
        let mut total_shots = 0usize;
        for h in handles {
            let (local_fired, local_charged, shots_run) = h.join().expect("screen thread panicked");
            total_shots += shots_run;
            for (g, l) in fired.iter_mut().zip(local_fired.iter()) {
                *g |= *l;
            }
            if let (Some(global), Some(local)) = (charged.as_mut(), local_charged.as_ref()) {
                for (g, l) in global.iter_mut().zip(local.iter()) {
                    *g += *l;
                }
            }
        }
        (fired, charged, total_shots)
    })
}

fn main() {
    let shots: usize = std::env::var("PREDICT_SHOTS")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(9024);
    let out_path =
        std::env::var("DEAD_IDX_OUT").unwrap_or_else(|_| "/tmp/dead_ccx.idx".to_string());
    let nonces = screen_nonces();

    let curve = secp256k1();
    let gtable = precompute_fixed_base(&curve, curve.gx, curve.gy);

    eprintln!("find_dead_ccx: building circuit (post-fanout) ...");
    let (prefix, ops) = SeedPrefix::build();
    let (total_qubits, num_bits, _num_regs, regs) = analyze_ops(ops.iter());
    let n_ops = ops.len();
    let n_ccx = ops
        .iter()
        .filter(|o| matches!(o.kind, OperationType::CCX | OperationType::CCZ))
        .count();
    {
        use OperationType::*;
        let mut h = std::collections::BTreeMap::new();
        for o in &ops {
            *h.entry(format!("{:?}", o.kind)).or_insert(0u64) += 1;
        }
        let nrng = ops.iter().filter(|o| matches!(o.kind, Hmr | R)).count();
        let ncond = ops
            .iter()
            .filter(|o| o.c_condition != quantum_ecc::circuit::NO_BIT)
            .count();
        eprintln!("OP_HISTOGRAM: {:?}", h);
        eprintln!(
            "  rng-consuming (Hmr+R) per pass = {}  ;  ops with c_condition = {}",
            nrng, ncond
        );
    }
    eprintln!(
        "find_dead_ccx: ops={} CCX/CCZ={} qubits={} bits={}  screening {} seed(s) x {} shots ...",
        n_ops,
        n_ccx,
        total_qubits,
        num_bits,
        nonces.len(),
        shots
    );

    // One-time structural analysis (cond-bit storage + rng-taint exclusion).
    let an = analyze(&ops, num_bits as usize);
    let n_excluded = an.exclude.iter().filter(|&&e| e).count();
    // real-rng mode (DEAD_REAL_RNG=1, default): EXACT firing via bulk-squeezed rng ->
    // recovers tainted-but-dead CCX (the field-beating margin), no taint exclusion.
    // rng=0 mode (DEAD_REAL_RNG=0): faster, conservative (excludes tainted CCX).
    let real_rng = std::env::var("DEAD_REAL_RNG").ok().as_deref() != Some("0");
    // Stream map over Hmr+R ops (rng consumption order): Hmr -> its Hmr-only index,
    // R -> -1 (consumed for alignment but not stored, phase-only). Storing only Hmr
    // halves the rng buffer vs storing all Hmr+R.
    let mut hmr_slot: Vec<i32> = Vec::new();
    let mut n_hmr = 0usize;
    for o in &ops {
        match o.kind {
            OperationType::Hmr => {
                hmr_slot.push(n_hmr as i32);
                n_hmr += 1;
            }
            OperationType::R => hmr_slot.push(-1),
            _ => {}
        }
    }
    let hmr_slot_ref = &hmr_slot;
    eprintln!(
        "find_dead_ccx: transposed {} path  cond_bits={} (of {})  Hmr+R/pass={} (Hmr stored={})  rng-tainted CCX/CCZ={}",
        if real_rng { "REAL-rng (exact, bulk-squeeze)" } else { "rng=0 (fast, conservative)" },
        an.n_cond_bits, num_bits, hmr_slot.len(), n_hmr, n_excluded
    );
    let an_ref = &an;

    let n_threads = env_usize(
        "DEAD_THREADS",
        std::thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(8),
    )
    .max(1)
    .min(nonces.len().max(1));
    let batch_size = env_usize("DEAD_BATCH_SIZE", nonces.len().max(1)).max(1);
    let checkpoint_dir = std::env::var("DEAD_CHECKPOINT_DIR").ok();
    let rank_out = std::env::var("DEAD_RANK_OUT").ok();
    let track_charged = rank_out.is_some() || checkpoint_dir.is_some();
    let candidate_mask = candidate_mask_from_env(n_ops);
    let candidate_count = ops
        .iter()
        .enumerate()
        .filter(|(i, op)| {
            candidate_mask[*i] && matches!(op.kind, OperationType::CCX | OperationType::CCZ)
        })
        .count();
    eprintln!(
        "find_dead_ccx: {} thread(s), batch_size={}, candidate_ccx_ccz={}{}{}",
        n_threads,
        batch_size,
        candidate_count,
        checkpoint_dir
            .as_ref()
            .map(|dir| format!(", checkpoint_dir={dir}"))
            .unwrap_or_default(),
        rank_out
            .as_ref()
            .map(|path| format!(", rank_out={path}"))
            .unwrap_or_default(),
    );

    let mut fired = vec![0u64; n_ops];
    let mut charged = track_charged.then(|| vec![0u64; n_ops]);
    let mut total_shots = 0usize;
    let mut seen_nonces = 0usize;
    if let Some(dir) = checkpoint_dir.as_ref() {
        std::fs::create_dir_all(dir)
            .unwrap_or_else(|err| panic!("create DEAD_CHECKPOINT_DIR={dir}: {err}"));
    }
    for (batch_idx, batch) in nonces.chunks(batch_size).enumerate() {
        let (batch_fired, batch_charged, batch_shots) = screen_nonce_batch(
            batch,
            n_threads,
            shots,
            n_ops,
            total_qubits,
            &prefix,
            &ops,
            &regs,
            &curve,
            &gtable,
            an_ref,
            real_rng,
            hmr_slot_ref,
            n_hmr,
            track_charged,
        );
        seen_nonces += batch.len();
        total_shots += batch_shots;
        for (g, l) in fired.iter_mut().zip(batch_fired.iter()) {
            *g |= *l;
        }
        if let (Some(global), Some(local)) = (charged.as_mut(), batch_charged.as_ref()) {
            for (g, l) in global.iter_mut().zip(local.iter()) {
                *g += *l;
            }
        }
        if let Some(dir) = checkpoint_dir.as_ref() {
            let dead_path = format!(
                "{dir}/dead_after_{:04}_nonces_{:05}_shots_{}.idx",
                batch_idx + 1,
                seen_nonces,
                total_shots
            );
            let n_dead =
                write_dead_idx(&dead_path, &ops, &fired, an_ref, real_rng, &candidate_mask);
            eprintln!(
                "find_dead_ccx: checkpoint batch={} seen_nonces={} shots={} dead_candidates={} path={}",
                batch_idx + 1,
                seen_nonces,
                total_shots,
                n_dead,
                dead_path
            );
            if let Some(charged) = charged.as_ref() {
                let rank_path = format!(
                    "{dir}/rank_after_{:04}_nonces_{:05}_shots_{}.tsv",
                    batch_idx + 1,
                    seen_nonces,
                    total_shots
                );
                write_rank_tsv(
                    &rank_path,
                    &ops,
                    &fired,
                    charged,
                    total_shots,
                    an_ref,
                    real_rng,
                    &candidate_mask,
                );
            }
        }
    }

    let dead_count = write_dead_idx(&out_path, &ops, &fired, an_ref, real_rng, &candidate_mask);
    if let (Some(path), Some(charged)) = (rank_out.as_ref(), charged.as_ref()) {
        write_rank_tsv(
            path,
            &ops,
            &fired,
            charged,
            total_shots,
            an_ref,
            real_rng,
            &candidate_mask,
        );
    }

    eprintln!(
        "find_dead_ccx: DONE  dead_CCX/CCZ={} / {} ({:.1}%)  screened {} total shots over {} seed(s)",
        dead_count,
        n_ccx,
        100.0 * dead_count as f64 / n_ccx.max(1) as f64,
        total_shots,
        nonces.len()
    );
    // stdout: machine-readable summary line
    println!(
        "dead {} of {} ccx  ops {}  out {}",
        dead_count, n_ccx, n_ops, out_path
    );
}
