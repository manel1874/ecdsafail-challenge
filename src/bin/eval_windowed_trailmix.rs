//! Sidecar evaluator for the QROM-lifted TrailMix windowed point-add circuit.
//!
//! This does not use the four-register challenge ABI. It validates the stronger
//! windowed contract:
//!
//! ```text
//! |a>|R> -> |a>|R + P_a>
//! ```
//!
//! The current builder intentionally supports small `k` values only. This is the
//! compatibility/fuzzing harness Craig suggested; `k = 16` should use an optimized
//! lookup implementation or an explicit Andre-style lookup cost model.

use alloy_primitives::U256;
use quantum_ecc::circuit::{analyze_ops, Op, QubitId, QubitOrBit};
use quantum_ecc::point_add::trailmix_ludicrous::windowed::{
    build_windowed_trailmix_ops, secp256k1_window_points,
};
use quantum_ecc::sim::Simulator;
use quantum_ecc::weierstrass_elliptic_curve::WeierstrassEllipticCurve;
use sha3::{
    digest::{ExtendableOutput, Update, XofReader},
    Shake256,
};

const BATCH: usize = 64;
const DEFAULT_CHALLENGE_SHOTS: usize = 9024;

#[derive(Clone, Copy)]
struct Case {
    addr: usize,
    target: (U256, U256),
    expected: (U256, U256),
}

#[derive(Clone, Copy, Debug, Default)]
struct EvalSummary {
    tested_cases: usize,
    stat_lanes: usize,
    failures: usize,
    phase_failures: usize,
    ancilla_failures: usize,
    toffoli_gates: u64,
    clifford_gates: u64,
}

impl EvalSummary {
    fn total_failures(&self) -> usize {
        self.failures + self.phase_failures + self.ancilla_failures
    }
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

fn has_arg(name: &str) -> bool {
    std::env::args().skip(1).any(|arg| arg == name)
}

fn has_arg_value(name: &str) -> bool {
    std::env::args()
        .skip(1)
        .any(|arg| arg == name || arg.starts_with(&format!("{name}=")))
}

fn parse_arg_usize(name: &str, default: usize) -> usize {
    let mut args = std::env::args().skip(1);
    while let Some(arg) = args.next() {
        if arg == name {
            return args
                .next()
                .unwrap_or_else(|| panic!("{name} requires a value"))
                .parse()
                .unwrap_or_else(|_| panic!("{name} expects usize"));
        }
        if let Some(rest) = arg.strip_prefix(&format!("{name}=")) {
            return rest
                .parse()
                .unwrap_or_else(|_| panic!("{name} expects usize"));
        }
    }
    default
}

fn parse_arg_u64(name: &str, default: u64) -> u64 {
    let mut args = std::env::args().skip(1);
    while let Some(arg) = args.next() {
        if arg == name {
            return args
                .next()
                .unwrap_or_else(|| panic!("{name} requires a value"))
                .parse()
                .unwrap_or_else(|_| panic!("{name} expects u64"));
        }
        if let Some(rest) = arg.strip_prefix(&format!("{name}=")) {
            return rest
                .parse()
                .unwrap_or_else(|_| panic!("{name} expects u64"));
        }
    }
    default
}

fn windowed_fiat_shamir_seed(
    ops: &[Op],
    k: usize,
    entries: usize,
    mode: &[u8],
) -> sha3::Shake256Reader {
    let mut hasher = Shake256::default();
    hasher.update(b"windowed-trailmix-fiat-shamir-v1");
    hasher.update(mode);
    hasher.update(&(k as u64).to_le_bytes());
    hasher.update(&(entries as u64).to_le_bytes());
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

fn deterministic_seed(k: usize) -> sha3::Shake256Reader {
    let mut hasher = Shake256::default();
    hasher.update(b"windowed-trailmix-eval");
    hasher.update(&(k as u64).to_le_bytes());
    hasher.finalize_xof()
}

fn point_at_infinity(p: (U256, U256)) -> bool {
    p.0.is_zero() && p.1.is_zero()
}

fn deterministic_cases(
    curve: &WeierstrassEllipticCurve,
    table: &[(U256, U256)],
    targets: usize,
) -> Vec<Case> {
    let entries = table.len();
    let mut cases = Vec::with_capacity(entries * targets);
    for t in 0..targets {
        let target = curve.mul(curve.gx, curve.gy, U256::from((1000 + t) as u64));
        for (addr, &q) in table.iter().enumerate() {
            let expected = curve.add(target.0, target.1, q.0, q.1);
            cases.push(Case {
                addr,
                target,
                expected,
            });
        }
    }
    cases
}

fn challenge_style_cases(
    curve: &WeierstrassEllipticCurve,
    table: &[(U256, U256)],
    shots: usize,
    include_zero_address: bool,
    xof: &mut impl XofReader,
) -> Vec<Case> {
    let entries = table.len();
    let mut cases = Vec::with_capacity(shots);
    while cases.len() < shots {
        let mut rb = [[0u8; 32]; 2];
        XofReader::read(xof, &mut rb[0]);
        XofReader::read(xof, &mut rb[1]);

        let target = curve.mul(curve.gx, curve.gy, U256::from_le_bytes(rb[0]));
        if point_at_infinity(target) {
            continue;
        }

        let addr_word = u64::from_le_bytes(rb[1][..8].try_into().unwrap()) as usize;
        let addr = if include_zero_address {
            addr_word & (entries - 1)
        } else {
            1 + (addr_word % (entries - 1))
        };
        let q = table[addr];

        // The affine mixed-add formula is undefined when x_R == x(P_i).
        // The challenge harness similarly rejects x-coordinate collisions.
        if target.0 == q.0 {
            continue;
        }

        let expected = curve.add(target.0, target.1, q.0, q.1);
        cases.push(Case {
            addr,
            target,
            expected,
        });
    }
    cases
}

fn set_batch_register(
    sim: &mut Simulator<impl sha3::digest::XofReader>,
    reg: &[QubitOrBit],
    vals: &[U256],
) {
    for (shot, &value) in vals.iter().enumerate() {
        sim.set_register(reg, value, shot);
    }
}

fn evaluate_cases(
    ops: &[Op],
    total_qubits: u64,
    num_bits: u64,
    regs: &[Vec<QubitOrBit>],
    cases: &[Case],
    mut xof: sha3::Shake256Reader,
    verbose_mismatches: bool,
) -> EvalSummary {
    let mut sim = Simulator::new(total_qubits as usize, num_bits as usize, &mut xof);

    let mut summary = EvalSummary {
        tested_cases: cases.len(),
        ..EvalSummary::default()
    };
    for chunk in cases.chunks(BATCH) {
        summary.stat_lanes += BATCH;
        sim.clear_for_shot();
        let mut addrs = Vec::with_capacity(chunk.len());
        let mut xs = Vec::with_capacity(chunk.len());
        let mut ys = Vec::with_capacity(chunk.len());
        for case in chunk {
            addrs.push(U256::from(case.addr as u64));
            xs.push(case.target.0);
            ys.push(case.target.1);
        }
        set_batch_register(&mut sim, &regs[0], &addrs);
        set_batch_register(&mut sim, &regs[1], &xs);
        set_batch_register(&mut sim, &regs[2], &ys);

        sim.apply_iter(ops.iter());

        let cond_mask = if chunk.len() == 64 {
            u64::MAX
        } else {
            (1u64 << chunk.len()) - 1
        };
        if sim.phase & cond_mask != 0 {
            summary.phase_failures += 1;
        }

        for (shot, case) in chunk.iter().enumerate() {
            let got_addr = sim.get_register(&regs[0], shot);
            let got_x = sim.get_register(&regs[1], shot);
            let got_y = sim.get_register(&regs[2], shot);
            if got_addr != U256::from(case.addr as u64)
                || got_x != case.expected.0
                || got_y != case.expected.1
            {
                summary.failures += 1;
                if verbose_mismatches {
                    eprintln!(
                        "mismatch shot={shot} addr={}: got addr={got_addr:#x} ({got_x:#x}, {got_y:#x}) expected ({:#x}, {:#x})",
                        case.addr, case.expected.0, case.expected.1
                    );
                }
            }
        }

        for reg in regs {
            for item in reg {
                if let QubitOrBit::Qubit(q) = *item {
                    *sim.qubit_mut(q) = 0;
                }
            }
        }
        for q in 0..total_qubits {
            if sim.qubit(QubitId(q)) & cond_mask != 0 {
                summary.ancilla_failures += 1;
                if verbose_mismatches {
                    eprintln!(
                        "ancilla garbage: q{q} = {:#018x}",
                        sim.qubit(QubitId(q)) & cond_mask
                    );
                }
                break;
            }
        }
    }

    summary.toffoli_gates = sim.stats.toffoli_gates;
    summary.clifford_gates = sim.stats.clifford_gates;
    summary
}

fn main() {
    let k = parse_arg_usize("--k", 1);
    let targets = parse_arg_usize("--targets", 2);
    let mut challenge_shots = parse_arg_usize("--challenge-shots", 0);
    let nonce_arg = has_arg_value("--nonce");
    let nonce = parse_arg_u64("--nonce", 0);
    let hunt_nonces = parse_arg_usize("--hunt-nonces", 0);
    let hunt_start = parse_arg_u64("--hunt-start", 0);
    let hunt_shots = parse_arg_usize(
        "--hunt-shots",
        if challenge_shots == 0 { BATCH } else { challenge_shots },
    );
    if has_arg("--challenge-style") && challenge_shots == 0 {
        challenge_shots = DEFAULT_CHALLENGE_SHOTS;
    }
    if hunt_nonces != 0 {
        challenge_shots = hunt_shots;
    }
    let include_zero_address = has_arg("--include-zero-address");
    let hunt_stop_on_clean = has_arg("--hunt-stop-on-clean");
    assert!(k <= 6, "small-k evaluator supports k <= 6");
    if challenge_shots != 0 {
        assert_eq!(
            challenge_shots % BATCH,
            0,
            "challenge-style shot count must be a multiple of {BATCH}"
        );
    }
    if hunt_nonces != 0 {
        std::env::set_var("DIALOG_TAIL_NONCE", hunt_start.to_string());
    } else if nonce_arg {
        std::env::set_var("DIALOG_TAIL_NONCE", nonce.to_string());
    }

    println!("=== windowed TrailMix compatibility eval ===");
    println!("  address bits : {k}");
    if challenge_shots == 0 {
        println!("  mode         : deterministic exhaustive addresses");
        println!("  target count : {targets}");
    } else {
        println!(
            "  mode         : {}",
            if hunt_nonces == 0 {
                "challenge-style random"
            } else {
                "challenge-style nonce hunt"
            }
        );
        println!("  random shots : {challenge_shots}");
        println!(
            "  zero address : {}",
            if include_zero_address {
                "included"
            } else {
                "excluded, matching the challenge infinity-offset skip"
            }
        );
        if hunt_nonces != 0 {
            println!("  hunt start   : {hunt_start}");
            println!("  hunt nonces  : {hunt_nonces}");
        }
    }

    let (ops, layout) = build_windowed_trailmix_ops(k);
    let (total_qubits, num_bits, _num_regs, regs) = analyze_ops(ops.iter());
    let entries = layout.table_entries;
    assert_eq!(entries, 1usize << k);
    assert_eq!(regs.len(), 3);

    println!("  emitted ops  : {}", ops.len());
    println!("  qubits       : {total_qubits}");
    println!("  bits         : {num_bits}");
    println!("  registers    : {}", regs.len());
    println!("  lookup path  : {:?}", layout.lookup_kind);
    println!(
        "  tail nonce   : {}",
        std::env::var("DIALOG_TAIL_NONCE").unwrap_or_else(|_| "<none>".to_string())
    );

    let curve = secp256k1();
    let table = secp256k1_window_points(entries);

    if hunt_nonces != 0 {
        let mut best: Option<(u64, EvalSummary)> = None;
        for offset in 0..hunt_nonces {
            let candidate_nonce = hunt_start + offset as u64;
            std::env::set_var("DIALOG_TAIL_NONCE", candidate_nonce.to_string());
            let (candidate_ops, candidate_layout) = build_windowed_trailmix_ops(k);
            let (candidate_total_qubits, candidate_num_bits, _candidate_num_regs, candidate_regs) =
                analyze_ops(candidate_ops.iter());
            assert_eq!(candidate_layout.table_entries, entries);
            let mut xof =
                windowed_fiat_shamir_seed(&candidate_ops, k, entries, b"challenge-style");
            let cases = challenge_style_cases(
                &curve,
                &table,
                challenge_shots,
                include_zero_address,
                &mut xof,
            );
            let summary = evaluate_cases(
                &candidate_ops,
                candidate_total_qubits,
                candidate_num_bits,
                &candidate_regs,
                &cases,
                xof,
                false,
            );
            println!(
                "HUNT nonce={} shots={} mismatches={} phase_batches={} ancilla_batches={} total={}",
                candidate_nonce,
                summary.tested_cases,
                summary.failures,
                summary.phase_failures,
                summary.ancilla_failures,
                summary.total_failures()
            );
            if best
                .as_ref()
                .map(|(_, best_summary)| summary.total_failures() < best_summary.total_failures())
                .unwrap_or(true)
            {
                best = Some((candidate_nonce, summary));
            }
            if hunt_stop_on_clean && summary.total_failures() == 0 {
                break;
            }
        }

        if let Some((best_nonce, summary)) = best {
            println!(
                "HUNT_BEST nonce={} shots={} mismatches={} phase_batches={} ancilla_batches={} total={}",
                best_nonce,
                summary.tested_cases,
                summary.failures,
                summary.phase_failures,
                summary.ancilla_failures,
                summary.total_failures()
            );
        }
        return;
    }

    let mut xof = if challenge_shots == 0 {
        deterministic_seed(k)
    } else {
        windowed_fiat_shamir_seed(&ops, k, entries, b"challenge-style")
    };
    let cases = if challenge_shots == 0 {
        deterministic_cases(&curve, &table, targets)
    } else {
        challenge_style_cases(
            &curve,
            &table,
            challenge_shots,
            include_zero_address,
            &mut xof,
        )
    };
    let summary = evaluate_cases(&ops, total_qubits, num_bits, &regs, &cases, xof, true);

    let shots = summary.stat_lanes.max(1) as f64;
    println!("  tested shots : {}", cases.len());
    println!(
        "  avg Toffoli  : {:.3}",
        summary.toffoli_gates as f64 / shots
    );
    println!(
        "  avg Clifford : {:.3}",
        summary.clifford_gates as f64 / shots
    );
    println!("  mismatches   : {}", summary.failures);
    println!("  phase fails  : {}", summary.phase_failures);
    println!("  ancilla fails: {}", summary.ancilla_failures);

    let core_t = summary.toffoli_gates as f64 / shots;
    let q_full = total_qubits + 16;
    let t_full = core_t + 3.0 * ((1u64 << 16) as f64);
    let t_shor = 28.0 * t_full;
    println!("  w=16 analytic:");
    println!("    Q_full = Q_core + 16 = {q_full}");
    println!("    T_full = T_core + 3*2^16 = {t_full:.3}");
    println!("    T_Shor ~= 28*T_full = {t_shor:.3}");

    if summary.total_failures() != 0 {
        std::process::exit(1);
    }
    println!("=== windowed TrailMix OK ===");
}
