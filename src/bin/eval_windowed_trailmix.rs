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
use quantum_ecc::circuit::{analyze_ops, QubitId, QubitOrBit};
use quantum_ecc::point_add::trailmix_ludicrous::windowed::{
    build_windowed_trailmix_ops, secp256k1_window_points,
};
use quantum_ecc::sim::Simulator;
use quantum_ecc::weierstrass_elliptic_curve::WeierstrassEllipticCurve;
use sha3::{
    digest::{ExtendableOutput, Update},
    Shake256,
};

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

fn set_batch_register(
    sim: &mut Simulator<impl sha3::digest::XofReader>,
    reg: &[QubitOrBit],
    vals: &[U256],
) {
    for (shot, &value) in vals.iter().enumerate() {
        sim.set_register(reg, value, shot);
    }
}

fn main() {
    let k = parse_arg_usize("--k", 1);
    let targets = parse_arg_usize("--targets", 2);
    assert!(k <= 6, "small-k evaluator supports k <= 6");

    println!("=== windowed TrailMix compatibility eval ===");
    println!("  address bits : {k}");
    println!("  target count : {targets}");

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

    let curve = secp256k1();
    let table = secp256k1_window_points(entries);

    let mut cases = Vec::with_capacity(entries * targets);
    for t in 0..targets {
        let target = curve.mul(curve.gx, curve.gy, U256::from((1000 + t) as u64));
        for addr in 0..entries {
            let q = table[addr];
            let expected = curve.add(target.0, target.1, q.0, q.1);
            cases.push((addr, target, expected));
        }
    }

    let mut hasher = Shake256::default();
    hasher.update(b"windowed-trailmix-eval");
    hasher.update(&(k as u64).to_le_bytes());
    let mut xof = hasher.finalize_xof();
    let mut sim = Simulator::new(total_qubits as usize, num_bits as usize, &mut xof);

    const BATCH: usize = 64;
    let mut failures = 0usize;
    let mut phase_failures = 0usize;
    let mut ancilla_failures = 0usize;
    let mut stat_lanes = 0usize;
    for chunk in cases.chunks(BATCH) {
        stat_lanes += BATCH;
        sim.clear_for_shot();
        let mut addrs = Vec::with_capacity(chunk.len());
        let mut xs = Vec::with_capacity(chunk.len());
        let mut ys = Vec::with_capacity(chunk.len());
        for &(addr, target, _) in chunk {
            addrs.push(U256::from(addr as u64));
            xs.push(target.0);
            ys.push(target.1);
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
            phase_failures += 1;
        }

        for (shot, &(addr, _, expected)) in chunk.iter().enumerate() {
            let got_addr = sim.get_register(&regs[0], shot);
            let got_x = sim.get_register(&regs[1], shot);
            let got_y = sim.get_register(&regs[2], shot);
            if got_addr != U256::from(addr as u64) || got_x != expected.0 || got_y != expected.1 {
                failures += 1;
                eprintln!(
                    "mismatch shot={shot} addr={addr}: got addr={got_addr:#x} ({got_x:#x}, {got_y:#x}) expected ({:#x}, {:#x})",
                    expected.0, expected.1
                );
            }
        }

        for reg in &regs {
            for item in reg {
                if let QubitOrBit::Qubit(q) = *item {
                    *sim.qubit_mut(q) = 0;
                }
            }
        }
        for q in 0..total_qubits {
            if sim.qubit(QubitId(q)) & cond_mask != 0 {
                ancilla_failures += 1;
                eprintln!(
                    "ancilla garbage: q{q} = {:#018x}",
                    sim.qubit(QubitId(q)) & cond_mask
                );
                break;
            }
        }
    }

    let shots = stat_lanes.max(1) as f64;
    println!("  tested shots : {}", cases.len());
    println!(
        "  avg Toffoli  : {:.3}",
        sim.stats.toffoli_gates as f64 / shots
    );
    println!(
        "  avg Clifford : {:.3}",
        sim.stats.clifford_gates as f64 / shots
    );
    println!("  mismatches   : {failures}");
    println!("  phase fails  : {phase_failures}");
    println!("  ancilla fails: {ancilla_failures}");

    let core_t = sim.stats.toffoli_gates as f64 / shots;
    let q_full = total_qubits + 16;
    let t_full = core_t + 3.0 * ((1u64 << 16) as f64);
    let t_shor = 28.0 * t_full;
    println!("  w=16 analytic:");
    println!("    Q_full = Q_core + 16 = {q_full}");
    println!("    T_full = T_core + 3*2^16 = {t_full:.3}");
    println!("    T_Shor ~= 28*T_full = {t_shor:.3}");

    if failures != 0 || phase_failures != 0 || ancilla_failures != 0 {
        std::process::exit(1);
    }
    println!("=== windowed TrailMix OK ===");
}
