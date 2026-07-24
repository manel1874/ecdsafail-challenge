//! Build the submitted point-addition circuit and export a six-phase resource profile.
//!
//! The output is intentionally binned: a row represents one equal-duration slice
//! within a phase, while `index` remains an operation index in the final circuit.
//! This keeps the CSV compact enough for pgfplots without hiding the very different
//! raw operation spans of the six algorithmic phases.

use alloy_primitives::U256;
use quantum_ecc::circuit::{analyze_ops, Op, OperationType, QubitOrBit, NO_BIT};
use quantum_ecc::point_add::{
    self, CircuitResourceTrace, RESOURCE_PHASE_COUNT, RESOURCE_PHASE_NAMES,
};
use quantum_ecc::sim::Simulator;
use quantum_ecc::weierstrass_elliptic_curve::WeierstrassEllipticCurve;
use sha3::{
    digest::{ExtendableOutput, Update, XofReader},
    Shake256,
};
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

const SHOTS_PER_TRIAL: usize = 64;
const HARNESS_TESTS: usize = 9024;
const DEFAULT_TRIALS: usize = HARNESS_TESTS / SHOTS_PER_TRIAL;
const DEFAULT_BINS_PER_PHASE: usize = 96;
const DEFAULT_DENSE_BINS_PER_PHASE: usize = 3 * DEFAULT_BINS_PER_PHASE;
const DEFAULT_OUTPUT: &str = "figures/circuit_profile.csv";

#[derive(Debug)]
struct Args {
    output: PathBuf,
    trials: usize,
    bins_per_phase: usize,
    dense_bins_per_phase: usize,
}

fn parse_positive(flag: &str, value: Option<String>) -> Result<usize, String> {
    let raw = value.ok_or_else(|| format!("{flag} requires a value"))?;
    let parsed = raw
        .parse::<usize>()
        .map_err(|_| format!("{flag} expects a positive integer, got {raw:?}"))?;
    if parsed == 0 {
        return Err(format!("{flag} must be greater than zero"));
    }
    Ok(parsed)
}

fn parse_args() -> Result<Args, String> {
    let mut output = PathBuf::from(DEFAULT_OUTPUT);
    let mut trials = DEFAULT_TRIALS;
    let mut bins_per_phase = DEFAULT_BINS_PER_PHASE;
    let mut dense_bins_per_phase = DEFAULT_DENSE_BINS_PER_PHASE;
    let mut args = std::env::args().skip(1);

    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--output" => {
                output = PathBuf::from(
                    args.next()
                        .ok_or_else(|| "--output requires a path".to_owned())?,
                );
            }
            "--trials" => trials = parse_positive("--trials", args.next())?,
            "--bins-per-phase" => bins_per_phase = parse_positive("--bins-per-phase", args.next())?,
            "--dense-bins-per-phase" => {
                dense_bins_per_phase = parse_positive("--dense-bins-per-phase", args.next())?
            }
            "--help" | "-h" => {
                println!(
                    "Usage: profile_circuit [--output PATH] [--trials N] \
                     [--bins-per-phase N] [--dense-bins-per-phase N]\n\
                     Defaults: --output {DEFAULT_OUTPUT} --trials {DEFAULT_TRIALS} \
                     --bins-per-phase {DEFAULT_BINS_PER_PHASE} \
                     --dense-bins-per-phase {DEFAULT_DENSE_BINS_PER_PHASE}\n\
                     Dense bins apply to inverse, square, and forward-multiply.\n\
                     One trial is one 64-shot batch from the challenge's Fiat--Shamir stream."
                );
                std::process::exit(0);
            }
            _ if arg.starts_with("--output=") => {
                output = PathBuf::from(arg.trim_start_matches("--output="));
            }
            _ if arg.starts_with("--trials=") => {
                trials = parse_positive(
                    "--trials",
                    Some(arg.trim_start_matches("--trials=").to_owned()),
                )?;
            }
            _ if arg.starts_with("--bins-per-phase=") => {
                bins_per_phase = parse_positive(
                    "--bins-per-phase",
                    Some(arg.trim_start_matches("--bins-per-phase=").to_owned()),
                )?;
            }
            _ if arg.starts_with("--dense-bins-per-phase=") => {
                dense_bins_per_phase = parse_positive(
                    "--dense-bins-per-phase",
                    Some(arg.trim_start_matches("--dense-bins-per-phase=").to_owned()),
                )?;
            }
            _ => return Err(format!("unknown argument {arg:?}; use --help")),
        }
    }

    Ok(Args {
        output,
        trials,
        bins_per_phase,
        dense_bins_per_phase,
    })
}

fn phase_bins(base: usize, dense: usize) -> [usize; RESOURCE_PHASE_COUNT + 1] {
    let mut bins = [base; RESOURCE_PHASE_COUNT + 1];
    bins[0] = 0;
    for phase in [2usize, 4, 5] {
        bins[phase] = dense;
    }
    bins
}

fn phase_bin_offsets(
    bins: &[usize; RESOURCE_PHASE_COUNT + 1],
) -> [usize; RESOURCE_PHASE_COUNT + 2] {
    let mut offsets = [0usize; RESOURCE_PHASE_COUNT + 2];
    for phase in 1..=RESOURCE_PHASE_COUNT {
        offsets[phase + 1] = offsets[phase] + bins[phase];
    }
    offsets
}

fn secp256k1() -> WeierstrassEllipticCurve {
    WeierstrassEllipticCurve {
        modulus: U256::from_str_radix(
            "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F",
            16,
        )
        .unwrap(),
        a: U256::ZERO,
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

fn fiat_shamir_seed(ops: &[Op]) -> sha3::Shake256Reader {
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

type AffinePoint = (U256, U256);

fn make_harness_inputs(
    curve: &WeierstrassEllipticCurve,
    xof: &mut sha3::Shake256Reader,
) -> (Vec<AffinePoint>, Vec<AffinePoint>, Vec<AffinePoint>) {
    let mut targets = Vec::with_capacity(HARNESS_TESTS);
    let mut offsets = Vec::with_capacity(HARNESS_TESTS);
    let mut expected = Vec::with_capacity(HARNESS_TESTS);

    for _ in 0..HARNESS_TESTS {
        let mut random = [[0u8; 32]; 2];
        XofReader::read(xof, &mut random[0]);
        XofReader::read(xof, &mut random[1]);
        let target = curve.mul(curve.gx, curve.gy, U256::from_le_bytes(random[0]));
        let offset = curve.mul(curve.gx, curve.gy, U256::from_le_bytes(random[1]));
        if target.0 == offset.0
            || (target.0.is_zero() && target.1.is_zero())
            || (offset.0.is_zero() && offset.1.is_zero())
        {
            continue;
        }
        let result = curve.add(target.0, target.1, offset.0, offset.1);
        targets.push(target);
        offsets.push(offset);
        expected.push(result);
    }
    (targets, offsets, expected)
}

fn evaluator_toffoli_total(
    ops: &[Op],
    registers: &[Vec<QubitOrBit>],
    num_qubits: u64,
    num_bits: u64,
    trials: usize,
) -> u64 {
    let curve = secp256k1();
    let mut xof = fiat_shamir_seed(ops);
    let (targets, offsets, _) = make_harness_inputs(&curve, &mut xof);
    let mut simulator = Simulator::new(num_qubits as usize, num_bits as usize, &mut xof);

    for trial in 0..trials {
        simulator.clear_for_shot();
        let batch_start = trial * SHOTS_PER_TRIAL;
        for shot in 0..SHOTS_PER_TRIAL {
            let input = batch_start + shot;
            simulator.set_register(&registers[0], targets[input].0, shot);
            simulator.set_register(&registers[1], targets[input].1, shot);
            simulator.set_register(&registers[2], offsets[input].0, shot);
            simulator.set_register(&registers[3], offsets[input].1, shot);
        }
        simulator.apply_iter(ops.iter());
    }
    simulator.stats.toffoli_gates
}

fn phase_bin(
    phase: u8,
    seen: &mut [usize; RESOURCE_PHASE_COUNT + 1],
    phase_counts: &[usize; RESOURCE_PHASE_COUNT + 1],
    bins_per_phase: &[usize; RESOURCE_PHASE_COUNT + 1],
    bin_offsets: &[usize; RESOURCE_PHASE_COUNT + 2],
) -> Option<usize> {
    let phase_index = usize::from(phase);
    if phase_index == 0 || phase_index > RESOURCE_PHASE_COUNT {
        return None;
    }
    let ordinal = seen[phase_index];
    seen[phase_index] += 1;
    let phase_bins = bins_per_phase[phase_index];
    let local_bin = (ordinal * phase_bins / phase_counts[phase_index]).min(phase_bins - 1);
    Some(bin_offsets[phase_index] + local_bin)
}

fn apply_profiled(
    sim: &mut Simulator<'_, sha3::Shake256Reader>,
    ops: &[Op],
    phases: &[u8],
    phase_counts: &[usize; RESOURCE_PHASE_COUNT + 1],
    bins_per_phase: &[usize; RESOURCE_PHASE_COUNT + 1],
    bin_offsets: &[usize; RESOURCE_PHASE_COUNT + 2],
    executed_toffoli: &mut [u64],
    executed_ccx_total: &mut u64,
    executed_ccz_total: &mut u64,
) {
    let mut condition_stack = Vec::new();
    let mut current_base_condition = u64::MAX;
    let mut seen = [0usize; RESOURCE_PHASE_COUNT + 1];

    for (op, &phase) in ops.iter().zip(phases) {
        let bin = phase_bin(phase, &mut seen, phase_counts, bins_per_phase, bin_offsets);
        let mut condition = current_base_condition;
        if op.c_condition != NO_BIT {
            condition &= sim.bit(op.c_condition);
        }
        if matches!(op.kind, OperationType::CCX | OperationType::CCZ) {
            let executed_shots = u64::from(condition.count_ones());
            match op.kind {
                OperationType::CCX => *executed_ccx_total += executed_shots,
                OperationType::CCZ => *executed_ccz_total += executed_shots,
                _ => unreachable!(),
            }
            if let Some(bin) = bin {
                executed_toffoli[bin] += executed_shots;
            }
        }

        match op.kind {
            OperationType::CCX => {
                let value = condition & sim.qubit(op.q_control1) & sim.qubit(op.q_control2);
                *sim.qubit_mut(op.q_target) ^= value;
            }
            OperationType::CX => {
                let value = condition & sim.qubit(op.q_control1);
                *sim.qubit_mut(op.q_target) ^= value;
            }
            OperationType::Swap => {
                let mut control = sim.qubit(op.q_control1);
                let mut target = sim.qubit(op.q_target);
                control ^= target;
                target ^= condition & control;
                control ^= target;
                *sim.qubit_mut(op.q_control1) = control;
                *sim.qubit_mut(op.q_target) = target;
            }
            OperationType::X => *sim.qubit_mut(op.q_target) ^= condition,
            OperationType::CCZ => {
                sim.phase ^= condition
                    & sim.qubit(op.q_target)
                    & sim.qubit(op.q_control1)
                    & sim.qubit(op.q_control2);
            }
            OperationType::CZ => {
                sim.phase ^= condition & sim.qubit(op.q_target) & sim.qubit(op.q_control1);
            }
            OperationType::Z => sim.phase ^= condition & sim.qubit(op.q_target),
            OperationType::Neg => sim.phase ^= condition,
            OperationType::Hmr => {
                let mut bytes = [0u8; 8];
                XofReader::read(sim.xof, &mut bytes);
                let random = u64::from_le_bytes(bytes);
                *sim.bit_mut(op.c_target) &= !condition;
                *sim.bit_mut(op.c_target) ^= random & condition;
                sim.phase ^= sim.qubit(op.q_target) & random & condition;
                *sim.qubit_mut(op.q_target) &= !condition;
            }
            OperationType::R => {
                let mut bytes = [0u8; 8];
                XofReader::read(sim.xof, &mut bytes);
                let random = u64::from_le_bytes(bytes);
                sim.phase ^= sim.qubit(op.q_target) & random & condition;
                *sim.qubit_mut(op.q_target) &= !condition;
            }
            OperationType::BitInvert => *sim.bit_mut(op.c_target) ^= condition,
            OperationType::BitStore0 => *sim.bit_mut(op.c_target) &= !condition,
            OperationType::BitStore1 => *sim.bit_mut(op.c_target) |= condition,
            OperationType::PushCondition => {
                condition_stack.push(current_base_condition);
                current_base_condition &= sim.bit(op.c_condition);
            }
            OperationType::PopCondition => {
                current_base_condition = condition_stack
                    .pop()
                    .expect("POP_CONDITION without PUSH_CONDITION");
            }
            OperationType::AppendToRegister
            | OperationType::Register
            | OperationType::DebugPrint => {}
        }
    }
    assert!(condition_stack.is_empty(), "unclosed condition stack");
}

fn verify_trial(
    sim: &mut Simulator<'_, sha3::Shake256Reader>,
    registers: &[Vec<QubitOrBit>],
    expected: &[AffinePoint],
) -> Result<(), String> {
    for (shot, &(want_x, want_y)) in expected.iter().enumerate() {
        let got_x = sim.get_register(&registers[0], shot);
        let got_y = sim.get_register(&registers[1], shot);
        if got_x != want_x || got_y != want_y {
            return Err(format!(
                "trial verification failed at shot {shot}: got ({got_x:#x}, {got_y:#x}), \
                 expected ({want_x:#x}, {want_y:#x})"
            ));
        }
    }
    if sim.phase != 0 {
        return Err(format!(
            "trial left non-zero global phase mask {:#018x}",
            sim.phase
        ));
    }

    for register in registers {
        for item in register {
            if let QubitOrBit::Qubit(qubit) = item {
                *sim.qubit_mut(*qubit) = 0;
            }
        }
    }
    if let Some((index, value)) = sim
        .qubits
        .iter()
        .copied()
        .enumerate()
        .find(|(_, value)| *value != 0)
    {
        return Err(format!(
            "trial left ancilla garbage on qubit {index}: {value:#018x}"
        ));
    }
    Ok(())
}

fn phase_counts(phases: &[u8]) -> Result<[usize; RESOURCE_PHASE_COUNT + 1], String> {
    let mut counts = [0usize; RESOURCE_PHASE_COUNT + 1];
    for &phase in phases {
        let index = usize::from(phase);
        if index > RESOURCE_PHASE_COUNT {
            return Err(format!("invalid resource phase id {phase}"));
        }
        counts[index] += 1;
    }
    for phase in 1..=RESOURCE_PHASE_COUNT {
        if counts[phase] == 0 {
            return Err(format!(
                "resource phase {} has no final operations",
                RESOURCE_PHASE_NAMES[phase]
            ));
        }
    }
    Ok(counts)
}

fn final_bin_geometry(
    phases: &[u8],
    counts: &[usize; RESOURCE_PHASE_COUNT + 1],
    bins_per_phase: &[usize; RESOURCE_PHASE_COUNT + 1],
    bin_offsets: &[usize; RESOURCE_PHASE_COUNT + 2],
) -> (Vec<u64>, Vec<u64>) {
    let total_bins = bin_offsets[RESOURCE_PHASE_COUNT + 1];
    let mut index_sum = vec![0u64; total_bins];
    let mut op_count = vec![0u64; total_bins];
    let mut seen = [0usize; RESOURCE_PHASE_COUNT + 1];

    for (index, &phase) in phases.iter().enumerate() {
        if let Some(bin) = phase_bin(phase, &mut seen, counts, bins_per_phase, bin_offsets) {
            index_sum[bin] += index as u64;
            op_count[bin] += 1;
        }
    }
    (index_sum, op_count)
}

fn source_phase_ranges(trace: &CircuitResourceTrace) -> Result<[(usize, usize); 7], String> {
    let mut ranges = [(0usize, 0usize); 7];
    for (position, &(start, phase)) in trace.source_phase_transitions.iter().enumerate() {
        let phase = usize::from(phase);
        if phase == 0 || phase > RESOURCE_PHASE_COUNT {
            continue;
        }
        let end = trace
            .source_phase_transitions
            .get(position + 1)
            .map(|transition| transition.0)
            .ok_or_else(|| {
                format!(
                    "phase {} has no closing transition",
                    RESOURCE_PHASE_NAMES[phase]
                )
            })?;
        ranges[phase] = (start, end);
    }
    for phase in 1..=RESOURCE_PHASE_COUNT {
        if ranges[phase].0 >= ranges[phase].1 {
            return Err(format!(
                "invalid source range for phase {}: {:?}",
                RESOURCE_PHASE_NAMES[phase], ranges[phase]
            ));
        }
    }
    Ok(ranges)
}

fn active_at(timeline: &[(usize, u32)], index: usize) -> (usize, u32) {
    let event = timeline.partition_point(|&(event_index, _)| event_index <= index);
    if event == 0 {
        (0, 0)
    } else {
        (event, timeline[event - 1].1)
    }
}

fn active_bins(
    trace: &CircuitResourceTrace,
    bins_per_phase: &[usize; RESOURCE_PHASE_COUNT + 1],
) -> Result<Vec<u32>, String> {
    let ranges = source_phase_ranges(trace)?;
    let mut values = Vec::with_capacity(bins_per_phase.iter().sum());
    for phase in 1..=RESOURCE_PHASE_COUNT {
        let (start, end) = ranges[phase];
        let length = end - start;
        let phase_bins = bins_per_phase[phase];
        for bin in 0..phase_bins {
            let midpoint = start + length * (2 * bin + 1) / (2 * phase_bins);
            values.push(active_at(&trace.active_timeline, midpoint.min(end - 1)).1);
        }
    }
    Ok(values)
}

fn write_csv(
    path: &Path,
    trials: usize,
    bins_per_phase: &[usize; RESOURCE_PHASE_COUNT + 1],
    bin_offsets: &[usize; RESOURCE_PHASE_COUNT + 2],
    index_sum: &[u64],
    op_count: &[u64],
    active: &[u32],
    executed_toffoli: &[u64],
) -> Result<(), String> {
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent)
            .map_err(|error| format!("create {}: {error}", parent.display()))?;
    }
    let file = File::create(path).map_err(|error| format!("create {}: {error}", path.display()))?;
    let mut writer = BufWriter::new(file);
    writeln!(
        writer,
        "index,active_qubit_count_per_index,avg_toffoli_gates_per_run_in_bin,phase_name"
    )
    .map_err(|error| error.to_string())?;

    let shots = (trials * SHOTS_PER_TRIAL) as f64;
    for phase in 1..=RESOURCE_PHASE_COUNT {
        for local_bin in 0..bins_per_phase[phase] {
            let bin = bin_offsets[phase] + local_bin;
            if op_count[bin] == 0 {
                return Err(format!("empty output bin {bin}"));
            }
            let index = index_sum[bin] / op_count[bin];
            let mean_toffoli_per_run = executed_toffoli[bin] as f64 / shots;
            writeln!(
                writer,
                "{index},{},{:.6},{}",
                active[bin], mean_toffoli_per_run, RESOURCE_PHASE_NAMES[phase]
            )
            .map_err(|error| error.to_string())?;
        }
    }
    writer.flush().map_err(|error| error.to_string())
}

fn run(args: Args) -> Result<(), String> {
    std::env::set_var("PROFILE_CIRCUIT_RESOURCES", "1");
    println!("building final challenge circuit with resource tracing enabled...");
    let ops = point_add::build();
    let trace = point_add::take_last_resource_trace()
        .ok_or_else(|| "resource trace was not captured".to_owned())?;
    if trace.final_op_phases.len() != ops.len() {
        return Err(format!(
            "final phase trace length {} does not match operation count {}",
            trace.final_op_phases.len(),
            ops.len()
        ));
    }

    let counts = phase_counts(&trace.final_op_phases)?;
    let bins_per_phase = phase_bins(args.bins_per_phase, args.dense_bins_per_phase);
    let bin_offsets = phase_bin_offsets(&bins_per_phase);
    let total_bins = bin_offsets[RESOURCE_PHASE_COUNT + 1];
    let (num_qubits, num_bits, _, registers) = analyze_ops(ops.iter());
    let peak_active = trace
        .active_timeline
        .iter()
        .map(|&(_, active)| active)
        .max()
        .unwrap_or(0);
    if u64::from(peak_active) != num_qubits {
        return Err(format!(
            "allocator peak {peak_active} does not match analyzed circuit width {num_qubits}"
        ));
    }
    if registers.len() < 4 {
        return Err(format!(
            "expected four circuit registers, found {}",
            registers.len()
        ));
    }
    let (index_sum, op_count) = final_bin_geometry(
        &trace.final_op_phases,
        &counts,
        &bins_per_phase,
        &bin_offsets,
    );
    let active = active_bins(&trace, &bins_per_phase)?;
    let mut executed_toffoli = vec![0u64; total_bins];
    let mut executed_ccx_total = 0u64;
    let mut executed_ccz_total = 0u64;
    let curve = secp256k1();
    let mut xof = fiat_shamir_seed(&ops);
    let (targets, offsets, expected) = make_harness_inputs(&curve, &mut xof);
    let available_trials = targets.len() / SHOTS_PER_TRIAL;
    if args.trials > available_trials {
        return Err(format!(
            "--trials {} exceeds the {} complete 64-shot batches in the harness stream",
            args.trials, available_trials
        ));
    }
    let mut simulator = Simulator::new(num_qubits as usize, num_bits as usize, &mut xof);

    println!(
        "profiling {} Fiat--Shamir batches × {} shots over {} operations and {} qubits...",
        args.trials,
        SHOTS_PER_TRIAL,
        ops.len(),
        num_qubits
    );
    for trial in 0..args.trials {
        simulator.clear_for_shot();
        let batch_start = trial * SHOTS_PER_TRIAL;
        for shot in 0..SHOTS_PER_TRIAL {
            let input = batch_start + shot;
            simulator.set_register(&registers[0], targets[input].0, shot);
            simulator.set_register(&registers[1], targets[input].1, shot);
            simulator.set_register(&registers[2], offsets[input].0, shot);
            simulator.set_register(&registers[3], offsets[input].1, shot);
        }
        apply_profiled(
            &mut simulator,
            &ops,
            &trace.final_op_phases,
            &counts,
            &bins_per_phase,
            &bin_offsets,
            &mut executed_toffoli,
            &mut executed_ccx_total,
            &mut executed_ccz_total,
        );
        verify_trial(
            &mut simulator,
            &registers,
            &expected[batch_start..batch_start + SHOTS_PER_TRIAL],
        )?;
        println!("  trial {}/{} verified", trial + 1, args.trials);
    }

    let binned_toffoli_total: u64 = executed_toffoli.iter().sum();
    let executed_toffoli_total = executed_ccx_total + executed_ccz_total;
    if binned_toffoli_total != executed_toffoli_total {
        return Err(format!(
            "phase bins contain {binned_toffoli_total} executed Toffolis, \
             but the simulator-equivalent full-run counter contains {executed_toffoli_total}"
        ));
    }
    let evaluator_toffoli_total =
        evaluator_toffoli_total(&ops, &registers, num_qubits, num_bits, args.trials);
    if evaluator_toffoli_total != executed_toffoli_total {
        return Err(format!(
            "profile counter contains {executed_toffoli_total} executed Toffolis, \
             but Simulator::apply_iter contains {evaluator_toffoli_total}"
        ));
    }
    let profiled_runs = args.trials * SHOTS_PER_TRIAL;
    println!("  allocator peak cross-check: {peak_active} active = {num_qubits} analyzed qubits");
    println!(
        "  Toffoli cross-check: {binned_toffoli_total} binned = \
         {evaluator_toffoli_total} Simulator::apply_iter"
    );
    println!(
        "  executed-gate split: CCX={executed_ccx_total}, CCZ={executed_ccz_total}, \
         combined={executed_toffoli_total}"
    );
    for phase in 1..=RESOURCE_PHASE_COUNT {
        let start = bin_offsets[phase];
        let end = bin_offsets[phase + 1];
        let phase_total: u64 = executed_toffoli[start..end].iter().sum();
        println!(
            "  phase {phase} {:<26} mean executed Toffoli/run = {:.6}",
            RESOURCE_PHASE_NAMES[phase],
            phase_total as f64 / profiled_runs as f64
        );
    }
    println!(
        "  all-phase mean executed Toffoli/run = {:.6} ({executed_toffoli_total} / {profiled_runs})",
        executed_toffoli_total as f64 / profiled_runs as f64
    );

    write_csv(
        &args.output,
        args.trials,
        &bins_per_phase,
        &bin_offsets,
        &index_sum,
        &op_count,
        &active,
        &executed_toffoli,
    )?;
    println!("wrote {} rows to {}", total_bins, args.output.display());
    Ok(())
}

fn main() {
    let result = parse_args().and_then(run);
    if let Err(error) = result {
        eprintln!("profile_circuit: {error}");
        std::process::exit(2);
    }
}
