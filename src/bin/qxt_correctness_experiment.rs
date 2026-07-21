//! Evaluate a serialized point-addition circuit on a fixed, circuit-independent
//! corpus and append the result to a history CSV.
//!
//! This is intentionally separate from the challenge evaluator. The official
//! evaluator derives its test inputs from the op stream; this experiment uses
//! the same 50,000 point pairs for every historical circuit so nonce hunting
//! cannot change the measured correctness probability.

use alloy_primitives::U256;
use quantum_ecc::circuit::{
    analyze_ops, BitId, Op, OperationType, QubitId, QubitOrBit, RegisterId, NO_BIT,
};
use quantum_ecc::sim::Simulator;
use quantum_ecc::weierstrass_elliptic_curve::WeierstrassEllipticCurve;
use sha3::{
    digest::{ExtendableOutput, Update, XofReader},
    Shake256,
};
use std::fs::{File, OpenOptions};
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Instant;

const RAW_OPS_MAGIC: &[u8; 8] = b"QECCOPS1";
const ZSTD_OPS_MAGIC: &[u8; 8] = b"QECCOPSZ";
const CORPUS_MAGIC: &[u8; 8] = b"QXTPTS1\0";
const OP_BYTES: usize = 56;
const POINT_VALUES: usize = 6;
const POINT_BYTES: usize = 32;
const MAX_OPS: u64 = 4_000_000_000;
const ZSTD_WINDOW_LOG_MAX: u32 = 27;

#[derive(Debug)]
struct Args {
    ops: PathBuf,
    corpus: PathBuf,
    csv: PathBuf,
    points: usize,
    seed: String,
    threads: usize,
    accepted_index: usize,
    commit: String,
    commit_date: String,
    submission_id: String,
}

fn usage() -> ! {
    eprintln!(
        "usage: qxt_correctness_experiment \\\n+         --ops PATH --corpus PATH --csv PATH --points N --seed TEXT \\\n+         --threads N --accepted-index N --commit HASH \\\n+         --commit-date ISO8601 --submission-id UUID"
    );
    std::process::exit(2);
}

fn parse_args() -> Args {
    let mut args = std::env::args().skip(1);
    let mut get = |wanted: &str| -> String {
        loop {
            let Some(flag) = args.next() else { usage() };
            if flag == wanted {
                return args.next().unwrap_or_else(|| usage());
            }
            eprintln!("unexpected argument {flag:?}; expected {wanted:?}");
            usage();
        }
    };

    let ops = PathBuf::from(get("--ops"));
    let corpus = PathBuf::from(get("--corpus"));
    let csv = PathBuf::from(get("--csv"));
    let points = get("--points").parse().unwrap_or_else(|_| usage());
    let seed = get("--seed");
    let threads = get("--threads").parse().unwrap_or_else(|_| usage());
    let accepted_index = get("--accepted-index").parse().unwrap_or_else(|_| usage());
    let commit = get("--commit");
    let commit_date = get("--commit-date");
    let submission_id = get("--submission-id");
    if args.next().is_some() || points == 0 || threads == 0 {
        usage();
    }
    Args {
        ops,
        corpus,
        csv,
        points,
        seed,
        threads,
        accepted_index,
        commit,
        commit_date,
        submission_id,
    }
}

fn op_kind_from_u32(v: u32) -> Option<OperationType> {
    Some(match v {
        0 => OperationType::Neg,
        1 => OperationType::Register,
        2 => OperationType::AppendToRegister,
        3 => OperationType::BitInvert,
        4 => OperationType::BitStore0,
        5 => OperationType::BitStore1,
        6 => OperationType::X,
        7 => OperationType::Z,
        8 => OperationType::CX,
        9 => OperationType::CZ,
        10 => OperationType::Swap,
        11 => OperationType::R,
        12 => OperationType::Hmr,
        13 => OperationType::CCX,
        14 => OperationType::CCZ,
        15 => OperationType::PushCondition,
        16 => OperationType::PopCondition,
        17 => OperationType::DebugPrint,
        _ => return None,
    })
}

fn read_u64(bytes: &[u8], offset: usize) -> u64 {
    u64::from_le_bytes(bytes[offset..offset + 8].try_into().unwrap())
}

fn read_ops_body(mut reader: impl Read, count: usize) -> Result<Vec<Op>, String> {
    let mut ops = Vec::with_capacity(count);
    let mut record = [0u8; OP_BYTES];
    for index in 0..count {
        reader
            .read_exact(&mut record)
            .map_err(|e| format!("op {index}: short read: {e}"))?;
        let kind_raw = u32::from_le_bytes(record[0..4].try_into().unwrap());
        let kind = op_kind_from_u32(kind_raw)
            .ok_or_else(|| format!("op {index}: unknown kind {kind_raw}"))?;
        let op = Op {
            kind,
            q_control2: QubitId(read_u64(&record, 8)),
            q_control1: QubitId(read_u64(&record, 16)),
            q_target: QubitId(read_u64(&record, 24)),
            c_target: BitId(read_u64(&record, 32)),
            c_condition: BitId(read_u64(&record, 40)),
            r_target: RegisterId(read_u64(&record, 48)),
        };
        let validated = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| op.validate()));
        if validated.is_err() {
            return Err(format!("op {index}: validation failed"));
        }
        ops.push(op);
    }
    let mut extra = [0u8; 1];
    match reader.read(&mut extra) {
        Ok(0) => Ok(ops),
        Ok(_) => Err(format!("trailing data after {count} ops")),
        Err(e) => Err(format!("checking trailing data: {e}")),
    }
}

fn load_ops(path: &Path) -> Result<Vec<Op>, String> {
    let mut file = File::open(path).map_err(|e| format!("open {}: {e}", path.display()))?;
    let mut header = [0u8; 16];
    file.read_exact(&mut header)
        .map_err(|e| format!("read {} header: {e}", path.display()))?;
    let magic: &[u8; 8] = header[0..8].try_into().unwrap();
    let count = u64::from_le_bytes(header[8..16].try_into().unwrap());
    if count > MAX_OPS {
        return Err(format!("op count {count} exceeds cap {MAX_OPS}"));
    }
    let count = count as usize;
    if magic == RAW_OPS_MAGIC {
        read_ops_body(BufReader::new(file), count)
    } else if magic == ZSTD_OPS_MAGIC {
        let mut decoder = zstd::stream::read::Decoder::new(BufReader::new(file))
            .map_err(|e| format!("zstd init: {e}"))?;
        decoder
            .window_log_max(ZSTD_WINDOW_LOG_MAX)
            .map_err(|e| format!("zstd window cap: {e}"))?;
        read_ops_body(decoder, count)
    } else {
        Err(format!("{}: unsupported ops magic", path.display()))
    }
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

#[derive(Clone, Copy)]
struct PointCase {
    target_x: U256,
    target_y: U256,
    offset_x: U256,
    offset_y: U256,
    expected_x: U256,
    expected_y: U256,
}

impl PointCase {
    fn values(self) -> [U256; POINT_VALUES] {
        [
            self.target_x,
            self.target_y,
            self.offset_x,
            self.offset_y,
            self.expected_x,
            self.expected_y,
        ]
    }

    fn from_values(values: [U256; POINT_VALUES]) -> Self {
        Self {
            target_x: values[0],
            target_y: values[1],
            offset_x: values[2],
            offset_y: values[3],
            expected_x: values[4],
            expected_y: values[5],
        }
    }
}

fn corpus_xof(seed: &str) -> sha3::Shake256Reader {
    let mut hasher = Shake256::default();
    hasher.update(b"quantum-ecc-qxt-correctness-corpus-v1");
    hasher.update(&(seed.len() as u64).to_le_bytes());
    hasher.update(seed.as_bytes());
    hasher.finalize_xof()
}

fn generate_corpus(count: usize, seed: &str) -> Vec<PointCase> {
    let curve = secp256k1();
    let mut xof = corpus_xof(seed);
    let mut cases = Vec::with_capacity(count);
    while cases.len() < count {
        let mut scalar_bytes = [[0u8; 32]; 2];
        XofReader::read(&mut xof, &mut scalar_bytes[0]);
        XofReader::read(&mut xof, &mut scalar_bytes[1]);
        let target = curve.mul(curve.gx, curve.gy, U256::from_le_bytes(scalar_bytes[0]));
        let offset = curve.mul(curve.gx, curve.gy, U256::from_le_bytes(scalar_bytes[1]));
        if target.0 == offset.0
            || (target.0.is_zero() && target.1.is_zero())
            || (offset.0.is_zero() && offset.1.is_zero())
        {
            continue;
        }
        let expected = curve.add(target.0, target.1, offset.0, offset.1);
        cases.push(PointCase {
            target_x: target.0,
            target_y: target.1,
            offset_x: offset.0,
            offset_y: offset.1,
            expected_x: expected.0,
            expected_y: expected.1,
        });
    }
    cases
}

fn write_corpus(path: &Path, cases: &[PointCase]) -> Result<(), String> {
    let temporary = path.with_extension("tmp");
    let mut writer = BufWriter::new(
        File::create(&temporary).map_err(|e| format!("create {}: {e}", temporary.display()))?,
    );
    writer
        .write_all(CORPUS_MAGIC)
        .and_then(|_| writer.write_all(&(cases.len() as u64).to_le_bytes()))
        .map_err(|e| format!("write {} header: {e}", temporary.display()))?;
    for case in cases {
        for value in case.values() {
            writer
                .write_all(&value.to_le_bytes::<POINT_BYTES>())
                .map_err(|e| format!("write {}: {e}", temporary.display()))?;
        }
    }
    writer
        .flush()
        .map_err(|e| format!("flush {}: {e}", temporary.display()))?;
    drop(writer);
    std::fs::rename(&temporary, path)
        .map_err(|e| format!("rename {} to {}: {e}", temporary.display(), path.display()))
}

fn read_corpus(path: &Path, expected_count: usize) -> Result<Vec<PointCase>, String> {
    let mut reader = BufReader::new(
        File::open(path).map_err(|e| format!("open corpus {}: {e}", path.display()))?,
    );
    let mut header = [0u8; 16];
    reader
        .read_exact(&mut header)
        .map_err(|e| format!("read corpus header: {e}"))?;
    if &header[0..8] != CORPUS_MAGIC {
        return Err("corpus has bad magic".into());
    }
    let count = u64::from_le_bytes(header[8..16].try_into().unwrap()) as usize;
    if count != expected_count {
        return Err(format!(
            "corpus contains {count} cases, expected {expected_count}"
        ));
    }
    let mut cases = Vec::with_capacity(count);
    for index in 0..count {
        let mut values = [U256::ZERO; POINT_VALUES];
        for value in &mut values {
            let mut bytes = [0u8; POINT_BYTES];
            reader
                .read_exact(&mut bytes)
                .map_err(|e| format!("read corpus case {index}: {e}"))?;
            *value = U256::from_le_bytes(bytes);
        }
        cases.push(PointCase::from_values(values));
    }
    let mut extra = [0u8; 1];
    if reader.read(&mut extra).map_err(|e| e.to_string())? != 0 {
        return Err("corpus has trailing data".into());
    }
    Ok(cases)
}

fn load_or_create_corpus(path: &Path, count: usize, seed: &str) -> Result<Vec<PointCase>, String> {
    if path.exists() {
        read_corpus(path, count)
    } else {
        let cases = generate_corpus(count, seed);
        write_corpus(path, &cases)?;
        Ok(cases)
    }
}

fn simulation_xof(seed: &str, batch: usize) -> sha3::Shake256Reader {
    let mut hasher = Shake256::default();
    hasher.update(b"quantum-ecc-qxt-correctness-simulator-v1");
    hasher.update(&(seed.len() as u64).to_le_bytes());
    hasher.update(seed.as_bytes());
    hasher.update(&(batch as u64).to_le_bytes());
    hasher.finalize_xof()
}

/// Equivalent to `Simulator::apply_iter`, except inactive lanes in a partial
/// final batch are masked out. This makes gate averages exact for 50,000 shots,
/// which is not divisible by the simulator's 64-lane packing width.
fn apply_ops_masked<R: XofReader>(sim: &mut Simulator<'_, R>, ops: &[Op], active_mask: u64) {
    let mut condition_stack = Vec::new();
    let mut current_base_condition = active_mask;
    for op in ops {
        let mut condition = current_base_condition;
        if op.c_condition != NO_BIT {
            condition &= sim.bit(op.c_condition);
        }
        let executed_shots = condition.count_ones() as u64;
        match op.kind {
            OperationType::CCZ | OperationType::CCX => {
                sim.stats.toffoli_gates += executed_shots;
            }
            OperationType::CX
            | OperationType::CZ
            | OperationType::Swap
            | OperationType::R
            | OperationType::Hmr => {
                sim.stats.clifford_gates += executed_shots;
            }
            _ => {}
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
                XofReader::read(&mut *sim.xof, &mut bytes);
                let random = u64::from_le_bytes(bytes);
                *sim.bit_mut(op.c_target) &= !condition;
                *sim.bit_mut(op.c_target) ^= random & condition;
                sim.phase ^= sim.qubit(op.q_target) & random & condition;
                *sim.qubit_mut(op.q_target) &= !condition;
            }
            OperationType::R => {
                let mut bytes = [0u8; 8];
                XofReader::read(&mut *sim.xof, &mut bytes);
                let random = u64::from_le_bytes(bytes);
                sim.phase ^= sim.qubit(op.q_target) & random & condition;
                *sim.qubit_mut(op.q_target) &= !condition;
            }
            OperationType::BitInvert => *sim.bit_mut(op.c_target) ^= condition,
            OperationType::BitStore0 => *sim.bit_mut(op.c_target) &= !condition,
            OperationType::BitStore1 => *sim.bit_mut(op.c_target) |= condition,
            OperationType::AppendToRegister
            | OperationType::Register
            | OperationType::DebugPrint => {}
            OperationType::PushCondition => {
                condition_stack.push(current_base_condition);
                current_base_condition &= sim.bit(op.c_condition);
            }
            OperationType::PopCondition => {
                current_base_condition = condition_stack
                    .pop()
                    .expect("condition-stack underflow in validated circuit");
            }
        }
    }
    assert!(condition_stack.is_empty(), "unbalanced condition stack");
}

#[derive(Default)]
struct RunReport {
    correct: usize,
    toffoli: u64,
    clifford: u64,
    phase_garbage_batches: usize,
    ancilla_garbage_batches: usize,
}

fn validate_layout(regs: &[Vec<QubitOrBit>]) -> Result<(), String> {
    if regs.len() != 4 {
        return Err(format!("expected 4 registers, got {}", regs.len()));
    }
    for (index, register) in regs.iter().enumerate() {
        if register.len() != 256 {
            return Err(format!(
                "register {index} has width {}, expected 256",
                register.len()
            ));
        }
    }
    if !regs[0..2]
        .iter()
        .flatten()
        .all(|item| matches!(item, QubitOrBit::Qubit(_)))
    {
        return Err("target registers must contain qubits".into());
    }
    if !regs[2..4]
        .iter()
        .flatten()
        .all(|item| matches!(item, QubitOrBit::Bit(_)))
    {
        return Err("offset registers must contain classical bits".into());
    }
    Ok(())
}

fn run_corpus(
    ops: &[Op],
    cases: &[PointCase],
    seed: &str,
    threads: usize,
) -> Result<(RunReport, u64), String> {
    let (qubits, bits, _register_count, regs) = analyze_ops(ops.iter());
    validate_layout(&regs)?;
    const BATCH: usize = 64;
    let batches = cases.len().div_ceil(BATCH);
    let next_batch = AtomicUsize::new(0);
    let worker_count = threads.min(batches);
    let reports: Vec<RunReport> = std::thread::scope(|scope| {
        let mut handles = Vec::with_capacity(worker_count);
        for _ in 0..worker_count {
            handles.push(scope.spawn(|| {
                let mut report = RunReport::default();
                loop {
                    let batch = next_batch.fetch_add(1, Ordering::Relaxed);
                    if batch >= batches {
                        break;
                    }
                    let start = batch * BATCH;
                    let end = (start + BATCH).min(cases.len());
                    let batch_cases = &cases[start..end];
                    let active_mask = if batch_cases.len() == 64 {
                        u64::MAX
                    } else {
                        (1u64 << batch_cases.len()) - 1
                    };
                    let mut xof = simulation_xof(seed, batch);
                    let mut sim = Simulator::new(qubits as usize, bits as usize, &mut xof);
                    for (shot, case) in batch_cases.iter().enumerate() {
                        sim.set_register(&regs[0], case.target_x, shot);
                        sim.set_register(&regs[1], case.target_y, shot);
                        sim.set_register(&regs[2], case.offset_x, shot);
                        sim.set_register(&regs[3], case.offset_y, shot);
                    }
                    apply_ops_masked(&mut sim, ops, active_mask);
                    for (shot, case) in batch_cases.iter().enumerate() {
                        if sim.get_register(&regs[0], shot) == case.expected_x
                            && sim.get_register(&regs[1], shot) == case.expected_y
                        {
                            report.correct += 1;
                        }
                    }
                    if sim.phase & active_mask != 0 {
                        report.phase_garbage_batches += 1;
                    }
                    for register in &regs {
                        for item in register {
                            if let QubitOrBit::Qubit(qubit) = item {
                                *sim.qubit_mut(*qubit) = 0;
                            }
                        }
                    }
                    if (0..qubits).any(|q| sim.qubit(QubitId(q)) & active_mask != 0) {
                        report.ancilla_garbage_batches += 1;
                    }
                    report.toffoli += sim.stats.toffoli_gates;
                    report.clifford += sim.stats.clifford_gates;
                }
                report
            }));
        }
        handles
            .into_iter()
            .map(|handle| handle.join().expect("experiment worker panicked"))
            .collect()
    });
    let mut combined = RunReport::default();
    for report in reports {
        combined.correct += report.correct;
        combined.toffoli += report.toffoli;
        combined.clifford += report.clifford;
        combined.phase_garbage_batches += report.phase_garbage_batches;
        combined.ancilla_garbage_batches += report.ancilla_garbage_batches;
    }
    Ok((combined, qubits))
}

fn csv_string(value: &str) -> String {
    format!("\"{}\"", value.replace('"', "\"\""))
}

fn append_csv(
    args: &Args,
    report: &RunReport,
    qubits: u64,
    ops: usize,
    elapsed_seconds: f64,
) -> Result<(), String> {
    let needs_header = std::fs::metadata(&args.csv)
        .map(|metadata| metadata.len() == 0)
        .unwrap_or(true);
    let mut output = OpenOptions::new()
        .create(true)
        .append(true)
        .open(&args.csv)
        .map_err(|e| format!("open CSV {}: {e}", args.csv.display()))?;
    if needs_header {
        writeln!(
            output,
            "accepted_submission_index,commit,commit_short,commit_date,submission_id,tested_point_pairs,correct_point_pairs,incorrect_point_pairs,correctness_probability,correctness_percent,avg_executed_toffoli,peak_qubits,qxt,qxt_over_p,emitted_ops,avg_executed_clifford,phase_garbage_batches,ancilla_garbage_batches,corpus_seed,elapsed_seconds"
        )
        .map_err(|e| format!("write CSV header: {e}"))?;
    }
    let tested = args.points;
    let incorrect = tested - report.correct;
    let probability = report.correct as f64 / tested as f64;
    let avg_toffoli = report.toffoli as f64 / tested as f64;
    let avg_clifford = report.clifford as f64 / tested as f64;
    let qxt = avg_toffoli * qubits as f64;
    let qxt_over_p = if probability == 0.0 {
        f64::INFINITY
    } else {
        qxt / probability
    };
    let short = &args.commit[..args.commit.len().min(7)];
    writeln!(
        output,
        "{},{},{},{},{},{},{},{},{:.12},{:.9},{:.9},{},{:.9},{:.9},{},{:.9},{},{},{},{:.3}",
        args.accepted_index,
        csv_string(&args.commit),
        csv_string(short),
        csv_string(&args.commit_date),
        csv_string(&args.submission_id),
        tested,
        report.correct,
        incorrect,
        probability,
        probability * 100.0,
        avg_toffoli,
        qubits,
        qxt,
        qxt_over_p,
        ops,
        avg_clifford,
        report.phase_garbage_batches,
        report.ancilla_garbage_batches,
        csv_string(&args.seed),
        elapsed_seconds,
    )
    .map_err(|e| format!("append CSV row: {e}"))
}

fn main() {
    let args = parse_args();
    let started = Instant::now();
    let ops = load_ops(&args.ops).unwrap_or_else(|e| {
        eprintln!("failed to load ops: {e}");
        std::process::exit(1);
    });
    let cases = load_or_create_corpus(&args.corpus, args.points, &args.seed).unwrap_or_else(|e| {
        eprintln!("failed to load/create corpus: {e}");
        std::process::exit(1);
    });
    let (report, qubits) = run_corpus(&ops, &cases, &args.seed, args.threads).unwrap_or_else(|e| {
        eprintln!("failed to evaluate circuit: {e}");
        std::process::exit(1);
    });
    let elapsed = started.elapsed().as_secs_f64();
    append_csv(&args, &report, qubits, ops.len(), elapsed).unwrap_or_else(|e| {
        eprintln!("failed to write result: {e}");
        std::process::exit(1);
    });

    let probability = report.correct as f64 / args.points as f64;
    let avg_toffoli = report.toffoli as f64 / args.points as f64;
    let qxt = avg_toffoli * qubits as f64;
    println!(
        "RESULT index={} commit={} correct={}/{} p={:.9} Q={} T={:.6} QxT={:.3} QxT/p={:.3} elapsed={:.1}s",
        args.accepted_index,
        &args.commit[..args.commit.len().min(7)],
        report.correct,
        args.points,
        probability,
        qubits,
        avg_toffoli,
        qxt,
        qxt / probability,
        elapsed,
    );
}
