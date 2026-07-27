//! Reproducible Monte Carlo experiment for the jump-2 GCD iteration budget.
//!
//! This models the exact integer macro-step, without the circuit's scheduled
//! register truncation or approximate comparison windows. It therefore
//! measures the stopping-time distribution behind the original 258-step
//! schedule and the current iteration budget; it does not measure the
//! additional failure probability of those approximations.

use alloy_primitives::U256;
use quantum_ecc::point_add::trailmix_ludicrous::schedule::{BAKED_ITERS, ITERS};
use sha3::{Digest, Keccak256};
use std::collections::BTreeMap;
use std::env;
use std::process::ExitCode;

const DEFAULT_SAMPLES: u64 = 1_000_000;
const DEFAULT_SEED: u64 = 0xECDA_0258;
const DEFAULT_BUDGET: usize = BAKED_ITERS;
const SAMPLE_DOMAIN: &[u8] = b"ecdsafail/jump-2-justify-258/v1";
const MODULUS_HEX: &str = "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F";

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Config {
    samples: u64,
    seed: u64,
    budget: usize,
}

fn usage(program: &str) -> String {
    format!(
        "Usage: {program} [--samples N] [--seed N|0xHEX] [--budget N]\n\
         Defaults: --samples {DEFAULT_SAMPLES} --seed 0x{DEFAULT_SEED:X} \
         --budget {DEFAULT_BUDGET}"
    )
}

fn parse_u64(value: &str, flag: &str) -> Result<u64, String> {
    let parsed = if let Some(hex) = value
        .strip_prefix("0x")
        .or_else(|| value.strip_prefix("0X"))
    {
        u64::from_str_radix(hex, 16)
    } else {
        value.parse()
    };
    parsed.map_err(|_| format!("invalid value for {flag}: {value}"))
}

fn parse_args() -> Result<Option<Config>, String> {
    let mut config = Config {
        samples: DEFAULT_SAMPLES,
        seed: DEFAULT_SEED,
        budget: DEFAULT_BUDGET,
    };
    let mut args = env::args();
    let program = args
        .next()
        .unwrap_or_else(|| "jump2_justify_258".to_owned());
    let mut rest = args.peekable();

    while let Some(flag) = rest.next() {
        if flag == "-h" || flag == "--help" {
            println!("{}", usage(&program));
            return Ok(None);
        }
        let value = rest
            .next()
            .ok_or_else(|| format!("missing value for {flag}\n{}", usage(&program)))?;
        match flag.as_str() {
            "--samples" => config.samples = parse_u64(&value, &flag)?,
            "--seed" => config.seed = parse_u64(&value, &flag)?,
            "--budget" => {
                config.budget = usize::try_from(parse_u64(&value, &flag)?)
                    .map_err(|_| format!("value for {flag} does not fit usize: {value}"))?;
            }
            _ => return Err(format!("unknown option: {flag}\n{}", usage(&program))),
        }
    }

    if config.samples == 0 {
        return Err("--samples must be positive".to_owned());
    }
    Ok(Some(config))
}

fn secp256k1_modulus() -> U256 {
    U256::from_str_radix(MODULUS_HEX, 16).expect("valid secp256k1 modulus")
}

/// Generate a deterministic uniform sample in `[1, q)` by rejection sampling.
fn sample_nonzero_field_element(seed: u64, counter: &mut u64, q: U256) -> U256 {
    loop {
        let mut hasher = Keccak256::new();
        hasher.update(SAMPLE_DOMAIN);
        hasher.update(seed.to_le_bytes());
        hasher.update(counter.to_le_bytes());
        *counter = counter
            .checked_add(1)
            .expect("sample counter exhausted its u64 domain");

        let digest = hasher.finalize();
        let candidate = U256::from_be_slice(&digest);
        if candidate != U256::ZERO && candidate < q {
            return candidate;
        }
    }
}

/// Return the first macro-step boundary at which `(u, v) = (1, 0)`.
///
/// This is the classical counterpart of the circuit's exact jump-2 control
/// flow. At every boundary after step zero, `v` is even. The first shift of a
/// later macro-step is therefore unconditional; the second shift is selected
/// when the quotient is still even.
fn jump2_stopping_time(x: U256, q: U256) -> usize {
    assert!(x != U256::ZERO && x < q);

    let mut u = q;
    let mut v = x;
    let mut step = 0usize;

    while u != U256::from(1) || v != U256::ZERO {
        if step == 0 {
            if !v.bit(0) {
                v >>= 1;
            }
        } else {
            assert!(!v.bit(0), "v must be even at a macro-step boundary");
            v >>= 1;
        }

        if !v.bit(0) {
            v >>= 1;
        }

        let subtract = v.bit(0);
        let swap = subtract && (step == 0 || v < u);
        if swap {
            std::mem::swap(&mut u, &mut v);
        }
        if subtract {
            v -= u;
        }

        assert!(!v.bit(0), "each macro-step must leave v even");
        step += 1;
    }

    step
}

fn quantile(
    histogram: &BTreeMap<usize, u64>,
    samples: u64,
    numerator: u64,
    denominator: u64,
) -> usize {
    let rank = samples
        .saturating_mul(numerator)
        .saturating_add(denominator - 1)
        / denominator;
    let mut cumulative = 0u64;
    for (&steps, &count) in histogram {
        cumulative += count;
        if cumulative >= rank {
            return steps;
        }
    }
    unreachable!("histogram contains every sample")
}

fn wilson_interval_95(successes: u64, trials: u64) -> (f64, f64) {
    let z = 1.959_963_984_540_054_f64;
    let n = trials as f64;
    let p = successes as f64 / n;
    let z2 = z * z;
    let denominator = 1.0 + z2 / n;
    let center = (p + z2 / (2.0 * n)) / denominator;
    let half_width = z * (p * (1.0 - p) / n + z2 / (4.0 * n * n)).sqrt() / denominator;
    (center - half_width, center + half_width)
}

fn samples_over_budget(histogram: &BTreeMap<usize, u64>, budget: usize) -> u64 {
    histogram
        .range((budget.saturating_add(1))..)
        .map(|(_, count)| count)
        .sum()
}

fn run(config: Config) {
    let q = secp256k1_modulus();
    let mut counter = 0u64;
    let mut sum = 0u128;
    let mut sum_squares = 0u128;
    let mut minimum = usize::MAX;
    let mut maximum = 0usize;
    let mut exceeded = 0u64;
    let mut histogram = BTreeMap::<usize, u64>::new();

    for _ in 0..config.samples {
        let x = sample_nonzero_field_element(config.seed, &mut counter, q);
        let steps = jump2_stopping_time(x, q);
        sum += steps as u128;
        sum_squares += (steps * steps) as u128;
        minimum = minimum.min(steps);
        maximum = maximum.max(steps);
        exceeded += u64::from(steps > config.budget);
        *histogram.entry(steps).or_default() += 1;
    }

    let n = config.samples as f64;
    let mean = sum as f64 / n;
    let variance = sum_squares as f64 / n - mean * mean;
    let population_stddev = variance.max(0.0).sqrt();
    let budget_z = (config.budget as f64 - mean) / population_stddev;
    let tail_rate = exceeded as f64 / n;
    let (tail_low, tail_high) = wilson_interval_95(exceeded, config.samples);

    println!("jump-2 exact stopping-time Monte Carlo");
    println!("samples: {}", config.samples);
    println!("seed: 0x{:016X}", config.seed);
    println!("generator: Keccak-256 counter mode, rejection sampled into [1,q)");
    println!("budget: {}", config.budget);
    println!("baked_schedule_budget: {BAKED_ITERS}");
    println!("current_circuit_budget: {ITERS}");
    println!("mean: {mean:.6}");
    println!("population_stddev: {population_stddev:.6}");
    println!("budget_z_score: {budget_z:.6}");
    println!("minimum: {minimum}");
    println!("median: {}", quantile(&histogram, config.samples, 1, 2));
    println!("p90: {}", quantile(&histogram, config.samples, 9, 10));
    println!("p99: {}", quantile(&histogram, config.samples, 99, 100));
    println!(
        "p99.9: {}",
        quantile(&histogram, config.samples, 999, 1_000)
    );
    println!(
        "p99.99: {}",
        quantile(&histogram, config.samples, 9_999, 10_000)
    );
    println!("maximum: {maximum}");
    println!("samples_over_budget: {exceeded}");
    println!("tail_rate: {tail_rate:.9}");
    println!("tail_rate_wilson_95: [{tail_low:.9}, {tail_high:.9}]");
    println!("hash_candidates_consumed: {counter}");
    println!("budget_tail_curve:");
    let last_curve_budget = config.budget.saturating_add(7).max(ITERS);
    for budget in config.budget..=last_curve_budget {
        let count = samples_over_budget(&histogram, budget);
        let rate = count as f64 / n;
        println!("  {budget}: {count} ({rate:.9})");
    }
    println!("tail_histogram:");
    for (&steps, &count) in histogram.range(config.budget.saturating_sub(4)..) {
        println!("  {steps}: {count}");
    }
}

fn main() -> ExitCode {
    match parse_args() {
        Ok(Some(config)) => {
            run(config);
            ExitCode::SUCCESS
        }
        Ok(None) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("{error}");
            ExitCode::FAILURE
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn known_stopping_times_match_the_exact_macro_step() {
        let q = secp256k1_modulus();
        assert_eq!(jump2_stopping_time(U256::from(1), q), 252);
        assert_eq!(jump2_stopping_time(U256::from(2), q), 252);
        assert_eq!(jump2_stopping_time(U256::from(3), q), 242);
        assert_eq!(jump2_stopping_time(U256::from(5), q), 191);
        assert_eq!(jump2_stopping_time(q - U256::from(1), q), 353);
    }

    #[test]
    fn deterministic_sampler_is_reproducible_and_in_range() {
        let q = secp256k1_modulus();
        let mut counter_a = 0;
        let mut counter_b = 0;
        for _ in 0..32 {
            let a = sample_nonzero_field_element(DEFAULT_SEED, &mut counter_a, q);
            let b = sample_nonzero_field_element(DEFAULT_SEED, &mut counter_b, q);
            assert_eq!(a, b);
            assert!(a != U256::ZERO && a < q);
        }
        assert_eq!(counter_a, counter_b);
    }
}
