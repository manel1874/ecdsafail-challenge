# Circuit resource profile

`circuit_resource_profile.tex` draws the six point-addition phases from
`circuit_profile.csv`. The first and final coordinate phases use one-third the
original facet width and remain independently scaled. The four central phases
share one operation-index scale: their displayed widths are proportional to
their final operation counts. The CSV is intentionally compact: every row is
one phase-local bin, while `index` is the mean final-circuit operation index of
the bin.

Regenerate the measured data:

```sh
cargo run --release --bin profile_circuit -- \
  --trials 141 \
  --bins-per-phase 96 \
  --dense-bins-per-phase 288 \
  --output figures/circuit_profile.csv
```

The inverse, square, and forward-multiply phases use the dense count (288
buckets each). The two coordinate phases and `coord_add3x` use 96 buckets each.

One trial is one 64-shot batch from the same Fiat--Shamir stream used by the
challenge evaluator, so 141 trials cover all 9,024 evaluator runs. Each batch
is checked for output correctness, zero global phase, and clean ancillas before
it contributes to the CSV.

This is a new runtime profile and does not reuse the original CCX-only phase
diagnostics. Here, "Toffoli" deliberately means the evaluator's combined
executed count of both `CCX` and `CCZ`.

The profiler also enforces two independent consistency checks:

- The maximum allocator liveness must equal the qubit width reported by
  `analyze_ops`.
- The sum of all phase/bin Toffoli counts must exactly equal the total produced
  by the evaluator's `Simulator::apply_iter`.

The active-qubit value is the exact allocator live count at the midpoint of the
phase-local bin. The Toffoli value is:

```text
total CCX/CCZ executions in this bin across all runs / number of runs
```

It is not divided by the number of circuit operations.

The TikZ figure displays the Toffoli series on a zero-safe
`log10(1 + mean)` axis. This preserves the 250 bins whose measured Toffoli
mean is exactly zero while retaining logarithmic spacing for positive values.

Compile the standalone TikZ figure from the repository root:

```sh
mkdir -p output/circuit_resource_profile
tectonic figures/circuit_resource_profile.tex \
  --outdir output/circuit_resource_profile
```

The CSV columns are:

1. `index`
2. `active_qubit_count_per_index`
3. `avg_toffoli_gates_per_run_in_bin`
4. `phase_name`
