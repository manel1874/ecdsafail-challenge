# Jump-2: empirical justification for 258 macro-steps

This experiment measures how many **exact** jump-2 macro-steps are needed to
reach the stable GCD state

```text
(u, v) = (1, 0)
```

from `(u, v) = (q, x)`, where `q` is the secp256k1 field modulus and `x` is
sampled uniformly from `[1, q)`.

The experiment is implemented by
[`src/bin/jump2_justify_258.rs`](../src/bin/jump2_justify_258.rs). Run it with:

```sh
cargo run --release --bin jump2_justify_258 -- \
  --samples 10000000 \
  --seed 0xECDA0258 \
  --budget 258
```

The defaults use one million samples, the seed above, and budget 258. The
sampler hashes a domain-separated `(seed, counter)` with Keccak-256 and uses
rejection sampling, making the input stream deterministic and uniform over the
nonzero field elements.

On the pulled `main` used for this experiment, `BAKED_ITERS = 258` identifies
the iteration count against which the stored schedules and dead-gate
certificates were originally fitted, while the live circuit has subsequently
moved to `ITERS = 261`. The experiment reports both constants and evaluates the
tail at every budget from 258 through 265.

## Macro-step being measured

The executable directly models the control flow used by the jump-2 GCD:

```text
if this is step 0:
    if is_even(v):
        v /= 2
else:
    # v is even at every later macro-step boundary
    v /= 2

if is_even(v):
    v /= 2

sub = is_odd(v)
swp = sub if this is step 0 else sub and (v < u)

if swp:
    swap(u, v)
if sub:
    v -= u
```

The stopping time is the first macro-step boundary at which `(u, v) = (1, 0)`.
Inputs needing at most 258 steps fit the fixed circuit budget; inputs needing
more do not.

## Ten-million-sample result

The command above produced:

| Statistic | Result |
| --- | ---: |
| Samples | 10,000,000 |
| Mean stopping time | 241.291096 |
| Population standard deviation | 4.987193 |
| Position of budget 258 | 3.350363 standard deviations above the mean |
| Median | 241 |
| 90th percentile | 248 |
| 99th percentile | 253 |
| 99.9th percentile | 257 |
| 99.99th percentile | 260 |
| Minimum / maximum observed | 216 / 268 |
| Samples requiring more than 258 steps | 2,839 |
| Empirical probability of exceeding 258 | 0.000283900 |
| 95% Wilson interval for that probability | [0.000273649, 0.000294535] |

The same sample gives the following budget curve:

| Budget | Samples requiring more steps | Empirical probability |
| ---: | ---: | ---: |
| 258 | 2,839 | 0.000283900 |
| 259 | 1,371 | 0.000137100 |
| 260 | 618 | 0.000061800 |
| 261 | 286 | 0.000028600 |
| 262 | 116 | 0.000011600 |
| 263 | 41 | 0.000004100 |
| 264 | 18 | 0.000001800 |
| 265 | 6 | 0.000000600 |

The tail around the budget was:

```text
steps  samples
254    31352
255    18436
256    10359
257     5612
258     3005
259     1468
260      753
261      332
262      170
263       75
264       23
265       12
266        5
267        0
268        1
```

This confirms the previously quoted approximations:

```text
mean ≈ 241.3
standard deviation ≈ 5.0
258 is ≈ 3.35 standard deviations above the mean
```

It also quantifies the trade-off more clearly: under the uniform-input model,
258 is a probabilistic resource choice with an observed exact-GCD
non-convergence rate of approximately `2.84e-4`, or about one input in 3,522.
The current 261-step setting reduces this exact stopping-time tail by about a
factor of 9.9, to `2.86e-5`, before accounting for other circuit
approximations. Neither value is a worst-case bound. For example, the
deterministic test case `x = q - 1` takes 353 macro-steps.

## Scope and limitations

This experiment intentionally uses full-width integers and exact comparisons.
It isolates the stopping-time question that justifies the loop count. It does
**not** include:

- truncation according to `SCHED_J2`;
- approximate comparisons according to `GAP_J2`;
- later environment-controlled schedule adjustments;
- the distribution of operands produced by a particular point-addition
  workload; or
- quantum-circuit simulation.

Those mechanisms have additional, workload-dependent failure probabilities and
must be evaluated separately. In particular, the result above should not be
described as a proof that the complete circuit succeeds for every field input.

## Checks

Run the focused unit tests with:

```sh
cargo test --bin jump2_justify_258
```

The tests verify deterministic sampling, the sample range, macro-step boundary
invariants, and several independently calculated stopping times, including the
353-step `x = q - 1` case.
