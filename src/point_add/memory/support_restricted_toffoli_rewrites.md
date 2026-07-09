# Support-Restricted Toffoli Rewrite Plan

This benchmark does not score approximate correctness. Any rewrite must still
pass all verifier shots with zero classical mismatches, zero phase garbage, and
zero ancilla garbage. The useful interpretation of "approximation" is therefore
support restriction: replace a nonlinear gate by a cheaper Clifford expression
only when the input pattern where the replacement differs is unreachable on the
tested support.

## One-Missing-Pattern CCX

For `CCX(a,b,t)`, the exact delta is `a & b`. Clifford-only deltas are affine.
Each affine replacement is wrong on exactly one two-bit pattern:

| Missing pattern `(a,b)` | Replacement |
| --- | --- |
| `00` | `X(t); CX(a,t); CX(b,t)` |
| `01` | `CX(b,t)` |
| `10` | `CX(a,t)` |
| `11` | drop |

The replacement is exact if the listed pattern is absent whenever the gate's
effective classical condition is true.

## One-Missing-Pattern CCZ

For `CCZ(a,b,c)`, the exact phase bit is `a & b & c`. If `000` is absent, use:

```text
NEG
Z a
Z b
Z c
CZ a b
CZ a c
CZ b c
```

This quadratic phase is wrong only at `000`. For any other missing pattern,
temporarily `X` the controls whose missing-pattern bit is `1`, apply the `000`
rewrite, then undo those `X` gates.

## Classical-Control Replacement

If a nonlinear quantum control is known to be an actual classical bit `m`, the
rewrite is exact:

```text
CCX(a, m, t) -> PUSH_CONDITION m; CX(a,t); POP_CONDITION
CCZ(a, b, m) -> PUSH_CONDITION m; CZ(a,b); POP_CONDITION
```

When the original op also has a `c_condition`, keep that condition on the inner
`CX`/`CZ`. The pushed bit composes with the surrounding condition stack.

## Implementation

`support_rewrites.rs` is an opt-in post-build pass applied after the existing
drop-dead lists. Final op indices from the profiler therefore match the stream
sent to `eval_circuit`.

Useful knobs:

```text
SUPPORT_REWRITE_PROFILE=1
SUPPORT_REWRITE_PROFILE_SHOTS=9024
SUPPORT_REWRITE_PROFILE_TOP=40

SUPPORT_REWRITE_CCX_MISSING="idx:pattern,..."
SUPPORT_REWRITE_CCZ_MISSING="idx:pattern,..."
SUPPORT_REWRITE_CCX_BITCTRL="idx:slot:bit,..."
SUPPORT_REWRITE_CCZ_BITCTRL="idx:slot:bit,..."
```

`*_FILE` variants are also supported, with comments introduced by `#`.

Validation workflow:

1. Build with `SUPPORT_REWRITE_PROFILE=1` and collect candidates.
2. Apply one or a small batch with the relevant `SUPPORT_REWRITE_*` knob.
3. Run the trusted evaluator.
4. Keep only rewrites that pass all 9024 shots. Any mismatch means the support
   claim was not stable after the op-stream reseed.
