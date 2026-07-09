# Toffoli-TODD Phase-Island Findings

Date: 2026-07-01

Goal: test the Vandaele/TODD-style idea on Toffoli/CCZ count by merging shared
2-face cubic phase terms.

## Corrected Model

Literal `CCZ(a,b,c)` terms can merge as `CCZ(a,b,c1) + CCZ(a,b,c2) =
CCZ(a,b,c1 xor c2)` only inside one phase island where:

- the two shared face qubits are the same logical values,
- the third-slot qubits are simultaneously available or unchanged across the
  parity compute/uncompute interval,
- the classical condition stack is the same, or the conditions are explicitly
  folded into a classically gated parity ancilla, and
- no write/reset/HMR/swap touches the face or already-collected third terms
  between collection and uncompute.

Ignoring `PushCondition`/`PopCondition` makes a large but invalid upper bound:
the biggest physical bucket `(768,1151)` had 313 apparent members because
q1151 is a recycled carry lane. With condition-stack context included, each
member is under a different HMR bit, so literal shared-face savings are zero.

## Analyzer Added

`src/bin/todd_face_report.rs` now reports:

- strict condition-aware literal CCZ shared-face savings,
- conditioned-parity write-safe segment candidates,
- diagonal phase-run CCZ counts,
- face-factor segments where a shared product `t = a&b` could replace several
  `CCZ(a,b,c_i)` by `CZ(t,c_i)`,
- targeted face inspection via `TODD_FACE_A/B`,
- source-site counts via `TODD_SITE_COUNTS=1`.

Representative commands:

```bash
SKIP_ALT_SEED_CHECKS=1 target/release/todd_face_report
SKIP_ALT_SEED_CHECKS=1 TODD_FACE_A=768 TODD_FACE_B=1151 TODD_FACE_LIMIT=60 target/release/todd_face_report
SKIP_ALT_SEED_CHECKS=1 TRACE_OP_SITES=1 TODD_SITE_COUNTS=1 target/release/todd_face_report
SKIP_ALT_SEED_CHECKS=1 TRACE_OP_SITES=1 TODD_FACE_FACTOR_SITE='arith.rs:1351' target/release/todd_face_report
SKIP_ALT_SEED_CHECKS=1 TRACE_TLM_ERASE_CARRY=1 target/release/todd_face_report
```

## Results On The d44cad3/71f5115 Local Base

Baseline remains clean:

```text
qubits=1152
avg executed Toffoli=1364229.770
errors=0 classical / 0 phase / 0 ancilla
```

Condition-aware CCZ face report:

```text
literal_ccz_shared_face_exact_diagonal_run_save=0
literal_ccz_shared_face_strict_condition_global_save=0
conditioned_parity_write_safe_save=0
diagonal_ccz_runs multi_run_count=0
face_factor_segments filter=arith.rs:1351 raw_save=0 stack_neutral_save=0
```

Top physical bucket `(768,1151)`:

```text
members=313
unique_thirds=34
cond_contexts=313
source=src/point_add/trailmix_ludicrous/arith.rs:1351
barrier_summary:
  face_write_or_noncommuting=19214
  support_write=41259
  condition_transitions=123488
  conservative_segment_save_sum=0
```

CCZ source-site count:

```text
src/point_add/trailmix_ludicrous/arith.rs:1351 -> 1632 CCZ
```

Measured-carry source trace for that site:

```text
TLM_ERASE_CARRY rows=1632
active range=817..1151
phase counts:
  tlm_apply_forward_mod_add_fold -> 745 rows, active 1074..1150
  tlm_apply_inverse_mod_sub_fold -> 705 rows, active 1074..1151
  square_b_hi_apply_f_times_sub  -> 114 rows, active 945..947
dominant zero-constant chunks:
  width=2 coff=50 ones=0 -> 413 rows
  width=3 coff=47 ones=0 -> 413 rows
```

The tempting source-level cache would specialize all-zero chunks. For width `s`
and zero constant bits, the HMR phase is:

```text
h * ctrl * cin * product_i(!a_i)
```

If several erasures shared the same `(ctrl, a[coff..coff+s])`, we could compute
`p = ctrl * product_i(!a_i)` once with `s` Toffolis, replace each carry phase by
`CZ(p, cin)` under its HMR bit, and HMR-clear the product chain with only
Clifford feedback. That would save roughly `s*(k-1)` Toffoli-class gates for a
stable group of `k` erasures.

The current circuit does not expose such a group cleanly:

- emitted face-factor stable segments are zero even at `arith.rs:1351`,
- the large repeated `(phase, ctrl, top)` groups are across separate fold
  applications where the accumulator changes, not multiple erasures of the same
  stable zero-test inside one suffix add, and
- the hot apply-fold rows are peak-adjacent (`active` up to 1151), so batching
  predicates live would probably spend peak qubits unless a lane can be borrowed.

## Prototype: Zero-Chunk Phase Cache

Added a default-off source prototype:

```text
TLM_ZERO_CHUNK_PHASE_CACHE=start:end[:width:coff:top]
```

For selected all-zero measured-carry erasers, it computes
`p = ctrl * product_i(!a_i)` once, replaces the original comparator phase with
`CZ(p, cin)` under the carry HMR bit, and clears the product chain by HMR +
negative-control `CZ` feedback.

Results:

```text
TLM_ZERO_CHUNK_PHASE_CACHE=1588:1622:2:50:307
  qubits=1152
  errors=19 classical / 140 phase / 0 ancilla

TLM_ZERO_CHUNK_PHASE_CACHE=1622:1622:2:50:307
  qubits=1152
  errors=19 classical / 15 phase / 0 ancilla
```

The single-call failure means this is not merely cross-call logical instability;
the cached all-zero formula or its HMR-clear realization is not equivalent to
the existing comparator under this simulator/circuit model. Leave the prototype
default-off. It is a diagnostic scaffold, not a hunt candidate.

Support-restricted profiling with `SUPPORT_REWRITE_PROFILE_TOP=20000` found no
CCZ one-missing-pattern candidates and no rare-effect CCZ candidates in the
reported set.

## Interpretation

The simple post-pass Toffoli-TODD merge has no exact candidates in the final
literal-CCZ stream. The apparent global buckets are mostly physical allocator
reuse across HMR/reset boundaries, not stable logical phase faces.

The only live route left in this family is a source-level multi-output measured
carry phase compiler that finds a stable shared predicate before the allocator
recycles it. The all-zero chunk formula above is the best concrete shape found,
but no clean group exists in the current emitted/source trace. A per-call rewrite
that computes the final carry and replaces the top `CCZ` by a conditioned `CZ`
is exact but Toffoli-count neutral: it trades one `CCZ` for one top `CCX`.

Do not GPU-hunt this route until a source-level candidate validates clean at the
baked nonce. Current exact saved Toffoli from this route: 0.

No-env rehearsal with the prototype disabled remains clean:

```text
qubits=1152
avg executed Toffoli=1364229.770
errors=0 classical / 0 phase / 0 ancilla
```

## Follow-up: Clean-AND MBU Rewrite

Added `src/bin/and_island_report.rs` to attribute CCX ownership and search for
the strict Gidney temporary-AND shape:

```text
CCX(a,b,t); uses of t as a control; CCX(a,b,t)
```

with `t` a zero ancilla, no writes to `t` before the clear, and no writes to
`a,b` between compute and clear.

Report on the no-env stream:

```text
ccx=1441344 ccz=5341
conservative_mbu_pair_candidates=356
source=src/point_add/trailmix_ludicrous/mcx.rs:80
```

The first candidate window is a Khattar-Gidney prefix ladder fragment:

```text
37299: CCX q563 q1079 -> q1080
37300: CX q1080 -> q564
37301: CCX q563 q1079 -> q1080
37302: CX q1079 -> q563
37303: R q1080
```

Three rewrite variants were tested:

```text
immediate HMR+CZ at clear, keep later R:
  rewritten=356, emitted_tof 1446685 -> 1446329
  errors=15 classical / 14 phase / 0 ancilla

defer HMR+CZ to later reset, require controls stable through reset:
  rewritten=0
  no-env baseline clean

immediate HMR+CZ at clear, delete later R to preserve RNG stream:
  rewritten=254, emitted_tof 1446685 -> 1446431
  errors=18 classical / 16 phase / 0 ancilla
```

A ten-nonce sweep for the 254-save variant did not find a near-clean seed:

```text
2430844 -> 16 classical / 19 phase / 0 ancilla
2430845 -> 17 classical / 14 phase / 0 ancilla
2430846 -> 21 classical / 13 phase / 0 ancilla
2430847 -> 14 classical /  9 phase / 0 ancilla
2430848 -> 14 classical / 13 phase / 0 ancilla
2430849 -> 21 classical / 10 phase / 0 ancilla
2430850 -> 18 classical / 14 phase / 0 ancilla
2430851 -> 17 classical / 10 phase / 0 ancilla
2430852 -> 11 classical /  8 phase / 0 ancilla
2430853 -> 21 classical / 13 phase / 0 ancilla
```

Conclusion: the local MBU identity is a real Toffoli-count reduction, but the
current circuit's support-truncated arithmetic is sensitive to the op-stream
reseeding caused by this exact rewrite. The pass is therefore opt-in only:

```text
MBU_CLEAN_AND_ENABLE=1
MBU_CLEAN_AND_LIMIT=...
MBU_CLEAN_AND_TRACE=1
```

Default remains clean and unchanged. This is a viable GPU/nonce-hunt target only
if combined with a larger reroll search or with a source-level change that keeps
the existing hard-input support margin.

## 2026-07-01 Follow-up: Exact Confirm Batch And KG Source Variant

Added `ISLAND_EXACT_CONFIRM=/path/to/nonces.txt` to `island_search` under
`ISLAND_EXACT_SIM=1`. The file may contain bare nonces or lines like
`EXACT_SURVIVOR <nonce>`; each listed nonce is full exact-simulated and reported
as:

```text
EXACT nonce=<n> shots=9024 cls=<c> pha=<p> anc=<a> avg_tof=<t>
```

This avoids re-sweeping nonce ranges after a prefix screen.

Full-confirmed the `MBU_CLEAN_AND_ENABLE=1` 254-save post-pass route:

```text
1024-prefix range: [2435972, 2442372) partial stop
prefix survivors full-confirmed: 526
exact survivors: 0
best: nonce=2436478 cls=7 pha=4 anc=0 avg_tof=1363993.535
min_cls=6 min_pha=3 min_anc=0
mean cls/pha/anc = 16.211 / 11.498 / 0.000
mean avg-Toffoli = 1363993.180
```

Low-severity rows existed (`7/4/0`, `6/7/0`) but neither `cls` nor `pha`
reached zero anywhere in 526 full confirmations. This is not ready for a GPU
hunt with the current weak 1024-prefix screen.

Limit sweep on an 8-nonce exact panel:

```text
limit 0:   best 14/12/0, avg-Toffoli ~1364244..1364268
limit 16:  best 12/10/0
limit 32:  best 12/13/0
limit 64:  best 12/9/0
limit 96:  best 13/13/0
limit 128: best 13/7/0
limit 192: best 10/12/0 plus one 16/3/0
limit 254: best 7/4/0, rows clustered around 6..9 classical failures
```

The full 254 rewrite remains the best tail in this post-pass family; smaller
limits do not expose a cleaner candidate.

Implemented a default-off source-level KG prefix reverse variant:

```text
TLM_KG_PREFIX_MBU=1
```

It replaces coherent reverse clears of plain KG prefix AND temporaries with
`HMR + cz_if_bit`, while leaving complemented `X;CCX` temporaries coherent. An
aggressive attempt to MBU-clear complemented temporaries produced a structural
failure (`9024 classical / 141 phase / 0 ancilla` on every sampled nonce), so
that case is intentionally not used.

Conservative source-level result:

```text
emitted CCX: 1441344 -> 1440684  (-660)
ops: 10221377 -> 10222853
remaining conservative MBU pair candidates: 0
```

Short exact-prefix screen:

```text
TLM_KG_PREFIX_MBU=1
range [2443000, 2444000), 1024-shot exact prefix
prefix survivors: 42 / 1000
full exact survivors: 0 / 42
best: nonce=2443608 cls=7 pha=9 anc=0 avg_tof=1363601.258
min_cls=7 min_pha=6 min_anc=0
mean cls/pha/anc = 19.571 / 14.071 / 0.000
mean avg-Toffoli = 1363586.019
```

The source-level route has a stronger Toffoli cut than the 254 post-pass, but
its prefix-survivor residual distribution is worse. Do not GPU-hunt it without a
stronger exact/prefilter screen or a repair that lowers the classical residual
tail.

Default-off sanity anchor remains clean:

```text
DIALOG_TAIL_NONCE=50400005525597
errors=0 classical / 0 phase / 0 ancilla
avg executed Toffoli=1364229.770
```

Current best actionable route in this family is still the 254-save post-pass,
but it needs either a much stronger exact-prefix selector (4096/8192 shots or a
proper classical model for the re-rolled hard-input support) or a semantic
source compiler that targets only sites whose reroll distribution keeps
`cls`/`pha` near zero.

## Toffoli Fire-Rate Profiler

Added a default-off exact-sim profiler in `island_search`:

```text
TOF_FIRE_PROFILE=1
TOF_FIRE_PROFILE_TOP=<n>
TOF_FIRE_PROFILE_SITE_FILTER=<substring>
```

It reuses `TRACE_OP_SITES` and, for each `CCX`/`CCZ`, records scorer-charged
shots (`classical condition mask`) versus actual quantum support:

```text
CCX fire = cond & q_control1 & q_control2
CCZ fire = cond & q_control1 & q_control2 & q_target
```

Frontier sanity with profiling:

```text
DIALOG_TAIL_NONCE=50400005525597
EXACT cls=0 pha=0 anc=0 avg_tof=1364229.770
tof_ops=1446685 total_charged=12310809446 total_fire=2830541768
wasted=9480267678 zero_fire_ops=43791
```

High-waste base sites on the frontier hash:

```text
gidney.rs:859  wasted=1.861B charged=2.412B zero_ops=4781
gidney.rs:925  wasted=1.591B charged=2.427B zero_ops=4034
gcd.rs:900     wasted=0.518B charged=0.633B zero_ops=2201
gcd.rs:1232    wasted=0.513B charged=0.629B zero_ops=1685
fused.rs:2205  wasted=0.453B charged=0.591B zero_ops=82
fused.rs:2136  wasted=0.453B charged=0.591B zero_ops=112
```

Important result: frontier-hash zero support is not by itself a safe deletion
certificate because deleting gates re-hashes and re-rolls Fiat-Shamir inputs.
Source-level candidates must be validated under the new hash.

Tested reverse-GCD cswap diagonal extension:

```text
TLM_GCD_REVERSE_DIAGONAL_MIN=264
static: ccx 1441344 -> 1440926, ops 10221377 -> 10220123
frontier exact: cls=22 pha=19 anc=0 avg_tof=1363820.175
256-nonce exact-prefix scan: 0 survivors, best clean_prefix=2264/9024
```

Chunking the added diagonal (`step+bit=264`, steps 9..217) into 16-step bands
also failed 2048-shot exact-prefix checks; even step 217 alone failed
(`clean_prefix=5/2048`). This diagonal is not a clean candidate.

Tested dormant skip knobs:

```text
TLM_GCD_SKIP_EXACT_REVERSE_CSWAPS=1
clean but no count/avg effect on current route

TLM_GIDNEY_SKIP_TOP2_THREAD=1
static ccx 1441344 -> 1422437, but dirty at shot 0

TLM_GIDNEY_SKIP_FULLVENT_TOP2=1
static ccx 1441344 -> 1422963, but dirty at shot 0
```

Next useful direction from this profiler is not another broad diagonal rule.
Use it to drive narrow source hooks with semantic proof or exact-key
minimization, then validate each under reseeded 9024-shot exact sim before any
GPU hunt.
