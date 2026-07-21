# Fixed-corpus Q×T correctness experiment

This experiment measures historical accepted circuits against one deterministic,
circuit-independent corpus. It is deliberately separate from the official
Fiat–Shamir validation corpus, so changing or hunting an op-stream nonce cannot
change the points being tested.

The history runner selects accepted submissions 1, 11, 21, and so on, and always
adds the last accepted submission. The repository currently has 405 accepted
submissions on the first-parent history, producing 42 sampled commits.

For each sampled circuit:

- exactly 50,000 valid secp256k1 point-addition pairs are tested;
- `correctness_probability` (`p`) is the fraction with the expected classical
  `(x, y)` output;
- `qxt` is average executed Toffoli count (`T`) times peak qubits (`Q`);
- `qxt_over_p` is `qxt / correctness_probability`;
- phase and ancilla-garbage batch counts are included as diagnostics but are not
  folded into `p`.

Run (the CSV is resumable):

```bash
experiments/qxt-correctness/run_history.sh
```

Optional environment overrides are `QXT_POINTS`, `QXT_SEED`, and `QXT_THREADS`.
The default corpus seed is `ecdsafail-qxt-history-v1`.

## LaTeX plot

`qxt-history-plot.tex` is a standalone PGFPlots figure that reads the CSV
directly. Compile it from this directory so the relative CSV path resolves:

```bash
tectonic --outdir ../../output/pdf qxt-history-plot.tex
```

The main axis contains the requested `Q×T` and `Q×T/p` lines. Because their
separation is small at the full score scale, an inset shows the relative
correctness penalty in percent.

`qxt-history-section.tex` provides a complete discussion section with the plot
embedded. Its compiled preview is written to `output/pdf/qxt-history-section.pdf`:

```bash
tectonic --outdir ../../output/pdf qxt-history-section.tex
```
