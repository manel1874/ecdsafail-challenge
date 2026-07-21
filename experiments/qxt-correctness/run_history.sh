#!/usr/bin/env bash
set -euo pipefail

repo="$(git rev-parse --show-toplevel)"
output="${1:-${repo}/experiments/qxt-correctness/qxt-history.csv}"
points="${QXT_POINTS:-50000}"
seed="${QXT_SEED:-ecdsafail-qxt-history-v1}"
threads="${QXT_THREADS:-8}"
run_tmp="$(mktemp -d "${TMPDIR:-/tmp}/ecdsafail-qxt.XXXXXX")"
history_worktree="${run_tmp}/history-worktree"
scratch="${run_tmp}/circuit-output"
accepts="${run_tmp}/accepted.tsv"
sampled="${run_tmp}/sampled.tsv"
evaluator="${run_tmp}/qxt_correctness_experiment"
corpus="${run_tmp}/points.bin"

cleanup() {
  if [[ -d "${history_worktree}" ]]; then
    git -C "${repo}" worktree remove --force "${history_worktree}" >/dev/null 2>&1 || true
  fi
  if [[ -n "${run_tmp}" && -d "${run_tmp}" ]]; then
    rm -rf -- "${run_tmp}"
  fi
}
trap cleanup EXIT

mkdir -p "$(dirname "${output}")" "${scratch}"

git -C "${repo}" log --first-parent --reverse \
  --format='%H%x09%h%x09%cI%x09%s' \
  | awk -F '\t' '$4 ~ /^Accept submission /' > "${accepts}"
accept_count="$(wc -l < "${accepts}" | tr -d ' ')"
if [[ "${accept_count}" -eq 0 ]]; then
  echo "no first-parent 'Accept submission' commits found" >&2
  exit 1
fi

# Select 1, 11, 21, ... and always include the last accepted submission.
# With the current 405 accepted submissions this yields 42 graph points.
awk -F '\t' -v total="${accept_count}" \
  'NR == 1 || NR % 10 == 1 || NR == total { print NR "\t" $0 }' \
  "${accepts}" > "${sampled}"

compiler=""
for candidate in "${CC:-}" gcc cc clang; do
  if [[ -n "${candidate}" ]] && command -v "${candidate}" >/dev/null 2>&1; then
    compiler="$(command -v "${candidate}")"
    break
  fi
done
if [[ -z "${compiler}" ]]; then
  echo "no C compiler/linker found" >&2
  exit 1
fi

export CC="${compiler}"
export CARGO_NET_OFFLINE=true
export RUSTFLAGS="-C linker=${compiler} -Awarnings"

echo "Building fixed-corpus evaluator"
(cd "${repo}" && cargo build --release --locked --offline --bin qxt_correctness_experiment)
cp "${repo}/target/release/qxt_correctness_experiment" "${evaluator}"

first_commit="$(sed -n '1s/^[^	]*	\([^	]*\).*/\1/p' "${sampled}")"
git -C "${repo}" worktree add --detach "${history_worktree}" "${first_commit}" >/dev/null

while IFS=$'\t' read -r accepted_index commit short commit_date subject; do
  submission_id="${subject#Accept submission }"
  if [[ -f "${output}" ]] && rg -F -q "${commit}" "${output}"; then
    echo "Skipping accepted submission ${accepted_index}/${accept_count} (${short}); already in CSV"
    continue
  fi

  echo "Evaluating accepted submission ${accepted_index}/${accept_count} (${short})"
  git -C "${history_worktree}" checkout --detach --force --quiet "${commit}"
  (cd "${history_worktree}" && cargo build --quiet --release --locked --offline --bin build_circuit)

  rm -f -- "${scratch}/ops.bin" "${scratch}/ops.bin.tmp"
  (
    cd "${scratch}"
    SKIP_ALT_SEED_CHECKS=1 "${history_worktree}/target/release/build_circuit"
  ) > "${run_tmp}/build-${short}.log" 2>&1

  "${evaluator}" \
    --ops "${scratch}/ops.bin" \
    --corpus "${corpus}" \
    --csv "${output}" \
    --points "${points}" \
    --seed "${seed}" \
    --threads "${threads}" \
    --accepted-index "${accepted_index}" \
    --commit "${commit}" \
    --commit-date "${commit_date}" \
    --submission-id "${submission_id}"
done < "${sampled}"

echo "Wrote ${output}"
