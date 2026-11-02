#!/usr/bin/env bash
# Submit F26 remedy benchmarks (DISP-F26-024 / F26-BM-009).
#
# Part 1 closes F26-H02: full-depth historical oracle replay for the four
# frozen E. coli datasets (R1-001, R1-003, R1-001-Val, R1-003-Val) with
# fastbact at replicate 1.
# Part 2 closes F26-H01: uniform Stage G determinism, three fresh fastbact
# replicates (1, 2, 3) for CS-SA01, CS-BS01, CS-PA01, CS-SE01 so every N=3
# group is pooled from one Git HEAD.
#
# One qfree job per dataset and replicate. Concurrent jobs do not share an
# outdir: the harness merges summary.json without a lock.
#   benchmarks/f26_results/raw/${dataset}_fastbact_r${replicate}
# Collate finished runs with scripts/collate_benchmark_results.py.
#
# The sbatch script already requests this partition shape. The command
# line repeats it so a header drift cannot change the submitted job:
#   partition=qfree  cpus-per-task=8  mem=14G  time=03:00:00
#
#   bash scripts/submit_f26_remedy_benchmarks.sh --self-test
#   bash scripts/submit_f26_remedy_benchmarks.sh --dry-run
#   bash scripts/submit_f26_remedy_benchmarks.sh --oracles-only
#   bash scripts/submit_f26_remedy_benchmarks.sh --stage-g-only

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ROOT=$(cd "${SCRIPT_DIR}/.." && pwd)
MANIFEST="${ROOT}/F26_DATASET_MANIFEST.tsv"
SBATCH_SCRIPT="${ROOT}/scripts/slurm_benchmark_harness.sbatch"

PARTITION="qfree"
CPUS="8"
MEMORY="14G"
WALLTIME="03:00:00"
TOOL="fastbact"

ORACLE_DATASETS=(R1-001 R1-003 R1-001-Val R1-003-Val)
STAGE_G_DATASETS=(CS-SA01 CS-BS01 CS-PA01 CS-SE01)
STAGE_G_REPLICATES=(1 2 3)

DRY_RUN=0
ORACLES_ONLY=0
STAGE_G_ONLY=0

LOADED_IDS=()
TARGETS=()

usage() {
  cat <<'EOF'
Usage: submit_f26_remedy_benchmarks.sh [options]

Submit the F26 remedy benchmarks for DISP-F26-024 / F26-BM-009.
Default: submit both parts, 16 jobs.

Part 1 (F26-H02 oracle replay): datasets R1-001, R1-003, R1-001-Val,
R1-003-Val with fastbact at replicate 1.
Part 2 (F26-H01 uniform Stage G determinism): datasets CS-SA01, CS-BS01,
CS-PA01, CS-SE01 with fastbact at replicates 1, 2, and 3.

  --manifest PATH    Dataset manifest. Default: F26_DATASET_MANIFEST.tsv
  --oracles-only     Submit only Part 1 (F26-H02 oracle replay).
  --stage-g-only     Submit only Part 2 (F26-H01 Stage G determinism).
  --dry-run          Print sbatch commands and the status table; do not submit.
  --self-test        Run built-in option, command, and manifest target tests.
  -h, --help         Show this help.

Every target dataset must be present in the manifest or nothing is
submitted.

Each job runs scripts/slurm_benchmark_harness.sbatch on partition qfree
with 8 CPUs, 14G memory, and a 03:00:00 wall clock. The harness arguments
are --manifest, --dataset, --tool fastbact, --replicate, and
--outdir benchmarks/f26_results/raw/${dataset}_fastbact_r${replicate}.

The status table columns are:
  status  dataset  tool  replicate  job_id  outdir
EOF
}

die() {
  echo "submit_f26_remedy_benchmarks.sh: $*" >&2
  exit 2
}

is_listed() {
  local needle="$1"
  shift
  local item
  for item in "$@"; do
    if [[ "${item}" == "${needle}" ]]; then
      return 0
    fi
  done
  return 1
}

shell_join() {
  local part
  local out=""
  for part in "$@"; do
    out+=$(printf '%q' "${part}")
    out+=" "
  done
  printf '%s' "${out% }"
}

parse_job_id() {
  local text="$1"
  local line
  line="${text%%$'\n'*}"
  line="${line%$'\r'}"
  if [[ "${line}" =~ ^([0-9]+)(;.*)?$ ]]; then
    printf '%s' "${BASH_REMATCH[1]}"
    return 0
  fi
  if [[ "${line}" =~ Submitted\ batch\ job\ ([0-9]+) ]]; then
    printf '%s' "${BASH_REMATCH[1]}"
    return 0
  fi
  return 1
}

require_value() {
  local flag="$1"
  local value="$2"
  [[ -n "${value}" ]] || die "${flag} requires a value"
  [[ "${value}" != --* ]] || die "${flag} requires a value"
}

resolve_manifest_path() {
  local path="$1"
  if [[ "${path}" != /* ]]; then
    if [[ -f "${path}" ]]; then
      path="$(cd "$(dirname "${path}")" && pwd)/$(basename "${path}")"
    else
      path="${ROOT}/${path}"
    fi
  fi
  printf '%s' "${path}"
}

load_manifest() {
  local path="$1"
  local header line id
  local -A seen=()
  LOADED_IDS=()
  [[ -f "${path}" ]] || die "manifest not found: ${path}"
  IFS= read -r header < "${path}" || die "manifest is empty: ${path}"
  header="${header%$'\r'}"
  [[ "${header%%$'\t'*}" == "dataset_id" ]] || die "manifest first column must be dataset_id: ${path}"
  while IFS= read -r line || [[ -n "${line}" ]]; do
    line="${line%$'\r'}"
    [[ -n "${line}" ]] || continue
    [[ "${line}" == "${header}" ]] && continue
    id="${line%%$'\t'*}"
    [[ -n "${id}" ]] || die "manifest row missing dataset_id"
    [[ -z "${seen[${id}]+x}" ]] || die "duplicate dataset_id: ${id}"
    seen["${id}"]=1
    LOADED_IDS+=("${id}")
  done < "${path}"
  [[ ${#LOADED_IDS[@]} -gt 0 ]] || die "manifest has no dataset rows: ${path}"
}

plan_targets() {
  local ds rep
  TARGETS=()
  if [[ "${STAGE_G_ONLY}" -eq 0 ]]; then
    for ds in "${ORACLE_DATASETS[@]}"; do
      TARGETS+=("$(printf 'oracles\t%s\t1' "${ds}")")
    done
  fi
  if [[ "${ORACLES_ONLY}" -eq 0 ]]; then
    for ds in "${STAGE_G_DATASETS[@]}"; do
      for rep in "${STAGE_G_REPLICATES[@]}"; do
        TARGETS+=("$(printf 'stage-g\t%s\t%s' "${ds}" "${rep}")")
      done
    done
  fi
}

check_targets_listed() {
  local target part dataset replicate
  local -A reported=()
  local -a missing=()
  for target in "${TARGETS[@]}"; do
    IFS=$'\t' read -r part dataset replicate <<< "${target}"
    if is_listed "${dataset}" "${LOADED_IDS[@]}"; then
      continue
    fi
    [[ -z "${reported[${dataset}]+x}" ]] || continue
    reported["${dataset}"]=1
    missing+=("${dataset}")
  done
  if [[ ${#missing[@]} -gt 0 ]]; then
    die "manifest is missing target datasets: ${missing[*]}"
  fi
}

row_outdir() {
  local dataset="$1"
  local replicate="$2"
  printf '%s' "${ROOT}/benchmarks/f26_results/raw/${dataset}_${TOOL}_r${replicate}"
}

build_sbatch_cmd() {
  local dataset="$1"
  local replicate="$2"
  local manifest="$3"
  local outdir="$4"
  SBATCH_CMD=(
    sbatch
    --parsable
    --partition="${PARTITION}"
    --cpus-per-task="${CPUS}"
    --mem="${MEMORY}"
    --time="${WALLTIME}"
    --job-name="f26remedy_${dataset}_r${replicate}"
    --chdir="${ROOT}"
    "${SBATCH_SCRIPT}"
    --manifest "${manifest}"
    --dataset "${dataset}"
    --tool "${TOOL}"
    --replicate "${replicate}"
    --outdir "${outdir}"
  )
}

submit_selected() {
  local manifest="$1"
  local target part dataset replicate outdir submission submit_status job_id
  local failed=0
  local -a command_lines=()
  local -a rows=()
  for target in "${TARGETS[@]}"; do
    [[ "${failed}" -eq 0 ]] || break
    IFS=$'\t' read -r part dataset replicate <<< "${target}"
    outdir="$(row_outdir "${dataset}" "${replicate}")"
    build_sbatch_cmd "${dataset}" "${replicate}" "${manifest}" "${outdir}"
    if [[ "${DRY_RUN}" -eq 1 ]]; then
      command_lines+=("$(shell_join "${SBATCH_CMD[@]}")")
      rows+=("$(printf 'DRY_RUN\t%s\t%s\t%s\t-\t%s' "${dataset}" "${TOOL}" "${replicate}" "${outdir}")")
      continue
    fi
    set +e
    submission="$("${SBATCH_CMD[@]}" 2>&1)"
    submit_status=$?
    set -e
    if [[ "${submit_status}" -ne 0 ]]; then
      echo "submit_f26_remedy_benchmarks.sh: sbatch failed for ${dataset} r${replicate}: ${submission}" >&2
      rows+=("$(printf 'FAILED\t%s\t%s\t%s\t-\t%s' "${dataset}" "${TOOL}" "${replicate}" "${outdir}")")
      failed=1
      break
    fi
    if ! job_id="$(parse_job_id "${submission}")"; then
      echo "submit_f26_remedy_benchmarks.sh: unrecognized sbatch output for ${dataset} r${replicate}: ${submission}" >&2
      rows+=("$(printf 'FAILED\t%s\t%s\t%s\t-\t%s' "${dataset}" "${TOOL}" "${replicate}" "${outdir}")")
      failed=1
      break
    fi
    rows+=("$(printf 'SUBMITTED\t%s\t%s\t%s\t%s\t%s' "${dataset}" "${TOOL}" "${replicate}" "${job_id}" "${outdir}")")
  done

  if [[ "${DRY_RUN}" -eq 1 ]]; then
    printf '%s\n' "${command_lines[@]}"
    printf '\n'
  fi
  printf 'status\tdataset\ttool\treplicate\tjob_id\toutdir\n'
  if [[ ${#rows[@]} -gt 0 ]]; then
    printf '%s\n' "${rows[@]}"
  fi
  if [[ "${failed}" -ne 0 ]]; then
    exit 1
  fi
}

parse_args() {
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --manifest)
        [[ $# -ge 2 ]] || die "--manifest requires a value"
        require_value "--manifest" "$2"
        MANIFEST="$2"
        shift 2
        ;;
      --manifest=*)
        MANIFEST="${1#--manifest=}"
        require_value "--manifest" "${MANIFEST}"
        shift
        ;;
      --oracles-only)
        [[ "${STAGE_G_ONLY}" -eq 0 ]] || die "--oracles-only and --stage-g-only are mutually exclusive"
        ORACLES_ONLY=1
        shift
        ;;
      --stage-g-only)
        [[ "${ORACLES_ONLY}" -eq 0 ]] || die "--oracles-only and --stage-g-only are mutually exclusive"
        STAGE_G_ONLY=1
        shift
        ;;
      --dry-run)
        DRY_RUN=1
        shift
        ;;
      -h|--help)
        usage
        exit 0
        ;;
      *)
        die "unknown argument: $1"
        ;;
    esac
  done
}

main() {
  parse_args "$@"
  plan_targets
  MANIFEST="$(resolve_manifest_path "${MANIFEST}")"
  [[ -f "${SBATCH_SCRIPT}" ]] || die "sbatch script not found: ${SBATCH_SCRIPT}"
  load_manifest "${MANIFEST}"
  check_targets_listed
  if [[ "${DRY_RUN}" -eq 0 ]]; then
    command -v sbatch >/dev/null 2>&1 || die "sbatch is not on PATH"
    mkdir -p "${ROOT}/logs"
  fi
  submit_selected "${MANIFEST}"
}

self_test() {
  local checks=0
  local failures=0
  local tmp fixture fixture_missing fixture_dup fixture_bad bindir out status expected
  local rows order tools submitted_rows submitted_ids failed_rows header_rows
  tmp="$(mktemp -d)"
  trap 'rm -rf "${tmp}"' EXIT
  fixture="${tmp}/manifest.tsv"
  fixture_missing="${tmp}/manifest_missing_target.tsv"
  fixture_dup="${tmp}/manifest_duplicate.tsv"
  fixture_bad="${tmp}/manifest_bad_header.tsv"
  bindir="${tmp}/bin"
  mkdir -p "${bindir}"
  ln -s "$(command -v dirname)" "${bindir}/dirname"

  cat > "${fixture}" <<'EOF'
dataset_id	species
R1-001	Escherichia coli
R1-003	Escherichia coli
R1-001-Val	Escherichia coli
R1-003-Val	Escherichia coli
CS-SA01	Staphylococcus aureus
CS-BS01	Bacillus subtilis
CS-PA01	Pseudomonas aeruginosa
CS-SE01	Salmonella enterica
EOF
  grep -v $'^R1-001-Val\t' "${fixture}" > "${fixture_missing}"
  printf 'dataset_id\tspecies\nAAA\tX\nAAA\tX\n' > "${fixture_dup}"
  printf 'id\tspecies\nR1-001\tX\n' > "${fixture_bad}"

  pass() {
    checks=$((checks + 1))
  }
  fail() {
    echo "SELF_TEST_FAIL: $*" >&2
    failures=$((failures + 1))
  }
  run_expect() {
    local name="$1"
    local expected_status="$2"
    shift 2
    set +e
    out="$("$@" 2>&1)"
    status=$?
    set -e
    if [[ "${status}" -ne "${expected_status}" ]]; then
      fail "${name}: exit ${status}, expected ${expected_status}: ${out}"
      return 0
    fi
    pass
    return 0
  }

  # --- option handling ---
  run_expect "unknown argument" 2 bash "${BASH_SOURCE[0]}" --not-a-flag
  run_expect "stray positional" 2 bash "${BASH_SOURCE[0]}" oracles
  run_expect "manifest needs a value" 2 bash "${BASH_SOURCE[0]}" --manifest
  run_expect "empty manifest value" 2 bash "${BASH_SOURCE[0]}" --manifest= --dry-run
  run_expect "part flags are exclusive" 2 \
    bash "${BASH_SOURCE[0]}" --oracles-only --stage-g-only --dry-run
  run_expect "part flags are exclusive in reverse order" 2 \
    bash "${BASH_SOURCE[0]}" --stage-g-only --oracles-only
  run_expect "missing manifest file" 2 \
    bash "${BASH_SOURCE[0]}" --manifest "${tmp}/missing.tsv" --dry-run
  run_expect "manifest missing a target dataset" 2 \
    bash "${BASH_SOURCE[0]}" --manifest "${fixture_missing}" --dry-run
  if [[ "${out}" == *"missing target datasets"* && "${out}" == *"R1-001-Val"* ]]; then
    pass
  else
    fail "missing target message must name R1-001-Val: ${out}"
  fi
  run_expect "duplicate manifest id" 2 \
    bash "${BASH_SOURCE[0]}" --manifest "${fixture_dup}" --dry-run
  run_expect "bad manifest header" 2 \
    bash "${BASH_SOURCE[0]}" --manifest "${fixture_bad}" --dry-run
  run_expect "help" 0 bash "${BASH_SOURCE[0]}" --help
  if [[ "${out}" != *"--self-test"* || "${out}" != *"--oracles-only"* \
    || "${out}" != *"--stage-g-only"* || "${out}" != *"qfree"* ]]; then
    fail "help text is missing required options"
  else
    pass
  fi

  # --- manifest loading and target planning (in process) ---
  load_manifest "${fixture}"
  if [[ ${#LOADED_IDS[@]} -eq 8 && "${LOADED_IDS[0]}" == "R1-001" \
    && "${LOADED_IDS[2]}" == "R1-001-Val" && "${LOADED_IDS[3]}" == "R1-003-Val" \
    && "${LOADED_IDS[7]}" == "CS-SE01" ]]; then
    pass
  else
    fail "fixture manifest load: count=${#LOADED_IDS[@]} first=${LOADED_IDS[0]:-}"
  fi

  ORACLES_ONLY=0
  STAGE_G_ONLY=0
  plan_targets
  if [[ ${#TARGETS[@]} -eq 16 && "${TARGETS[0]}" == $'oracles\tR1-001\t1' \
    && "${TARGETS[15]}" == $'stage-g\tCS-SE01\t3' ]]; then
    pass
  else
    fail "default plan: count=${#TARGETS[@]} first=${TARGETS[0]:-} last=${TARGETS[15]:-}"
  fi

  ORACLES_ONLY=1
  STAGE_G_ONLY=0
  plan_targets
  if [[ ${#TARGETS[@]} -eq 4 && "${TARGETS[0]}" == $'oracles\tR1-001\t1' \
    && "${TARGETS[3]}" == $'oracles\tR1-003-Val\t1' ]]; then
    pass
  else
    fail "oracles-only plan: count=${#TARGETS[@]} first=${TARGETS[0]:-}"
  fi

  ORACLES_ONLY=0
  STAGE_G_ONLY=1
  plan_targets
  if [[ ${#TARGETS[@]} -eq 12 && "${TARGETS[0]}" == $'stage-g\tCS-SA01\t1' \
    && "${TARGETS[11]}" == $'stage-g\tCS-SE01\t3' ]]; then
    pass
  else
    fail "stage-g-only plan: count=${#TARGETS[@]} first=${TARGETS[0]:-}"
  fi
  STAGE_G_ONLY=0

  LOADED_IDS=("R1-001")
  plan_targets
  if ( check_targets_listed ) 2>/dev/null; then
    fail "incomplete manifest should fail the target check"
  else
    pass
  fi

  load_manifest "${fixture}"
  if ( check_targets_listed ); then
    pass
  else
    fail "complete fixture should pass the target check"
  fi

  load_manifest "${ROOT}/F26_DATASET_MANIFEST.tsv"
  if [[ ${#LOADED_IDS[@]} -eq 8 ]]; then
    pass
  else
    fail "real manifest row count: ${#LOADED_IDS[@]}"
  fi
  if ( check_targets_listed ); then
    pass
  else
    fail "real manifest must contain every target dataset"
  fi

  # --- harness resource header drift guard ---
  local harness
  harness="$(<"${SBATCH_SCRIPT}")"
  if [[ "${harness}" == *"#SBATCH --partition=${PARTITION}"* \
    && "${harness}" == *"#SBATCH --cpus-per-task=${CPUS}"* \
    && "${harness}" == *"#SBATCH --mem=${MEMORY}"* \
    && "${harness}" == *"#SBATCH -t ${WALLTIME}"* ]]; then
    pass
  else
    fail "harness resource header drifted from qfree/8/14G/03:00:00"
  fi

  # --- command construction via dry-run with a sentinel sbatch stub ---
  cat > "${bindir}/sbatch" <<'EOF'
#!/bin/bash
echo "sbatch should not run during --dry-run" >&2
exit 99
EOF
  chmod +x "${bindir}/sbatch"

  run_expect "oracles-only dry-run does not execute sbatch" 0 \
    env PATH="${bindir}:/usr/bin:/bin" \
    bash "${BASH_SOURCE[0]}" --manifest "${fixture}" --oracles-only --dry-run
  rows="$(printf '%s\n' "${out}" | awk -F'\t' '$1=="DRY_RUN"{c++} END{print c+0}')"
  order="$(printf '%s\n' "${out}" | awk -F'\t' '$1=="DRY_RUN"{print $2}' | paste -sd, -)"
  tools="$(printf '%s\n' "${out}" | awk -F'\t' '$1=="DRY_RUN"{print $3}' | sort -u | paste -sd, -)"
  expected="${ROOT}/benchmarks/f26_results/raw/R1-001-Val_fastbact_r1"
  if [[ "${rows}" -eq 4 \
    && "${order}" == "R1-001,R1-003,R1-001-Val,R1-003-Val" \
    && "${tools}" == "fastbact" \
    && "${out}" == *"sbatch --parsable --partition=qfree --cpus-per-task=8 --mem=14G --time=03:00:00 --job-name=f26remedy_R1-001_r1 --chdir=${ROOT}"* \
    && "${out}" == *"--manifest ${fixture} --dataset R1-001-Val --tool fastbact --replicate 1 --outdir ${expected}"* \
    && "${out}" == *"DRY_RUN"*$'\t'"R1-001-Val"$'\t'"fastbact"$'\t'"1"$'\t'"-"$'\t'"${expected}"* \
    && "${out}" != *"CS-SA01"* ]]; then
    pass
  else
    fail "oracles-only dry-run: ${out}"
  fi

  run_expect "stage-g-only dry-run" 0 \
    env PATH="${bindir}:/usr/bin:/bin" \
    bash "${BASH_SOURCE[0]}" --manifest "${fixture}" --stage-g-only --dry-run
  rows="$(printf '%s\n' "${out}" | awk -F'\t' '$1=="DRY_RUN"{c++} END{print c+0}')"
  order="$(printf '%s\n' "${out}" | awk -F'\t' '$1=="DRY_RUN"{print $2 ":" $4}' | paste -sd, -)"
  expected="${ROOT}/benchmarks/f26_results/raw/CS-SE01_fastbact_r3"
  if [[ "${rows}" -eq 12 \
    && "${order}" == "CS-SA01:1,CS-SA01:2,CS-SA01:3,CS-BS01:1,CS-BS01:2,CS-BS01:3,CS-PA01:1,CS-PA01:2,CS-PA01:3,CS-SE01:1,CS-SE01:2,CS-SE01:3" \
    && "${out}" == *"--dataset CS-SE01 --tool fastbact --replicate 3 --outdir ${expected}"* \
    && "${out}" == *"--job-name=f26remedy_CS-BS01_r2"* \
    && "${out}" != *"R1-001"* ]]; then
    pass
  else
    fail "stage-g-only dry-run: ${out}"
  fi

  run_expect "full dry-run covers both parts" 0 \
    env PATH="${bindir}:/usr/bin:/bin" \
    bash "${BASH_SOURCE[0]}" --manifest "${fixture}" --dry-run
  rows="$(printf '%s\n' "${out}" | awk -F'\t' '$1=="DRY_RUN"{c++} END{print c+0}')"
  header_rows="$(printf '%s\n' "${out}" | grep -c $'^status\tdataset\ttool\treplicate\tjob_id\toutdir$' || true)"
  if [[ "${rows}" -eq 16 && "${header_rows}" -eq 1 ]]; then
    pass
  else
    fail "full dry-run rows=${rows} header=${header_rows}"
  fi

  run_expect "real manifest dry-run" 0 \
    env PATH="${bindir}:/usr/bin:/bin" \
    bash "${BASH_SOURCE[0]}" --dry-run
  rows="$(printf '%s\n' "${out}" | awk -F'\t' '$1=="DRY_RUN"{c++} END{print c+0}')"
  if [[ "${rows}" -eq 16 \
    && "${out}" == *"--manifest ${ROOT}/F26_DATASET_MANIFEST.tsv --dataset R1-003-Val --tool fastbact --replicate 1 --outdir ${ROOT}/benchmarks/f26_results/raw/R1-003-Val_fastbact_r1"* \
    && "${out}" == *"--dataset CS-PA01 --tool fastbact --replicate 2 --outdir ${ROOT}/benchmarks/f26_results/raw/CS-PA01_fastbact_r2"* ]]; then
    pass
  else
    fail "real manifest dry-run: ${out}"
  fi

  # --- submission path with stub sbatch binaries ---
  cat > "${bindir}/sbatch" <<'EOF'
#!/bin/bash
echo "424242"
EOF
  chmod +x "${bindir}/sbatch"
  run_expect "oracles-only submission" 0 \
    env PATH="${bindir}:/usr/bin:/bin" \
    bash "${BASH_SOURCE[0]}" --manifest "${fixture}" --oracles-only
  submitted_rows="$(printf '%s\n' "${out}" | awk -F'\t' '$1=="SUBMITTED"{c++} END{print c+0}')"
  submitted_ids="$(printf '%s\n' "${out}" | awk -F'\t' '$1=="SUBMITTED"{print $5}' | sort -u | paste -sd, -)"
  if [[ "${submitted_rows}" -eq 4 && "${submitted_ids}" == "424242" && "${out}" != *"sbatch --parsable"* ]]; then
    pass
  else
    fail "oracles-only submission: ${out}"
  fi

  cat > "${bindir}/sbatch" <<'EOF'
#!/bin/bash
echo "Submitted batch job 77"
EOF
  chmod +x "${bindir}/sbatch"
  run_expect "stage-g-only submission with classic sbatch message" 0 \
    env PATH="${bindir}:/usr/bin:/bin" \
    bash "${BASH_SOURCE[0]}" --manifest "${fixture}" --stage-g-only
  submitted_rows="$(printf '%s\n' "${out}" | awk -F'\t' '$1=="SUBMITTED"{c++} END{print c+0}')"
  submitted_ids="$(printf '%s\n' "${out}" | awk -F'\t' '$1=="SUBMITTED"{print $5}' | sort -u | paste -sd, -)"
  if [[ "${submitted_rows}" -eq 12 && "${submitted_ids}" == "77" ]]; then
    pass
  else
    fail "stage-g-only submission: ${out}"
  fi

  cat > "${bindir}/sbatch" <<'EOF'
#!/bin/bash
echo "sbatch: allocation denied" >&2
exit 1
EOF
  chmod +x "${bindir}/sbatch"
  run_expect "sbatch failure stops the batch" 1 \
    env PATH="${bindir}:/usr/bin:/bin" \
    bash "${BASH_SOURCE[0]}" --manifest "${fixture}"
  failed_rows="$(printf '%s\n' "${out}" | awk -F'\t' 'NF>=6 && $1!="status"{c++} END{print c+0}')"
  if [[ "${failed_rows}" -eq 1 \
    && "${out}" == *"FAILED"*$'\t'"R1-001"$'\t'"fastbact"$'\t'"1"$'\t'"-"$'\t'* ]]; then
    pass
  else
    fail "failure stops the batch: ${out}"
  fi

  rm -f "${bindir}/sbatch"
  run_expect "sbatch missing from PATH" 2 \
    env PATH="${bindir}" \
    "$(command -v bash)" "${BASH_SOURCE[0]}" --manifest "${fixture}" --oracles-only
  if [[ "${out}" == *"sbatch is not on PATH"* ]]; then
    pass
  else
    fail "missing sbatch message: ${out}"
  fi

  if [[ -x "${BASH_SOURCE[0]}" ]]; then
    pass
  else
    fail "script is not executable"
  fi

  if [[ "${failures}" -ne 0 ]]; then
    echo "SELF_TEST_FAIL checks=${checks} failures=${failures}" >&2
    exit 1
  fi
  echo "SELF_TEST_PASS checks=${checks}"
  exit 0
}

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
  for arg in "$@"; do
    if [[ "${arg}" == "--self-test" ]]; then
      self_test
    fi
  done
  main "$@"
fi
