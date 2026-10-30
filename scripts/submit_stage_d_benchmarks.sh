#!/usr/bin/env bash
# Submit Stage D baseline benchmarks: one qfree job per dataset and tool.
#
# Concurrent jobs do not share an outdir. The harness merges summary.json
# without a lock, so each job writes:
#   benchmarks/f26_results/raw/${dataset}_${tool}_r1
# Collate those runs afterwards with scripts/collate_benchmark_results.py.
#
#   bash scripts/submit_stage_d_benchmarks.sh --dry-run
#   bash scripts/submit_stage_d_benchmarks.sh --dataset CS-SA01 --tool spades

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ROOT=$(cd "${SCRIPT_DIR}/.." && pwd)
MANIFEST="${ROOT}/F26_DATASET_MANIFEST.tsv"
SBATCH_SCRIPT="${ROOT}/scripts/slurm_benchmark_harness.sbatch"
REPLICATE=1

DATASETS=(CS-SA01 CS-BS01 CS-PA01 CS-SE01)
TOOLS=(fastbact spades skesa unicycler)

DATASET=""
TOOL=""
DRY_RUN=0

usage() {
  cat <<'EOF'
Usage: submit_stage_d_benchmarks.sh [--dataset ID] [--tool NAME] [--dry-run]

Submit one Slurm job per Stage D baseline cell. The default matrix is
datasets CS-SA01, CS-BS01, CS-PA01, CS-SE01 and tools fastbact, spades,
skesa, unicycler, all at --replicate 1.

  --dataset ID   Submit only this dataset.
  --tool NAME    Submit only this tool.
  --dry-run      Print the sbatch commands and a status table; do not submit.

Each job uses scripts/slurm_benchmark_harness.sbatch and writes to
benchmarks/f26_results/raw/${dataset}_${tool}_r1.
EOF
}

die() {
  echo "submit_stage_d_benchmarks.sh: $*" >&2
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
  if [[ "${text}" =~ ^([0-9]+)(;.*)?$ ]]; then
    printf '%s' "${BASH_REMATCH[1]}"
    return 0
  fi
  if [[ "${text}" =~ Submitted\ batch\ job\ ([0-9]+) ]]; then
    printf '%s' "${BASH_REMATCH[1]}"
    return 0
  fi
  return 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --dataset)
      [[ $# -ge 2 ]] || die "--dataset requires a value"
      [[ -z "${DATASET}" ]] || die "--dataset given more than once"
      DATASET="$2"
      shift 2
      ;;
    --dataset=*)
      [[ -z "${DATASET}" ]] || die "--dataset given more than once"
      DATASET="${1#--dataset=}"
      [[ -n "${DATASET}" ]] || die "--dataset requires a value"
      shift
      ;;
    --tool)
      [[ $# -ge 2 ]] || die "--tool requires a value"
      [[ -z "${TOOL}" ]] || die "--tool given more than once"
      TOOL="$2"
      shift 2
      ;;
    --tool=*)
      [[ -z "${TOOL}" ]] || die "--tool given more than once"
      TOOL="${1#--tool=}"
      [[ -n "${TOOL}" ]] || die "--tool requires a value"
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

if [[ -n "${DATASET}" ]] && ! is_listed "${DATASET}" "${DATASETS[@]}"; then
  die "unknown dataset: ${DATASET}"
fi
if [[ -n "${TOOL}" ]] && ! is_listed "${TOOL}" "${TOOLS[@]}"; then
  die "unknown tool: ${TOOL}"
fi
if [[ ! -f "${MANIFEST}" ]]; then
  die "manifest not found: ${MANIFEST}"
fi
if [[ ! -f "${SBATCH_SCRIPT}" ]]; then
  die "sbatch script not found: ${SBATCH_SCRIPT}"
fi
if [[ "${DRY_RUN}" -eq 0 ]] && ! command -v sbatch >/dev/null 2>&1; then
  die "sbatch is not on PATH"
fi

selected_datasets=()
selected_tools=()
if [[ -n "${DATASET}" ]]; then
  selected_datasets=("${DATASET}")
else
  selected_datasets=("${DATASETS[@]}")
fi
if [[ -n "${TOOL}" ]]; then
  selected_tools=("${TOOL}")
else
  selected_tools=("${TOOLS[@]}")
fi

if [[ "${DRY_RUN}" -eq 0 ]]; then
  mkdir -p "${ROOT}/logs"
fi

command_lines=()
rows=()
failed=0

for dataset in "${selected_datasets[@]}"; do
  for tool in "${selected_tools[@]}"; do
    outdir="${ROOT}/benchmarks/f26_results/raw/${dataset}_${tool}_r${REPLICATE}"
    cmd=(
      sbatch
      --parsable
      --job-name="f26d_${dataset}_${tool}"
      --chdir="${ROOT}"
      "${SBATCH_SCRIPT}"
      --manifest "${MANIFEST}"
      --dataset "${dataset}"
      --tool "${tool}"
      --replicate "${REPLICATE}"
      --outdir "${outdir}"
    )
    if [[ "${DRY_RUN}" -eq 1 ]]; then
      command_lines+=("$(shell_join "${cmd[@]}")")
      rows+=("$(printf 'DRY_RUN\t%s\t%s\t%s\t-\t%s' "${dataset}" "${tool}" "${REPLICATE}" "${outdir}")")
      continue
    fi
    set +e
    submission=$("${cmd[@]}" 2>&1)
    submit_status=$?
    set -e
    if [[ "${submit_status}" -ne 0 ]]; then
      echo "submit_stage_d_benchmarks.sh: sbatch failed for ${dataset} ${tool}: ${submission}" >&2
      rows+=("$(printf 'FAILED\t%s\t%s\t%s\t-\t%s' "${dataset}" "${tool}" "${REPLICATE}" "${outdir}")")
      failed=1
      break
    fi
    if ! job_id=$(parse_job_id "${submission}"); then
      echo "submit_stage_d_benchmarks.sh: unrecognized sbatch output for ${dataset} ${tool}: ${submission}" >&2
      rows+=("$(printf 'FAILED\t%s\t%s\t%s\t-\t%s' "${dataset}" "${tool}" "${REPLICATE}" "${outdir}")")
      failed=1
      break
    fi
    rows+=("$(printf 'SUBMITTED\t%s\t%s\t%s\t%s\t%s' "${dataset}" "${tool}" "${REPLICATE}" "${job_id}" "${outdir}")")
  done
  if [[ "${failed}" -ne 0 ]]; then
    break
  fi
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
