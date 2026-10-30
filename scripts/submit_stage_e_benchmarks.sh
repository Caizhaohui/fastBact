#!/usr/bin/env bash
# Submit Stage E coverage and Stage G replicate benchmarks.
#
# One qfree job per dataset, tool, and replicate. Concurrent jobs do not
# share an outdir: the harness merges summary.json without a lock.
#   benchmarks/f26_results/raw/${dataset}_${tool}_r${replicate}
# Collate finished runs with scripts/collate_benchmark_results.py.
#
# The sbatch script already requests this partition shape. The command
# line repeats it so a header drift cannot change the submitted job:
#   partition=qfree  cpus-per-task=8  mem=14G  time=03:00:00
#
#   bash scripts/submit_stage_e_benchmarks.sh --self-test
#   bash scripts/submit_stage_e_benchmarks.sh --dry-run --dataset CS-BS01_30x
#   bash scripts/submit_stage_e_benchmarks.sh --species 'Bacillus subtilis' --depth 30
#   bash scripts/submit_stage_e_benchmarks.sh --tool all --replicate 2 --dry-run

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ROOT=$(cd "${SCRIPT_DIR}/.." && pwd)
MANIFEST="${ROOT}/F26_STAGE_E_MANIFEST.tsv"
SBATCH_SCRIPT="${ROOT}/scripts/slurm_benchmark_harness.sbatch"

PARTITION="qfree"
CPUS="8"
MEMORY="14G"
WALLTIME="03:00:00"
TOOLS=(fastbact spades skesa unicycler)

DATASET=""
SPECIES=""
TOOL=""
DEPTH=""
REPLICATE="1"
DRY_RUN=0

LOADED_IDS=()
LOADED_SPECIES=()
SELECTED=()
SELECTED_TOOLS=()
DATASET_MODE="none"

usage() {
  cat <<'EOF'
Usage: submit_stage_e_benchmarks.sh [options]

Submit one Slurm job per Stage E manifest row, tool, and replicate.
Stage G uses the same script with --replicate 2 and --replicate 3.
The default tool is fastbact and the default replicate is 1. With no
dataset, species, or depth filter, every manifest row is submitted.

  --manifest PATH    Dataset manifest. Default: F26_STAGE_E_MANIFEST.tsv
  --dataset ID       Exact id, id prefix, substring, or shell glob.
                     CS-BS01_30x is one row. CS-BS01 is every CS-BS01_*
                     row. R1-001 is R1-001_* and not R1-001-Val_*.
                     SA01 is a substring. '*_30x' and 'R1-001_*' are globs.
  --species NAME     Case-insensitive substring of the species column.
  --organism NAME    Alias of --species.
  --tool NAME        fastbact (default), spades, skesa, unicycler, or all.
  --depth DEPTH      Coverage label: 30, 50, 100, 150, or the same with an x.
  --replicate N      Positive integer. Default: 1.
  --dry-run          Print sbatch commands and the status table; do not submit.
  --self-test        Run built-in parser, manifest, filter, and dry-run tests.
  -h, --help         Show this help.

Each job runs scripts/slurm_benchmark_harness.sbatch on partition qfree
with 8 CPUs, 14G memory, and a 03:00:00 wall clock. The harness arguments
are --manifest, --dataset, --tool, --replicate, and
--outdir benchmarks/f26_results/raw/${dataset}_${tool}_r${replicate}.

The status table columns are:
  status  dataset  tool  replicate  job_id  outdir
EOF
}

die() {
  echo "submit_stage_e_benchmarks.sh: $*" >&2
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

depth_of_id() {
  local id="$1"
  if [[ "${id}" =~ _([0-9]+)x$ ]]; then
    printf '%s' "${BASH_REMATCH[1]}"
    return 0
  fi
  return 1
}

normalize_depth() {
  local raw="$1"
  raw="${raw%[xX]}"
  if [[ ! "${raw}" =~ ^[1-9][0-9]*$ ]]; then
    return 1
  fi
  printf '%s' "${raw}"
}

species_matches() {
  local value="${1,,}"
  local filter="${2,,}"
  [[ "${value}" == *"${filter}"* ]]
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
  local header line id species rest
  local -A seen=()
  LOADED_IDS=()
  LOADED_SPECIES=()
  [[ -f "${path}" ]] || die "manifest not found: ${path}"
  IFS= read -r header < "${path}" || die "manifest is empty: ${path}"
  header="${header%$'\r'}"
  local -a header_cols=()
  IFS=$'\t' read -r -a header_cols <<< "${header}"
  [[ "${header_cols[0]:-}" == "dataset_id" ]] || die "manifest first column must be dataset_id: ${path}"
  [[ "${header_cols[1]:-}" == "species" ]] || die "manifest second column must be species: ${path}"
  while IFS= read -r line || [[ -n "${line}" ]]; do
    line="${line%$'\r'}"
    [[ -n "${line}" ]] || continue
    [[ "${line}" == "${header}" ]] && continue
    id="${line%%$'\t'*}"
    rest="${line#*$'\t'}"
    species="${rest%%$'\t'*}"
    [[ -n "${id}" && -n "${species}" ]] || die "manifest row missing dataset_id or species"
    [[ -z "${seen[${id}]+x}" ]] || die "duplicate dataset_id: ${id}"
    seen["${id}"]=1
    LOADED_IDS+=("${id}")
    LOADED_SPECIES+=("${species}")
  done < "${path}"
  [[ ${#LOADED_IDS[@]} -gt 0 ]] || die "manifest has no dataset rows: ${path}"
}

choose_dataset_mode() {
  local filter="$1"
  local id
  if [[ -z "${filter}" ]]; then
    DATASET_MODE="none"
    return 0
  fi
  if [[ "${filter}" == *'*'* || "${filter}" == *'?'* || "${filter}" == *'['* ]]; then
    DATASET_MODE="glob"
    return 0
  fi
  for id in "${LOADED_IDS[@]}"; do
    if [[ "${id}" == "${filter}" ]]; then
      DATASET_MODE="exact"
      return 0
    fi
  done
  for id in "${LOADED_IDS[@]}"; do
    if [[ "${id}" == "${filter}_"* ]]; then
      DATASET_MODE="prefix"
      return 0
    fi
  done
  DATASET_MODE="substring"
}

dataset_id_matches() {
  local id="$1"
  local filter="$2"
  case "${DATASET_MODE}" in
    none)
      return 0
      ;;
    exact)
      [[ "${id}" == "${filter}" ]]
      ;;
    prefix)
      [[ "${id}" == "${filter}_"* ]]
      ;;
    substring)
      [[ "${id}" == *"${filter}"* ]]
      ;;
    glob)
      [[ "${id}" == ${filter} ]]
      ;;
    *)
      return 1
      ;;
  esac
}

select_rows() {
  local i id species depth
  choose_dataset_mode "${DATASET}"
  SELECTED=()
  for i in "${!LOADED_IDS[@]}"; do
    id="${LOADED_IDS[$i]}"
    species="${LOADED_SPECIES[$i]}"
    if ! dataset_id_matches "${id}" "${DATASET}"; then
      continue
    fi
    if [[ -n "${SPECIES}" ]] && ! species_matches "${species}" "${SPECIES}"; then
      continue
    fi
    if [[ -n "${DEPTH}" ]]; then
      depth="$(depth_of_id "${id}")" || continue
      [[ "${depth}" == "${DEPTH}" ]] || continue
    fi
    SELECTED+=("${i}")
  done
}

resolve_tools() {
  SELECTED_TOOLS=()
  if [[ -z "${TOOL}" ]]; then
    SELECTED_TOOLS=(fastbact)
    return 0
  fi
  if [[ "${TOOL}" == "all" ]]; then
    SELECTED_TOOLS=("${TOOLS[@]}")
    return 0
  fi
  if ! is_listed "${TOOL}" "${TOOLS[@]}"; then
    die "unknown tool: ${TOOL}"
  fi
  SELECTED_TOOLS=("${TOOL}")
}

validate_replicate() {
  [[ "${REPLICATE}" =~ ^[1-9][0-9]*$ ]] || die "replicate must be a positive integer: ${REPLICATE}"
}

parse_args() {
  local species_seen=0
  local organism_seen=0
  local organism=""
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
      --dataset)
        [[ $# -ge 2 ]] || die "--dataset requires a value"
        [[ -z "${DATASET}" ]] || die "--dataset given more than once"
        require_value "--dataset" "$2"
        DATASET="$2"
        shift 2
        ;;
      --dataset=*)
        [[ -z "${DATASET}" ]] || die "--dataset given more than once"
        DATASET="${1#--dataset=}"
        require_value "--dataset" "${DATASET}"
        shift
        ;;
      --species)
        [[ $# -ge 2 ]] || die "--species requires a value"
        [[ "${species_seen}" -eq 0 ]] || die "--species given more than once"
        require_value "--species" "$2"
        SPECIES="$2"
        species_seen=1
        shift 2
        ;;
      --species=*)
        [[ "${species_seen}" -eq 0 ]] || die "--species given more than once"
        SPECIES="${1#--species=}"
        require_value "--species" "${SPECIES}"
        species_seen=1
        shift
        ;;
      --organism)
        [[ $# -ge 2 ]] || die "--organism requires a value"
        [[ "${organism_seen}" -eq 0 ]] || die "--organism given more than once"
        require_value "--organism" "$2"
        organism="$2"
        organism_seen=1
        shift 2
        ;;
      --organism=*)
        [[ "${organism_seen}" -eq 0 ]] || die "--organism given more than once"
        organism="${1#--organism=}"
        require_value "--organism" "${organism}"
        organism_seen=1
        shift
        ;;
      --tool)
        [[ $# -ge 2 ]] || die "--tool requires a value"
        [[ -z "${TOOL}" ]] || die "--tool given more than once"
        require_value "--tool" "$2"
        TOOL="$2"
        shift 2
        ;;
      --tool=*)
        [[ -z "${TOOL}" ]] || die "--tool given more than once"
        TOOL="${1#--tool=}"
        require_value "--tool" "${TOOL}"
        shift
        ;;
      --depth)
        [[ $# -ge 2 ]] || die "--depth requires a value"
        [[ -z "${DEPTH}" ]] || die "--depth given more than once"
        require_value "--depth" "$2"
        DEPTH="$2"
        shift 2
        ;;
      --depth=*)
        [[ -z "${DEPTH}" ]] || die "--depth given more than once"
        DEPTH="${1#--depth=}"
        require_value "--depth" "${DEPTH}"
        shift
        ;;
      --replicate)
        [[ $# -ge 2 ]] || die "--replicate requires a value"
        if [[ "${REPLICATE_SET:-0}" -eq 1 ]]; then
          die "--replicate given more than once"
        fi
        require_value "--replicate" "$2"
        REPLICATE="$2"
        REPLICATE_SET=1
        shift 2
        ;;
      --replicate=*)
        if [[ "${REPLICATE_SET:-0}" -eq 1 ]]; then
          die "--replicate given more than once"
        fi
        REPLICATE="${1#--replicate=}"
        require_value "--replicate" "${REPLICATE}"
        REPLICATE_SET=1
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
  if [[ "${species_seen}" -eq 1 && "${organism_seen}" -eq 1 ]]; then
    local left="${SPECIES,,}"
    local right="${organism,,}"
    [[ "${left}" == "${right}" ]] || die "--species and --organism disagree"
  elif [[ "${organism_seen}" -eq 1 ]]; then
    SPECIES="${organism}"
  fi
}

row_outdir() {
  local dataset="$1"
  local tool="$2"
  local replicate="$3"
  printf '%s' "${ROOT}/benchmarks/f26_results/raw/${dataset}_${tool}_r${replicate}"
}

build_sbatch_cmd() {
  local dataset="$1"
  local tool="$2"
  local replicate="$3"
  local manifest="$4"
  local outdir="$5"
  SBATCH_CMD=(
    sbatch
    --parsable
    --partition="${PARTITION}"
    --cpus-per-task="${CPUS}"
    --mem="${MEMORY}"
    --time="${WALLTIME}"
    --job-name="f26e_${dataset}_${tool}_r${replicate}"
    --chdir="${ROOT}"
    "${SBATCH_SCRIPT}"
    --manifest "${manifest}"
    --dataset "${dataset}"
    --tool "${tool}"
    --replicate "${replicate}"
    --outdir "${outdir}"
  )
}

submit_selected() {
  local manifest="$1"
  local idx dataset tool outdir submission submit_status job_id
  local failed=0
  local -a command_lines=()
  local -a rows=()
  for idx in "${SELECTED[@]}"; do
    [[ "${failed}" -eq 0 ]] || break
    dataset="${LOADED_IDS[$idx]}"
    for tool in "${SELECTED_TOOLS[@]}"; do
      outdir="$(row_outdir "${dataset}" "${tool}" "${REPLICATE}")"
      build_sbatch_cmd "${dataset}" "${tool}" "${REPLICATE}" "${manifest}" "${outdir}"
      if [[ "${DRY_RUN}" -eq 1 ]]; then
        command_lines+=("$(shell_join "${SBATCH_CMD[@]}")")
        rows+=("$(printf 'DRY_RUN\t%s\t%s\t%s\t-\t%s' "${dataset}" "${tool}" "${REPLICATE}" "${outdir}")")
        continue
      fi
      set +e
      submission="$("${SBATCH_CMD[@]}" 2>&1)"
      submit_status=$?
      set -e
      if [[ "${submit_status}" -ne 0 ]]; then
        echo "submit_stage_e_benchmarks.sh: sbatch failed for ${dataset} ${tool}: ${submission}" >&2
        rows+=("$(printf 'FAILED\t%s\t%s\t%s\t-\t%s' "${dataset}" "${tool}" "${REPLICATE}" "${outdir}")")
        failed=1
        break
      fi
      if ! job_id="$(parse_job_id "${submission}")"; then
        echo "submit_stage_e_benchmarks.sh: unrecognized sbatch output for ${dataset} ${tool}: ${submission}" >&2
        rows+=("$(printf 'FAILED\t%s\t%s\t%s\t-\t%s' "${dataset}" "${tool}" "${REPLICATE}" "${outdir}")")
        failed=1
        break
      fi
      rows+=("$(printf 'SUBMITTED\t%s\t%s\t%s\t%s\t%s' "${dataset}" "${tool}" "${REPLICATE}" "${job_id}" "${outdir}")")
    done
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

main() {
  local normalized
  parse_args "$@"
  validate_replicate
  if [[ -n "${DEPTH}" ]]; then
    if ! normalized="$(normalize_depth "${DEPTH}")"; then
      die "invalid depth: ${DEPTH}"
    fi
    DEPTH="${normalized}"
  fi
  resolve_tools
  MANIFEST="$(resolve_manifest_path "${MANIFEST}")"
  [[ -f "${SBATCH_SCRIPT}" ]] || die "sbatch script not found: ${SBATCH_SCRIPT}"
  load_manifest "${MANIFEST}"
  select_rows
  if [[ ${#SELECTED[@]} -eq 0 ]]; then
    die "no manifest rows matched the filters"
  fi
  if [[ "${DRY_RUN}" -eq 0 ]]; then
    command -v sbatch >/dev/null 2>&1 || die "sbatch is not on PATH"
    mkdir -p "${ROOT}/logs"
  fi
  submit_selected "${MANIFEST}"
}

self_test() {
  local checks=0
  local failures=0
  local tmp fixture bindir out status expected
  tmp="$(mktemp -d)"
  trap 'rm -rf "${tmp}"' EXIT
  fixture="${tmp}/manifest.tsv"
  bindir="${tmp}/bin"
  mkdir -p "${bindir}"
  ln -s "$(command -v dirname)" "${bindir}/dirname"

  cat > "${fixture}" <<'EOF'
dataset_id	species	strain	notes
AAA_30x	Bacillus subtilis	168	stage_e depth=30x
AAA_50x	Bacillus subtilis	168	stage_e depth=50x
BBB_30x	Staphylococcus aureus	MW2	stage_e depth=30x
R1-001_30x	Escherichia coli	K-12	stage_e depth=30x
R1-001-Val_30x	Escherichia coli	BW25113	stage_e depth=30x
CS-SA01_100x	Staphylococcus aureus	MW2	stage_e depth=100x
EOF

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

  run_expect "unknown argument" 2 bash "${BASH_SOURCE[0]}" --not-a-flag
  run_expect "dataset needs a value" 2 bash "${BASH_SOURCE[0]}" --dataset
  run_expect "duplicate dataset" 2 bash "${BASH_SOURCE[0]}" --dataset AAA --dataset BBB
  run_expect "unknown tool" 2 bash "${BASH_SOURCE[0]}" --tool megahit --dry-run
  run_expect "uppercase tool is not an alias" 2 bash "${BASH_SOURCE[0]}" --tool FASTBACT --dry-run
  run_expect "replicate zero" 2 bash "${BASH_SOURCE[0]}" --replicate 0 --dry-run
  run_expect "replicate word" 2 bash "${BASH_SOURCE[0]}" --replicate abc --dry-run
  run_expect "duplicate replicate" 2 bash "${BASH_SOURCE[0]}" --replicate 1 --replicate 2 --dry-run
  run_expect "invalid depth" 2 bash "${BASH_SOURCE[0]}" --depth 30xx --dry-run
  run_expect "species and organism disagree" 2 \
    bash "${BASH_SOURCE[0]}" --species 'Bacillus subtilis' --organism 'Escherichia coli' --dry-run
  run_expect "missing manifest" 2 \
    bash "${BASH_SOURCE[0]}" --manifest "${tmp}/missing.tsv" --dry-run
  run_expect "help" 0 bash "${BASH_SOURCE[0]}" --help
  if [[ "${out}" != *"--self-test"* || "${out}" != *"--organism"* || "${out}" != *"qfree"* ]]; then
    fail "help text is missing required options"
  else
    pass
  fi

  run_expect "no dataset match" 2 \
    bash "${BASH_SOURCE[0]}" --manifest "${fixture}" --dataset aureus --dry-run

  load_manifest "${fixture}"
  if [[ ${#LOADED_IDS[@]} -eq 6 && "${LOADED_IDS[0]}" == "AAA_30x" && "${LOADED_SPECIES[3]}" == "Escherichia coli" ]]; then
    pass
  else
    fail "fixture manifest load: count=${#LOADED_IDS[@]} first=${LOADED_IDS[0]:-}"
  fi

  DATASET="AAA_30x"
  SPECIES=""
  DEPTH=""
  select_rows
  if [[ ${#SELECTED[@]} -eq 1 && "${LOADED_IDS[${SELECTED[0]}]}" == "AAA_30x" && "${DATASET_MODE}" == "exact" ]]; then
    pass
  else
    fail "exact dataset filter"
  fi

  DATASET="AAA"
  select_rows
  if [[ ${#SELECTED[@]} -eq 2 && "${DATASET_MODE}" == "prefix" ]]; then
    pass
  else
    fail "prefix dataset filter mode=${DATASET_MODE} count=${#SELECTED[@]}"
  fi

  DATASET="R1-001"
  select_rows
  if [[ ${#SELECTED[@]} -eq 1 && "${LOADED_IDS[${SELECTED[0]}]}" == "R1-001_30x" && "${DATASET_MODE}" == "prefix" ]]; then
    pass
  else
    fail "R1-001 prefix must exclude R1-001-Val"
  fi

  DATASET="SA01"
  select_rows
  if [[ ${#SELECTED[@]} -eq 1 && "${LOADED_IDS[${SELECTED[0]}]}" == "CS-SA01_100x" && "${DATASET_MODE}" == "substring" ]]; then
    pass
  else
    fail "substring dataset filter"
  fi

  DATASET="*_30x"
  select_rows
  if [[ ${#SELECTED[@]} -eq 4 && "${DATASET_MODE}" == "glob" ]]; then
    pass
  else
    fail "glob depth pattern count=${#SELECTED[@]} mode=${DATASET_MODE}"
  fi

  DATASET="R1-001_*"
  select_rows
  if [[ ${#SELECTED[@]} -eq 1 && "${LOADED_IDS[${SELECTED[0]}]}" == "R1-001_30x" ]]; then
    pass
  else
    fail "glob R1-001_* must exclude the validation dataset"
  fi

  DATASET=""
  SPECIES="bacillus subtilis"
  DEPTH=""
  select_rows
  if [[ ${#SELECTED[@]} -eq 2 ]]; then
    pass
  else
    fail "case-insensitive species filter count=${#SELECTED[@]}"
  fi

  SPECIES="Staphylococcus aureus"
  DEPTH="30"
  select_rows
  if [[ ${#SELECTED[@]} -eq 1 && "${LOADED_IDS[${SELECTED[0]}]}" == "BBB_30x" ]]; then
    pass
  else
    fail "species and depth intersection"
  fi

  if [[ "$(normalize_depth 30x)" == "30" && "$(normalize_depth 150)" == "150" ]]; then
    pass
  else
    fail "depth normalization"
  fi
  if normalize_depth 0 >/dev/null 2>&1; then
    fail "depth zero should be rejected"
  else
    pass
  fi

  if [[ "$(parse_job_id 424242)" == "424242" \
    && "$(parse_job_id '424242;cluster')" == "424242" \
    && "$(parse_job_id 'Submitted batch job 77')" == "77" ]]; then
    pass
  else
    fail "job id parser accepted forms"
  fi
  if parse_job_id "sbatch: error" >/dev/null 2>&1; then
    fail "job id parser accepted an error line"
  else
    pass
  fi

  load_manifest "${ROOT}/F26_STAGE_E_MANIFEST.tsv"
  if [[ ${#LOADED_IDS[@]} -eq 32 ]]; then
    pass
  else
    fail "stage E manifest row count: ${#LOADED_IDS[@]}"
  fi
  local species_name species_count depth_value depth_count id notes_ok
  for species_name in \
    "Escherichia coli" \
    "Staphylococcus aureus" \
    "Bacillus subtilis" \
    "Pseudomonas aeruginosa" \
    "Salmonella enterica"
  do
    DATASET=""
    SPECIES="${species_name}"
    DEPTH=""
    select_rows
    species_count=${#SELECTED[@]}
    case "${species_name}" in
      "Escherichia coli") expected=16 ;;
      *) expected=4 ;;
    esac
    if [[ "${species_count}" -eq "${expected}" ]]; then
      pass
    else
      fail "species ${species_name} count ${species_count}, expected ${expected}"
    fi
  done
  SPECIES=""
  for depth_value in 30 50 100 150; do
    DEPTH="${depth_value}"
    select_rows
    depth_count=${#SELECTED[@]}
    if [[ "${depth_count}" -eq 8 ]]; then
      pass
    else
      fail "depth ${depth_value} count ${depth_count}, expected 8"
    fi
  done
  DATASET="CS-BS01"
  DEPTH=""
  select_rows
  if [[ ${#SELECTED[@]} -eq 4 && "${DATASET_MODE}" == "prefix" ]]; then
    pass
  else
    fail "real CS-BS01 prefix count=${#SELECTED[@]} mode=${DATASET_MODE}"
  fi
  DATASET="CS-BS01_30x"
  SPECIES="Bacillus subtilis"
  DEPTH="30"
  select_rows
  if [[ ${#SELECTED[@]} -eq 1 && "${LOADED_IDS[${SELECTED[0]}]}" == "CS-BS01_30x" ]]; then
    pass
  else
    fail "combined real filters"
  fi
  DATASET=""
  SPECIES="Salmonella enterica"
  DEPTH="150"
  select_rows
  if [[ ${#SELECTED[@]} -eq 1 && "${LOADED_IDS[${SELECTED[0]}]}" == "CS-SE01_150x" ]]; then
    pass
  else
    fail "Salmonella 150x filter"
  fi
  notes_ok=1
  while IFS=$'\t' read -r id _rest; do
    [[ "${id}" == "dataset_id" ]] && continue
    [[ -z "${id}" ]] && continue
    depth_value="$(depth_of_id "${id}")" || notes_ok=0
    if [[ "${_rest}" != *"depth=${depth_value}x"* ]]; then
      notes_ok=0
    fi
  done < "${ROOT}/F26_STAGE_E_MANIFEST.tsv"
  if [[ "${notes_ok}" -eq 1 ]]; then
    pass
  else
    fail "stage E dataset ids and depth notes disagree"
  fi

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

  cat > "${bindir}/sbatch" <<'EOF'
#!/bin/bash
echo "sbatch should not run during --dry-run" >&2
exit 99
EOF
  chmod +x "${bindir}/sbatch"

  run_expect "dry-run does not execute sbatch" 0 \
    env PATH="${bindir}:/usr/bin:/bin" \
    bash "${BASH_SOURCE[0]}" --manifest "${fixture}" --dataset AAA_30x --dry-run
  expected="${ROOT}/benchmarks/f26_results/raw/AAA_30x_fastbact_r1"
  if [[ "${out}" == *"sbatch --parsable --partition=qfree --cpus-per-task=8 --mem=14G --time=03:00:00 "* \
    && "${out}" == *"--manifest ${fixture} --dataset AAA_30x --tool fastbact --replicate 1 --outdir ${expected}"* \
    && "${out}" == *$'\nstatus\tdataset\ttool\treplicate\tjob_id\toutdir\n'DRY_RUN$'\t'AAA_30x$'\t'fastbact$'\t'1$'\t'-$'\t'"${expected}"* \
    && "${out}" != *"spades"* ]]; then
    pass
  else
    fail "dry-run formatting: ${out}"
  fi

  run_expect "default tool and replicate" 0 \
    env PATH="${bindir}:/usr/bin:/bin" \
    bash "${BASH_SOURCE[0]}" --manifest "${fixture}" --dataset BBB_30x --dry-run
  if [[ "${out}" == *"--tool fastbact --replicate 1 --outdir ${ROOT}/benchmarks/f26_results/raw/BBB_30x_fastbact_r1"* ]]; then
    pass
  else
    fail "defaults were not fastbact replicate 1"
  fi

  run_expect "unrelated organism matches nothing" 2 \
    env PATH="${bindir}:/usr/bin:/bin" \
    bash "${BASH_SOURCE[0]}" --manifest "${fixture}" --organism 'pseudomonas aeruginosa' --depth 100x --dry-run
  if [[ "${out}" == *"no manifest rows matched the filters"* ]]; then
    pass
  else
    fail "unrelated organism should match nothing on the fixture, got: ${out}"
  fi

  run_expect "organism staphylococcus depth" 0 \
    env PATH="${bindir}:/usr/bin:/bin" \
    bash "${BASH_SOURCE[0]}" --manifest "${fixture}" \
      --organism 'Staphylococcus aureus' --depth 100 --tool skesa --replicate 3 --dry-run
  expected="${ROOT}/benchmarks/f26_results/raw/CS-SA01_100x_skesa_r3"
  if [[ "${out}" == *"DRY_RUN"*$'\t'"CS-SA01_100x"$'\t'"skesa"$'\t'"3"$'\t'"-"$'\t'"${expected}" \
    && "${out}" == *"--tool skesa --replicate 3 --outdir ${expected}"* ]]; then
    pass
  else
    fail "organism/depth/tool/replicate dry-run: ${out}"
  fi

  run_expect "tool all order" 0 \
    env PATH="${bindir}:/usr/bin:/bin" \
    bash "${BASH_SOURCE[0]}" --manifest "${fixture}" --dataset AAA_50x --tool all --dry-run
  local tools_in_order
  tools_in_order="$(printf '%s\n' "${out}" | awk -F'\t' '$1=="DRY_RUN"{print $3}' | paste -sd, -)"
  if [[ "${tools_in_order}" == "fastbact,spades,skesa,unicycler" ]]; then
    pass
  else
    fail "tool all order: ${tools_in_order}"
  fi

  run_expect "real manifest dry-run cell" 0 \
    env PATH="${bindir}:/usr/bin:/bin" \
    bash "${BASH_SOURCE[0]}" --dataset CS-BS01_30x --dry-run
  expected="${ROOT}/benchmarks/f26_results/raw/CS-BS01_30x_fastbact_r1"
  if [[ "${out}" == *"--manifest ${ROOT}/F26_STAGE_E_MANIFEST.tsv --dataset CS-BS01_30x --tool fastbact --replicate 1 --outdir ${expected}"* \
    && "${out}" == *"DRY_RUN"*$'\t'"CS-BS01_30x"$'\t'"fastbact"$'\t'"1"$'\t'"-"$'\t'"${expected}" ]]; then
    pass
  else
    fail "default manifest dry-run: ${out}"
  fi

  local dry_lines data_lines
  dry_lines="$(printf '%s\n' "${out}" | awk -F'\t' '$1=="DRY_RUN"{c++} END{print c+0}')"
  data_lines="$(printf '%s\n' "${out}" | awk -F'\t' 'NF>=6 && $1!="status"{c++} END{print c+0}')"
  if [[ "${dry_lines}" -eq 1 && "${data_lines}" -eq 1 ]]; then
    pass
  else
    fail "summary table row count dry=${dry_lines} data=${data_lines}"
  fi

  cat > "${bindir}/sbatch" <<'EOF'
#!/bin/bash
joined="$*"
[[ "${joined}" == *"--partition=qfree"* ]] || exit 3
[[ "${joined}" == *"--cpus-per-task=8"* ]] || exit 3
[[ "${joined}" == *"--mem=14G"* ]] || exit 3
[[ "${joined}" == *"--time=03:00:00"* ]] || exit 3
[[ "${joined}" == *"--parsable"* ]] || exit 3
[[ "${joined}" == *"slurm_benchmark_harness.sbatch"* ]] || exit 3
echo "424242"
EOF
  chmod +x "${bindir}/sbatch"
  run_expect "parsed sbatch job id" 0 \
    env PATH="${bindir}:/usr/bin:/bin" \
    bash "${BASH_SOURCE[0]}" --manifest "${fixture}" --dataset AAA_30x --tool spades --replicate 2
  expected="${ROOT}/benchmarks/f26_results/raw/AAA_30x_spades_r2"
  if [[ "${out}" == *"SUBMITTED"*$'\t'"AAA_30x"$'\t'"spades"$'\t'"2"$'\t'"424242"$'\t'"${expected}" \
    && "${out}" != *"sbatch --parsable"* ]]; then
    pass
  else
    fail "submission table: ${out}"
  fi

  cat > "${bindir}/sbatch" <<'EOF'
#!/bin/bash
echo "Submitted batch job 77"
EOF
  chmod +x "${bindir}/sbatch"
  run_expect "classic sbatch message" 0 \
    env PATH="${bindir}:/usr/bin:/bin" \
    bash "${BASH_SOURCE[0]}" --manifest "${fixture}" --dataset BBB_30x
  if [[ "${out}" == *"SUBMITTED"*$'\t'"BBB_30x"$'\t'"fastbact"$'\t'"1"$'\t'"77"$'\t'* ]]; then
    pass
  else
    fail "classic sbatch message parse: ${out}"
  fi

  cat > "${bindir}/sbatch" <<'EOF'
#!/bin/bash
echo "sbatch: allocation denied" >&2
exit 1
EOF
  chmod +x "${bindir}/sbatch"
  run_expect "sbatch failure" 1 \
    env PATH="${bindir}:/usr/bin:/bin" \
    bash "${BASH_SOURCE[0]}" --manifest "${fixture}" --dataset AAA --tool all
  if [[ "${out}" == *"FAILED"*$'\t'"AAA_30x"$'\t'"fastbact"$'\t'"1"$'\t'"-"$'\t'* \
    && "${out}" != *"AAA_50x"* ]]; then
    pass
  else
    fail "failure stops the batch: ${out}"
  fi

  rm -f "${bindir}/sbatch"
  run_expect "sbatch missing from PATH" 2 \
    env PATH="${bindir}" \
    "$(command -v bash)" "${BASH_SOURCE[0]}" --manifest "${fixture}" --dataset AAA_30x
  if [[ "${out}" == *"sbatch is not on PATH"* ]]; then
    pass
  else
    fail "missing sbatch message: ${out}"
  fi

  local dup="${tmp}/duplicate.tsv"
  printf 'dataset_id\tspecies\nAAA_30x\tBacillus subtilis\nAAA_30x\tBacillus subtilis\n' > "${dup}"
  run_expect "duplicate dataset id" 2 \
    bash "${BASH_SOURCE[0]}" --manifest "${dup}" --dry-run

  local bad="${tmp}/bad-header.tsv"
  printf 'id\tspecies\nAAA\tBacillus subtilis\n' > "${bad}"
  run_expect "bad manifest header" 2 \
    bash "${BASH_SOURCE[0]}" --manifest "${bad}" --dry-run

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
