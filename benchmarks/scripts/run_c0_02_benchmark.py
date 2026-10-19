#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = []
# ///
# How to run: sbatch benchmarks/scripts/submit_c0_02r_final_gate.sh SOURCE_ROOT OUTPUT_DIR
from __future__ import annotations

import sys
from c0_02r.comparison import compare_fasta as compare_fasta


if __name__ == "__main__":
    from c0_02r.benchmark import main

    raise SystemExit(main(sys.argv[1:]))
