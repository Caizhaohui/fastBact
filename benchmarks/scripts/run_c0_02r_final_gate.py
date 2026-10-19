#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = []
# ///
# How to run: sbatch benchmarks/scripts/execute_c0_02r_benchmark.sh SOURCES_ROOT DATA_ROOT NEW_OUTPUT_DIR
from __future__ import annotations

import sys
from c0_02r.benchmark import main

if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
