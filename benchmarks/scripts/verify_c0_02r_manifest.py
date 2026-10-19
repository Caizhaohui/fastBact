#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = []
# ///
# How to run: srun -p qfree -w vnode1 -c 8 --mem=12G uv run verify_c0_02r_manifest.py MANIFEST
from __future__ import annotations

from pathlib import Path
import sys
from c0_02r.verification import verify_manifest

if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("Usage: verify_c0_02r_manifest.py MANIFEST")
    files, commands = verify_manifest(Path(sys.argv[1]))
    print(f"INDEPENDENT_INTEGRITY_PASS files={files} commands={commands}")
