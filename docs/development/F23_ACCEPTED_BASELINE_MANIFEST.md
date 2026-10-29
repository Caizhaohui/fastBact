# fastBact Phase F23 Accepted Baseline Manifest

- **Phase**: F23 (Provenance-Aware Repeat Resolution Engine)
- **Status**: ACCEPTED & IMMUTABLE FREEZE
- **Scientific Contract**: C0-02R: ACCEPTED
- **Freeze Tag**: `f23-accepted-2026-10-29`
- **Tag Object**: `0a74fc496d511526753c365b02cdea0ecb9754ae`
- **Target Peeled Commit**: `68b19f06d22537cf199e8fe124ea6b62bb0bb50a`
- **Remote Branch**: `origin/main` (`68b19f06d22537cf199e8fe124ea6b62bb0bb50a`)
- **Engine Source**: `src/resolve/exspander.rs`
- **Engine SHA-256**: `fd59e8bdc231dbb83468c711f3d87834c51638d75941f587eb0aa1c89298790c`
- **Verified Binary SHA-256**: `c2a5dc9cf2d88195a9b7dc2dca592f61e8601d0a51c73ecaa357fe1dc7189fef`

## Test Verification Signatures
- **Pre-Merge Deterministic Test**: Slurm Job 3062621 (97/97 tests pass, 19 F1-F9b, 42 T26-T39)
- **Release Regression Suite**: Slurm Job 3062677 (204/204 cargo test --release pass across all workspace crates)
- **Release Regression Log SHA-256**: `7b8554f730683e1c5e55bc3dbeab8a34cac097e5bf194b2a58dee95882364c43`

## Baseline QUAST Metrics
- **R1-001** (*E. coli* K-12):
  - Misassemblies: 5
  - N50: 28,137 bp
  - NGA50: 27,031 bp
  - Genome Fraction: 98.243%
  - AUD-01 2.1 Mb Chimera: ELIMINATED
- **R1-003** (*S. aureus*):
  - Misassemblies: 0
  - N50: 17,923 bp
  - NGA50: 17,509 bp
  - Indels per 100 kbp: 0.00
  - Genome Fraction: 99.390%
