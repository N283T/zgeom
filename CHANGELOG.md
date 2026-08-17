# Changelog

All notable changes to zgeom are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html).
Because zgeom is still below 1.0, minor releases may include intentional CLI or
output-schema changes; these are called out explicitly below.

## [Unreleased]

## [0.2.0] - 2026-08-17

### Added

- Residue-level `altloc` fields in TSV, CSV, and JSON backbone output so the
  selected conformer is auditable.
- Stable fixtures for coherent alternate conformers, modified amino acids,
  non-contiguous residue rows, corrupt gzip input, and 1CRN.
- An opt-in `zig build oracle-test` regression comparing all 135 defined 1CRN
  backbone torsions with PyMOL 3.x to within 0.001 degree.
- GitHub Actions CI with native Linux x86_64, macOS aarch64, and Windows x86_64
  tests plus ReleaseFast cross-compilation for x86_64/aarch64 Linux, macOS, and
  Windows.
- A cross-platform CI installer that downloads Zig 0.16.0 only from official
  `ziglang.org` archives and verifies fixed, official SHA-256 values before
  extraction.

### Changed

- Backbone analysis now selects one coherent named conformer per residue using
  greatest mean occupancy, with exact ties preferring `A` and then lexical
  order. Blank-altloc atoms remain shared across conformers.
- Backbone peptide-link checks now use the independently resolved conformer for
  each neighboring residue.
- **Breaking output change:** backbone tables and JSON objects gained an
  `altloc` field. This is why the release is 0.2.0 rather than 0.1.1.

### Fixed

- Prevented nonphysical backbone torsions assembled from atom-by-atom mixtures
  of different alternate conformers.
- Rejected protein residues whose atom rows reappear non-contiguously instead
  of silently splitting or partially scanning them.
- Classified corrupt gzip streams as malformed structure input with exit status
  3 rather than a generic failure.
- Kept Zig source files LF-normalized so formatting checks are stable on
  Windows runners.

## [0.1.0] - 2026-08-17

### Added

- Standalone Zig 0.16 CLI commands for distance, three-atom angle, signed
  four-atom dihedral, and protein backbone phi/psi/omega calculations.
- PDB, mmCIF, and transparent gzip input.
- Author-ID atom selection with model, insertion-code, chain, and alternate-
  location handling.
- TSV, CSV, and JSON output with deterministic ordering and documented exit
  behavior.
- Protein-backbone break detection using explicit PDB `TER`, mmCIF
  `label_asym_id`, chain identity, and a 2.0-angstrom C-N distance cutoff.
- Unit and CLI integration coverage for geometry, parsing, multiple models,
  alternate locations, insertion codes, blank chains, and segment boundaries.

[Unreleased]: https://github.com/N283T/zgeom/compare/v0.2.0...HEAD
[0.2.0]: https://github.com/N283T/zgeom/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/N283T/zgeom/releases/tag/v0.1.0
