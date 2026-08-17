# zgeom

[![CI](https://github.com/N283T/zgeom/actions/workflows/ci.yml/badge.svg)](https://github.com/N283T/zgeom/actions/workflows/ci.yml)

`zgeom` is a small, standalone Zig CLI for reproducible geometry measurements
in biomolecular structures. The first MVP calculates distances, three-atom
angles, signed four-atom dihedrals, and protein backbone `phi`/`psi`/`omega`.

It is an independent member of the ZBioKit family. It has no runtime or package
dependencies and does not depend on sibling ZBioKit repositories.

## MVP scope

Included now:

- PDB and mmCIF input, including transparent `.gz` decompression
- arbitrary 2-, 3-, and 4-atom measurements by author identifiers
- one-row TSV, CSV, or JSON output for arbitrary measurements
- one-row-per-residue TSV, CSV, or JSON output for backbone torsions
- explicit model, alternate-location, insertion-code, and chain handling
- deterministic output and documented exit status

Deliberately deferred until concrete use cases require them:

- side-chain `chi` definitions and residue-specific symmetry handling
- selection expressions, Cartesian products, and batch query files
- trajectories and periodic boundary conditions (the separate `ztraj` project
  already covers trajectory geometry)
- assemblies, symmetry mates, coordinate transforms, and bond inference
- BCIF, MMTF, SDF, GRO, and other formats

The rationale and follow-up decisions are recorded in
[`docs/design-mvp.md`](docs/design-mvp.md).

## Build and test

Requires Zig 0.16.0 or newer.

```bash
zig build -Doptimize=ReleaseFast
zig build test --summary all

# Optional: compare all 135 defined 1CRN backbone torsions with installed PyMOL
zig build oracle-test --summary all
```

The executable is written to `zig-out/bin/zgeom`.
`oracle-test` requires PyMOL 3.x but PyMOL is not a project dependency.

GitHub Actions runs the dependency-free test suite natively on Linux x86_64,
macOS aarch64, and Windows x86_64. It also cross-compiles release binaries for
x86_64/aarch64 Linux, macOS, and Windows. CI downloads Zig 0.16.0 directly from
the official `ziglang.org` release archive and verifies the platform-specific
SHA-256 published in the official download index before extraction; no
third-party Zig setup action is used.

## Quick start

```bash
# Distance in ångström
zgeom distance structure.cif A:42:N A:42:CA

# Angle in degrees (default) or radians
zgeom angle structure.pdb A:42:N A:42:CA A:42:C --unit deg
zgeom angle structure.pdb A:42:N A:42:CA A:42:C --unit rad

# Signed four-atom torsion
zgeom dihedral structure.cif A:42:N A:42:CA A:42:C A:43:N --format json

# phi, psi, and omega for every polymer residue in chain A
zgeom backbone structure.cif.gz --chain A --format tsv
```

Run `zgeom --help` for the compact command reference.

## Atom selection

An atom is specified as:

```text
CHAIN:AUTH_SEQ[^INS_CODE]:ATOM_NAME[@ALTLOC]
```

Examples:

| Selection | Meaning |
| --- | --- |
| `A:42:CA` | chain A, author residue 42, atom CA |
| `A:42^B:N` | insertion code B on author residue 42 |
| `A:42:CB@B` | explicitly select altloc B |

The CLI uses PDB author identifiers. For mmCIF these are `auth_asym_id`,
`auth_seq_id`, `auth_atom_id`, and `auth_comp_id`, falling back to the matching
`label_*` column only if an author column is absent. Blank chain identifiers
are represented by an empty chain field, for example `:1:CA`.

### Alternate locations

- An atom-level `@ALTLOC` is strict: `@B` only matches conformer B.
- `--altloc B` admits conformer B plus atoms with a blank altloc, which allows
  a complete conformer to be assembled from shared atoms.
- For scalar queries without either option, a blank atom site wins. If none is
  blank, the highest-occupancy atom wins; ties prefer `A`, then lexical order.
- For `backbone`, one named conformer is selected per residue and reused for
  `N`, `CA`, and `C`. The greatest mean occupancy across rows belonging to each
  named conformer wins; exact ties prefer `A`, then lexical order. Blank-altloc
  atoms are shared by every conformer and do not affect this score.
- Scalar atom labels and the backbone `altloc` column expose the resolved
  choice. Blank-only residues have a blank `altloc`, including when a global
  `--altloc` was requested, because they have no named conformer.

Selection must resolve to one logical atom. A missing model or atom is an
error; zgeom never silently substitutes another residue or model.

## Scientific definitions and units

- **Distance:** Euclidean Cartesian distance, output in ångström.
- **Angle `a-b-c`:** angle between `a-b` and `c-b`, in `[0, 180]` degrees by
  default.
- **Dihedral `a-b-c-d`:** signed IUPAC torsion from plane `a-b-c` to plane
  `b-c-d`, in `[-180, 180]` degrees by default.
- **Backbone phi:** `C(i-1)-N(i)-CA(i)-C(i)`.
- **Backbone psi:** `N(i)-CA(i)-C(i)-N(i+1)`.
- **Backbone omega:** `CA(i)-C(i)-N(i+1)-CA(i+1)`; reported on residue `i`.

Angles use `--unit deg` (default) or `--unit rad`. Degenerate geometry, such as
a zero-length angle arm or collinear torsion plane, is undefined. An arbitrary
measurement then exits with status 4. Backbone output reports undefined and
missing values as `NA` in TSV/CSV and `null` in JSON.

Backbone neighbors are consecutive residues in file order, in the same chain,
whose selected `C(i)-N(i+1)` distance is at most 2.0 Å. This prevents torsions
from being calculated across ordinary chain breaks or missing peptide links.
Residue numbering gaps alone do not imply a break, and insertion codes remain
distinct residues. A PDB `TER` record is always a hard break even if the chain
identifier is reused and coordinates happen to be close. In mmCIF, different
`label_asym_id` values are also hard segment boundaries.

All atom rows for one protein residue must be contiguous in the atom table.
If a residue reappears after another protein residue in the same model and
segment, `backbone` rejects the input with status 3 instead of silently
splitting or partially scanning that residue.

## Models and input records

Both `ATOM` and `HETATM` records are retained for arbitrary measurements.
Backbone rows are limited to a recognized protein-residue set (standard amino
acids plus ASX, GLX, UNK, MSE, SEC, and PYL), so nucleic acids, waters, and
ordinary ligands do not create `NA` rows. For mmCIF, a populated `label_seq_id`
is also required. Modified residues outside this explicit set are currently
excluded rather than guessed from atom names; extending the set requires a
reviewed fixture and an unambiguous standard-backbone interpretation.
`--model N` selects a PDB `MODEL` or mmCIF `pdbx_PDB_model_num`. Without it,
the numerically lowest model is used (normally model 1). One invocation never
mixes coordinates across models.

Input format is recognized from `.pdb`, `.cif`, `.pdb.gz`, or `.cif.gz`, with a
small content-based fallback for uncompressed input. Coordinates are read as
64-bit floating-point values; the precision available remains limited by the
source file.

## Output and exit status

TSV is the default. CSV uses the same columns, and JSON provides typed numbers
and `null`. Numeric values are printed to six decimal places. Output ordering
follows coordinate-file residue order and is deterministic.
Backbone output includes the selected residue-level `altloc` after
`residue_name`; JSON uses the same field name.

| Status | Meaning |
| --- | --- |
| `0` | success (including a downstream pipe closed normally) |
| `1` | ordinary I/O, allocation, or unexpected failure |
| `2` | command-line or atom-specification error |
| `3` | unsupported/malformed structure (including corrupt gzip), missing model/atom, or empty result |
| `4` | requested scalar geometry is mathematically undefined |

Diagnostics go to stderr; result data goes to stdout.

## Validation

Unit and CLI integration tests cover analytic distance/angle/dihedral cases,
PDB and mmCIF parsing, TSV/CSV/JSON, exit codes, coherent residue-level altloc
selection, non-contiguous residue rejection, insertion codes, multiple models,
blank chain IDs, and explicit PDB `TER` boundaries. Development smoke tests
additionally use read-only data under
`/Users/nagaet/pdb`:

- AlphaFold DB PDB and mmCIF structures for clean-input agreement
- experimental `1crn.cif.gz` for standard crystallographic mmCIF
- experimental multi-model `1d3z.cif.gz` for explicit model selection
- experimental `4hhb.cif.gz` for multi-chain backbone output
- experimental `1bna.cif.gz` to ensure nucleic acids are excluded
- experimental microheterogeneous `1ejg.cif.gz` to verify altloc-specific
  residue names

On an AlphaFold DB structure, distance, angle, and signed dihedral values were
also checked against PyMOL; differences were below `1e-4` in the reported
units.
The committed 1CRN oracle regression compares every defined residue-level
phi/psi/omega value (135 torsions) directly with PyMOL 3.x using a tolerance of
0.001 degree. The 1CRN fixture was converted from the local wwPDB mmCIF archive
with Gemmi 0.7.4 and contains deposited coordinate records only.

## License

MIT. See [`LICENSE`](LICENSE).
