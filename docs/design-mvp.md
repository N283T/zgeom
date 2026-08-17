# zgeom MVP design memo

## Investigation summary

The repository initially contained only `AGENTS.md` and a short README. The
following sibling projects were inspected before defining this MVP:

- **zsasa:** current Zig 0.16 build layout, explicit CLI/library separation,
  broad structure I/O, and deterministic scientific output conventions.
- **zdssp:** compact geometry primitives and protein-backbone parsing patterns.
- **zreduce:** author identifiers, insertion codes, altloc/model policy, and
  current Zig 0.16 `std.Io` usage.
- **ztraj:** distances, angles, signed torsions, protein torsion definitions,
  selections, and trajectory I/O.

`ztraj` already has broad trajectory geometry. Depending on it would add
unrelated analysis and I/O surface to this small static-structure tool, so the
MVP remains self-contained. Its public behavior and tests were useful design
references, but zgeom owns its parser, data model, CLI, and scalar geometry.

The local data tree contains clean AlphaFold DB PDB/mmCIF sets plus the PDB
archive as gzipped experimental mmCIF files. That led to including both PDB and
mmCIF plus gzip in the initial input boundary; PDB-only support would not meet
the stated validation workflow.

## Why these four commands

The smallest useful surface is:

1. `distance` for two explicitly named atoms;
2. `angle` for an explicit triplet;
3. `dihedral` for an explicit quartet;
4. `backbone` for phi, psi, and omega together.

A generic four-atom command alone is too cumbersome for residue-wide
Ramachandran analysis. Three separate `phi`, `psi`, and `omega` commands would
repeat parsing and make row correspondence harder. One combined backbone table
is smaller and more useful for a conference workflow.

Side-chain chi angles are deferred. They require a reviewed residue-specific
atom table, naming rules for modified residues, and symmetry conventions; a
partial generic pattern table would look complete while producing ambiguous
scientific results.

## Decisions

### Author IDs over label IDs

Users normally cite PDB chain and residue numbers. The selection surface thus
uses mmCIF author IDs, with label IDs only as column-level fallbacks. Mixing
author chain IDs with label sequence IDs within one selection is not allowed.

### One model per invocation

Defaulting to the lowest model matches common single-model structures while
remaining deterministic for NMR inputs. `--model` is explicit, and no
coordinate tuple may cross models. An eventual all-model mode should emit a
model column rather than average values implicitly.

### Altloc policy

Default selection favors a blank site, otherwise maximum occupancy. This is a
deterministic pragmatic default, not a claim that independently selected
maximum-occupancy atoms form a physical conformer. `--altloc` exists for a
coherent named conformer and still accepts shared blank atoms. Strict atom-level
`@ALT` is provided for inspecting a particular site. The resolved label is
printed for arbitrary measurements.

### Backbone breaks

File adjacency and chain equality are insufficient because missing segments
can leave distant residues adjacent. Residue-number continuity is also
insufficient because legitimate insertion codes and numbering gaps occur. The
MVP therefore requires a `C-N` distance no greater than 2.0 Å. The threshold is
intentionally exposed in documentation but not yet as a CLI option; adding a
flag without a demonstrated use case would enlarge the interface prematurely.

### Missing and degenerate data

Backbone tables retain one row per residue and mark undefined values, which
preserves correspondence for downstream analysis. A scalar query has no useful
row if an atom is missing or geometry is degenerate, so it exits nonzero.

### Output formats

TSV is the human-readable and Unix-friendly default. CSV supports common data
tools. JSON carries numeric types and nulls. JSONL and batch query files are
deferred until a real high-throughput query format is designed.

## Known MVP limitations / next decisions

1. Define whether batch queries need a simple table input or a selection
   expression language; do not introduce both initially.
2. Decide whether multi-model output should mean every model independently or
   include circular statistics across models.
3. Review residue-specific chi definitions and symmetric terminal groups before
   adding chi output.
4. If assembly geometry is requested, distinguish deposited coordinates,
   biological assemblies, and crystallographic symmetry explicitly.
5. Consider a narrow `ztraj` dependency only when trajectory support becomes a
   concrete requirement; static structure I/O should remain independent.
6. CSV currently targets conventional PDB/mmCIF identifiers without commas or
   newlines. If arbitrary quoted mmCIF identifiers become a requirement, add a
   fully RFC 4180-compliant field writer and a reversible atom-spec escaping
   rule together.
7. Protein-backbone classification uses a small recognized residue set
   (including MSE, SEC, and PYL), with mmCIF `label_seq_id` as an additional
   polymer signal. Extend this deliberately when a modified-residue fixture is
   added rather than treating nucleic acids or every ligand as protein.

## Hardening decisions after independent review

- PDB `TER` and mmCIF `label_asym_id` are retained as hard segment boundaries;
  geometry is never calculated across them.
- Backbone atom lookup is bounded by residue atom ranges, making the torsion
  pass linear rather than rescanning the whole structure for every residue.
- Microheterogeneous output uses the component ID of the selected CA (falling
  back to N or C), rather than the first component row at that author site.
- Non-finite coordinates/occupancies and incomplete final atom-site rows are
  rejected as malformed structures with exit status 3.
- CIF reserved words, category names, and tags are matched case-insensitively.
- Torsion collinearity is tested relative to bond-vector scale, so changing
  coordinate units does not change whether geometry is classified degenerate.
