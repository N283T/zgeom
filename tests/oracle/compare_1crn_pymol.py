"""Compare every defined 1CRN backbone torsion with PyMOL 3.x."""

import json
import os
from pathlib import Path
import subprocess
import traceback

from pymol import cmd


ROOT = Path(os.environ["ZGEOM_ORACLE_ROOT"])
FIXTURE = ROOT / "tests" / "fixtures" / "1crn.pdb"
ZGEOM = ROOT / "zig-out" / "bin" / "zgeom"
TOLERANCE_DEG = 0.001


def selection(resi: int, atom: str) -> str:
    return f"/crn//A/{resi}/{atom}"


def pymol_dihedral(residues: tuple[int, int, int, int], atoms: tuple[str, str, str, str]) -> float:
    return cmd.get_dihedral(*(selection(resi, atom) for resi, atom in zip(residues, atoms)))


def main() -> None:
    cmd.load(str(FIXTURE), "crn")
    rows = json.loads(subprocess.check_output([str(ZGEOM), "backbone", str(FIXTURE), "--format", "json"], text=True))
    if len(rows) != 46:
        raise AssertionError(f"expected 46 residues, got {len(rows)}")

    comparisons = 0
    max_difference = 0.0
    for row in rows:
        resi = int(row["residue_number"])
        expected = {
            "phi": None if resi == 1 else pymol_dihedral((resi - 1, resi, resi, resi), ("C", "N", "CA", "C")),
            "psi": None if resi == 46 else pymol_dihedral((resi, resi, resi, resi + 1), ("N", "CA", "C", "N")),
            "omega": None if resi == 46 else pymol_dihedral((resi, resi, resi + 1, resi + 1), ("CA", "C", "N", "CA")),
        }
        for name, oracle_value in expected.items():
            value = row[name]
            if oracle_value is None:
                if value is not None:
                    raise AssertionError(f"A:{resi} {name}: expected null, got {value}")
                continue
            if value is None:
                raise AssertionError(f"A:{resi} {name}: unexpectedly null")
            difference = abs(value - oracle_value)
            max_difference = max(max_difference, difference)
            comparisons += 1
            if difference > TOLERANCE_DEG:
                raise AssertionError(
                    f"A:{resi} {name}: zgeom={value:.6f}, PyMOL={oracle_value:.6f}, "
                    f"difference={difference:.6g} deg"
                )
    print(
        f"1CRN oracle regression: {comparisons} torsions agree with "
        f"PyMOL {cmd.get_version()[0]} within {TOLERANCE_DEG} deg "
        f"(max {max_difference:.6g} deg)"
    )


try:
    main()
except BaseException:
    traceback.print_exc()
    cmd.quit(1)
else:
    cmd.quit(0)
