"""Smoke test for the freshly-built openbabel wheel.

Runs in cibuildwheel's test phase. Exercises just enough of the API to
prove that:
- the SWIG _openbabel extension loads,
- BABEL_DATADIR/BABEL_LIBDIR are correctly resolved by __init__.py so
  the plugin format machinery is registered,
- a trivial format conversion (SMILES -> InChI) round-trips.
"""

from __future__ import annotations

import sys


def main() -> int:
    from openbabel import openbabel as ob

    conv = ob.OBConversion()
    if not conv.SetInAndOutFormats("smi", "inchi"):
        print("Failed to find SMILES/InChI format plugins", file=sys.stderr)
        return 1

    mol = ob.OBMol()
    if not conv.ReadString(mol, "c1ccccc1"):
        print("Failed to parse benzene SMILES", file=sys.stderr)
        return 1

    inchi = conv.WriteString(mol).strip()
    expected = "InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H"
    if inchi != expected:
        print(f"Unexpected InChI: {inchi!r} (expected {expected!r})", file=sys.stderr)
        return 1

    print(f"OK: benzene SMILES -> {inchi}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
