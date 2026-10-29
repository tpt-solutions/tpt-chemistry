#!/usr/bin/env python3
"""Write tests/data/openmm_4atoms.dcd with OpenMM's DCDFile writer (MIT licence).

OpenMM positions are in nm and the unit cell in nm. The values in the file are
the same frames as tests/dcd_fixture.rs, scaled by 0.1 to nm.

    PYTHONPATH=<site-packages with openmm> python tools/openmm_dcd_fixture.py tests/data
"""

import sys
from pathlib import Path

import numpy as np
import openmm.app as app
from openmm import Vec3, unit

FRAMES_ANGSTROM = np.array(
    [
        [[0.0, 0.0, 0.0], [1.5, -2.25, 3.0], [0.125, 0.25, 0.5], [-4.0, 4.0, 2.5]],
        [[0.5, 0.0, 0.0], [1.0, -2.0, 3.5], [0.0, 0.75, 0.5], [-3.5, 4.25, 2.0]],
        [[1.0, 0.5, -0.5], [0.5, -1.75, 3.25], [0.25, 1.0, 0.75], [-3.0, 4.5, 1.5]],
    ]
)
BOX_ANGSTROM = (10.0, 11.0, 12.0)


def main(out_dir: str) -> None:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    topology = app.Topology()
    chain = topology.addChain()
    residue = topology.addResidue("MOL", chain)
    for _ in range(4):
        topology.addAtom("X", app.Element.getByAtomicNumber(6), residue)
    # A periodic topology makes OpenMM write unit-cell records for each frame.
    topology.setPeriodicBoxVectors(
        [Vec3(BOX_ANGSTROM[0] * 0.1, 0, 0), Vec3(0, BOX_ANGSTROM[1] * 0.1, 0), Vec3(0, 0, BOX_ANGSTROM[2] * 0.1)]
        * unit.nanometer
    )
    path = out / "openmm_4atoms.dcd"
    with open(path, "wb") as handle:
        dcd = app.DCDFile(handle, topology, 1.0 * unit.picoseconds)
        for frame in FRAMES_ANGSTROM:
            positions = [Vec3(*p) * 0.1 for p in frame] * unit.nanometer
            dcd.writeModel(positions)
    print("wrote", path)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "tests/data")
