#!/usr/bin/env python3
"""Write OpenMM reference energies for the crate's electrostatics and LJ tests.

OpenMM (MIT) is used only as an independent reference implementation. It is
not a dependency of the crate.

Writes into tests/data/:
  periodic_ions.json   a neutral 8-ion periodic box (Coulomb only, PME)
  isolated_mixed.json  an isolated 6-atom cluster (LJ + Coulomb, no cutoff)

Each JSON file holds the configuration and OpenMM's potential energy (kJ/mol).
Run with OpenMM importable, for example:

    PYTHONPATH=<site-packages with openmm> python tools/openmm_reference.py tests/data
"""

import json
import sys
from pathlib import Path

import numpy as np
import openmm as mm

ONE_4PI_EPS0 = 138.935456  # kJ·nm/(mol·e²), OpenMM's Coulomb constant


def periodic_ions(rng: np.random.Generator):
    box = 10.0  # Å
    charges = np.array([1.0, -1.0] * 4)
    positions = []
    while len(positions) < 8:
        candidate = rng.uniform(0.0, box, 3)
        if all(np.linalg.norm(candidate - p) > 2.2 for p in positions):
            positions.append(candidate)
    positions = np.array(positions)
    sigma = np.full(8, 0.3)  # nm (unused: epsilon is zero)
    epsilon = np.zeros(8)

    system = mm.System()
    nb = mm.NonbondedForce()
    nb.setNonbondedMethod(mm.NonbondedForce.PME)
    nb.setCutoffDistance(0.45)  # nm
    nb.setEwaldErrorTolerance(1e-8)
    system.setDefaultPeriodicBoxVectors(
        mm.Vec3(box / 10, 0, 0), mm.Vec3(0, box / 10, 0), mm.Vec3(0, 0, box / 10)
    )
    for q, s, e in zip(charges, sigma, epsilon):
        system.addParticle(39.95)
        nb.addParticle(q, s, e)
    system.addForce(nb)
    return box, positions, charges, sigma, epsilon, system


def isolated_mixed(rng: np.random.Generator):
    positions = rng.uniform(-3.0, 3.0, (6, 3))
    charges = np.array([0.6, -0.3, -0.3, 0.4, -0.2, -0.2])
    sigma = np.full(6, 0.34)  # nm
    epsilon = np.full(6, 0.997)  # kJ/mol

    system = mm.System()
    nb = mm.NonbondedForce()
    nb.setNonbondedMethod(mm.NonbondedForce.NoCutoff)
    for q, s, e in zip(charges, sigma, epsilon):
        system.addParticle(39.95)
        nb.addParticle(q, s, e)
    system.addForce(nb)
    return None, positions, charges, sigma, epsilon, system


def energy(system, positions_angstrom):
    context = mm.Context(system, mm.VerletIntegrator(0.001), mm.Platform.getPlatformByName("Reference"))
    context.setPositions([mm.Vec3(*p) / 10.0 for p in positions_angstrom])  # Å → nm
    state = context.getState(getEnergy=True)
    return state.getPotentialEnergy().value_in_unit(mm.unit.kilojoule_per_mole)


def main(out_dir: str) -> None:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    for name, builder, seed in [
        ("periodic_ions", periodic_ions, 11),
        ("isolated_mixed", isolated_mixed, 23),
    ]:
        rng = np.random.default_rng(seed)
        box, positions, charges, sigma, epsilon, system = builder(rng)
        record = {
            "box_angstrom": box,
            "positions_angstrom": positions.tolist(),
            "charges": charges.tolist(),
            "sigma_angstrom": (sigma * 10.0).tolist(),
            "epsilon_kj_mol": epsilon.tolist(),
            "openmm_energy_kj_mol": energy(system, positions),
        }
        (out / f"{name}.json").write_text(json.dumps(record, indent=1), encoding="utf-8")
        print(name, record["openmm_energy_kj_mol"])


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "tests/data")
