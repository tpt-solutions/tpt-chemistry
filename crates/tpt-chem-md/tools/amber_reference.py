#!/usr/bin/env python3
"""Write N-methylacetamide (ACE–NME caps, AMBER ff14SB) and OpenMM's reference values.

OpenMM (MIT) and its amber14 parameter files are used as an independent reference.
Parameters come from OpenMM's amber14-all.xml. The output is converted to Å and
kJ/mol so the Rust side needs no unit handling.

Writes tests/data/nma_amber.json: atoms, bonded terms, nonbonded terms, OpenMM's
potential energy (kJ/mol) and per-atom forces (kJ/mol/Å).

    PYTHONPATH=<site-packages with openmm> python tools/amber_reference.py tests/data
"""

import json
import sys
from pathlib import Path

import numpy as np
import openmm as mm
from openmm import app, unit


def methyl_hydrogens(carbon, partner, rng):
    """Three hydrogens tetrahedrally placed on `carbon`, away from `partner`."""
    axis = carbon - partner
    axis /= np.linalg.norm(axis)
    helper = np.array([0.0, 0.0, 1.0]) if abs(axis[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
    u = np.cross(axis, helper)
    u /= np.linalg.norm(u)
    v = np.cross(axis, u)
    theta = np.radians(180.0 - 109.5)
    hs = []
    for k in range(3):
        phi = np.radians(120.0 * k + 30.0)
        direction = np.cos(theta) * axis + np.sin(theta) * (np.cos(phi) * u + np.sin(phi) * v)
        hs.append(carbon + 1.09 * direction)
    return hs


def build():
    # Approximate geometry (Å). Good enough for an energy/force comparison.
    c_carbonyl = np.array([0.0, 0.0, 0.0])
    o_carbonyl = np.array([0.0, 1.23, 0.0])
    c_acetyl = np.array([-1.52, -0.30, 0.0])
    n_amide = np.array([1.35, -0.70, 0.0])
    c_methyl = np.array([2.70, -0.20, 0.0])
    h_amide = n_amide + 1.01 * np.array([0.0, -1.0, 0.0])
    rng = np.random.default_rng(3)
    h_acetyl = methyl_hydrogens(c_acetyl, c_carbonyl, rng)
    h_nme = methyl_hydrogens(c_methyl, n_amide, rng)
    names_positions = [
        ("ACE", "HH31", "H", h_acetyl[0]),
        ("ACE", "CH3", "C", c_acetyl),
        ("ACE", "HH32", "H", h_acetyl[1]),
        ("ACE", "HH33", "H", h_acetyl[2]),
        ("ACE", "C", "C", c_carbonyl),
        ("ACE", "O", "O", o_carbonyl),
        ("NME", "N", "N", n_amide),
        ("NME", "H", "H", h_amide),
        ("NME", "CH3", "C", c_methyl),
        ("NME", "HH31", "H", h_nme[0]),
        ("NME", "HH32", "H", h_nme[1]),
        ("NME", "HH33", "H", h_nme[2]),
    ]
    topology = app.Topology()
    chain = topology.addChain()
    residues = {}
    atoms = {}
    for resname in ("ACE", "NME"):
        residues[resname] = topology.addResidue(resname, chain)
    for resname, name, element, _ in names_positions:
        atoms[(resname, name)] = topology.addAtom(name, app.Element.getBySymbol(element), residues[resname])
    bonds = [
        (("ACE", "HH31"), ("ACE", "CH3")), (("ACE", "CH3"), ("ACE", "HH32")),
        (("ACE", "CH3"), ("ACE", "HH33")), (("ACE", "CH3"), ("ACE", "C")),
        (("ACE", "C"), ("ACE", "O")),
        (("NME", "N"), ("NME", "H")), (("NME", "N"), ("NME", "CH3")),
        (("NME", "CH3"), ("NME", "HH31")), (("NME", "CH3"), ("NME", "HH32")),
        (("NME", "CH3"), ("NME", "HH33")),
        (("ACE", "C"), ("NME", "N")),
    ]
    for a, b in bonds:
        topology.addBond(atoms[a], atoms[b])
    positions_angstrom = [p for _, _, _, p in names_positions]
    return topology, positions_angstrom, [n for _, n, _, _ in names_positions]


def main(out_dir: str) -> None:
    topology, positions_angstrom, names = build()
    ff = app.ForceField("amber14/protein.ff14SB.xml")
    system = ff.createSystem(topology, nonbondedMethod=app.NoCutoff, constraints=None)
    forces = {type(f).__name__: f for f in system.getForces()}
    bond_f = forces["HarmonicBondForce"]
    angle_f = forces["HarmonicAngleForce"]
    torsion_f = forces["PeriodicTorsionForce"]
    nb_f = forces["NonbondedForce"]

    atoms = []
    for i in range(nb_f.getNumParticles()):
        q, sigma_nm, eps = nb_f.getParticleParameters(i)
        atoms.append({
            "name": names[i],
            "mass": system.getParticleMass(i).value_in_unit(unit.dalton),
            "charge": q.value_in_unit(unit.elementary_charge),
            "sigma": sigma_nm.value_in_unit(unit.nanometer) * 10.0,
            "epsilon": eps.value_in_unit(unit.kilojoule_per_mole),
        })
    bonds = []
    for b in range(bond_f.getNumBonds()):
        i, j, r0, k = bond_f.getBondParameters(b)
        # OpenMM: E = ½ k (r − r0)² with k in kJ/mol/nm² and r0 in nm.
        bonds.append({
            "i": i, "j": j,
            "r0": r0.value_in_unit(unit.nanometer) * 10.0,
            "k": k.value_in_unit(unit.kilojoule_per_mole / unit.nanometer**2) / 100.0,
        })
    angles = []
    for a in range(angle_f.getNumAngles()):
        i, j, k, theta0, kk = angle_f.getAngleParameters(a)
        angles.append({
            "i": i, "j": j, "k": k,
            "theta0": theta0.value_in_unit(unit.radian),
            "kk": kk.value_in_unit(unit.kilojoule_per_mole / unit.radian**2),
        })
    torsions = []
    for t in range(torsion_f.getNumTorsions()):
        i, j, k, l, periodicity, phase, kt = torsion_f.getTorsionParameters(t)
        torsions.append({
            "i": i, "j": j, "k": k, "l": l,
            "periodicity": periodicity,
            "phase": phase.value_in_unit(unit.radian),
            "kt": kt.value_in_unit(unit.kilojoule_per_mole),
        })

    context = mm.Context(system, mm.VerletIntegrator(0.001), mm.Platform.getPlatformByName("Reference"))
    context.setPositions([mm.Vec3(*p) / 10.0 for p in positions_angstrom])
    state = context.getState(getEnergy=True, getForces=True)
    energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    # OpenMM forces are kJ/mol/nm; convert to kJ/mol/Å.
    force = [list(f) for f in state.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole / unit.nanometer)]
    force = [[c / 10.0 for c in f] for f in force]

    record = {
        "positions_angstrom": [list(map(float, p)) for p in positions_angstrom],
        "atoms": atoms,
        "bonds": bonds,
        "angles": angles,
        "torsions": torsions,
        "scee": 1.2,
        "scnb": 2.0,
        "openmm_energy_kj_mol": energy,
        "openmm_forces_kj_mol_angstrom": force,
    }
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    (out / "nma_amber.json").write_text(json.dumps(record, indent=1), encoding="utf-8")
    print("wrote nma_amber.json: energy", energy, "bonds", len(bonds), "angles", len(angles), "torsions", len(torsions))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "tests/data")
