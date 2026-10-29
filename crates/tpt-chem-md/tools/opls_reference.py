#!/usr/bin/env python3
"""Write n-butane (OPLS-AA) and OpenMM's reference values to tests/data/butane_opls.json.

Parameters come from foyer's oplsaa.xml (MIT licence, https://github.com/mosdef-hub/foyer).
Each term is looked up by OPLS class, so this script contains no hand-entered
numbers. OpenMM (MIT) evaluates the same terms, and its energy and forces are the reference.

Atom types: CH3 = opls_135, CH2 = opls_136, H = opls_140.

    PYTHONPATH=<site-packages with openmm> python tools/opls_reference.py tests/data /path/to/oplsaa.xml
"""

import json
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import openmm as mm
from openmm import app, unit


def load_xml(path):
    root = ET.parse(path).getroot()
    types = {t.get("name"): t.get("class") for t in root.iter("Type")}
    atoms = {a.get("type"): a for a in root.iter("Atom")}
    bonds = {
        frozenset((b.get("class1"), b.get("class2"))): b for b in root.iter("Bond")
    }
    angles = {}
    for a in root.iter("Angle"):
        angles[(a.get("class1"), a.get("class2"), a.get("class3"))] = a
    torsions = {}
    for t in root.iter("Proper"):
        torsions[(t.get("class1"), t.get("class2"), t.get("class3"), t.get("class4"))] = t
    return types, atoms, bonds, angles, torsions


def first_present(table, *keys):
    # Element objects with no children are falsy, so test for None, not truthiness.
    for key in keys:
        if key in table:
            return table[key]
    return None


def find_angle(angles, c1, c2, c3):
    return first_present(angles, (c1, c2, c3), (c3, c2, c1))


def find_torsion(torsions, c):
    return first_present(torsions, tuple(c), tuple(reversed(c)))


def methyl_hydrogens(carbon, partner):
    axis = carbon - partner
    axis /= np.linalg.norm(axis)
    helper = np.array([0.0, 0.0, 1.0]) if abs(axis[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
    u = np.cross(axis, helper)
    u /= np.linalg.norm(u)
    v = np.cross(axis, u)
    theta = np.radians(180.0 - 109.5)
    return [
        carbon + 1.09 * (np.cos(theta) * axis + np.sin(theta) * (np.cos(np.radians(120 * k + 30)) * u + np.sin(np.radians(120 * k + 30)) * v))
        for k in range(3)
    ]


def build():
    """n-butane: C0–C1–C2–C3 in a zigzag, hydrogens placed approximately (Å)."""
    c = [np.array([0.0, 0.0, 0.0]), np.array([1.27, 0.78, 0.0]), np.array([2.54, 0.0, 0.0]), np.array([3.81, 0.78, 0.0])]
    elements = ["C", "C", "C", "C"]
    types = ["opls_135", "opls_136", "opls_136", "opls_135"]
    positions = list(c)
    bonds = [(0, 1), (1, 2), (2, 3)]
    for h in methyl_hydrogens(c[0], c[1]):
        positions.append(h); elements.append("H"); types.append("opls_140"); bonds.append((0, len(positions) - 1))
    for h in methyl_hydrogens(c[3], c[2]):
        positions.append(h); elements.append("H"); types.append("opls_140"); bonds.append((3, len(positions) - 1))
    # CH2 hydrogens: tetrahedral about each internal carbon.
    for centre, (a, b) in ((1, (0, 2)), (2, (1, 3))):
        to_a = (c[a] - c[centre]) / np.linalg.norm(c[a] - c[centre])
        to_b = (c[b] - c[centre]) / np.linalg.norm(c[b] - c[centre])
        away = -(to_a + to_b)
        away /= np.linalg.norm(away)
        normal = np.cross(to_a, to_b)
        normal /= np.linalg.norm(normal)
        for sign in (1.0, -1.0):
            direction = np.cos(np.radians(54.75)) * away + sign * np.sin(np.radians(54.75)) * normal
            positions.append(c[centre] + 1.09 * direction)
            elements.append("H")
            types.append("opls_140")
            bonds.append((centre, len(positions) - 1))
    return positions, elements, types, bonds


def main(out_dir: str, xml_path: str) -> None:
    types_to_class, atom_params, bond_params, angle_params, torsion_params = load_xml(xml_path)
    positions, elements, opls_types, bonds = build()
    classes = [types_to_class[t] for t in opls_types]

    topology = app.Topology()
    chain = topology.addChain()
    residue = topology.addResidue("BUT", chain)
    atoms = []
    for el, t in zip(elements, opls_types):
        atoms.append(topology.addAtom(t, app.Element.getBySymbol(el), residue))
    for i, j in bonds:
        topology.addBond(atoms[i], atoms[j])

    system = mm.System()
    bond_f = mm.HarmonicBondForce()
    angle_f = mm.HarmonicAngleForce()
    torsion_f = mm.RBTorsionForce()
    nb = mm.NonbondedForce()
    nb.setNonbondedMethod(mm.NonbondedForce.NoCutoff)
    for t in opls_types:
        system.addParticle(12.011 if t != "opls_140" else 1.008)
        a = atom_params[t]
        nb.addParticle(float(a.get("charge")), float(a.get("sigma")), float(a.get("epsilon")))
    system.addForce(nb)
    nb.createExceptionsFromBonds(bonds, 0.5, 0.5)  # OPLS: coulomb14scale = lj14scale = 0.5

    adjacency = {i: set() for i in range(len(atoms))}
    for i, j in bonds:
        adjacency[i].add(j)
        adjacency[j].add(i)
        b = bond_params[frozenset((classes[i], classes[j]))]
        bond_f.addBond(i, j, float(b.get("length")), float(b.get("k")))
    for j in adjacency:
        for i in adjacency[j]:
            for k in adjacency[j]:
                if i < k:
                    a = find_angle(angle_params, classes[i], classes[j], classes[k])
                    if a is None:
                        raise SystemExit(f"no OPLS angle for classes {classes[i]}-{classes[j]}-{classes[k]}")
                    angle_f.addAngle(i, j, k, float(a.get("angle")), float(a.get("k")))
    for j, k in bonds:
        for i in adjacency[j]:
            if i == k:
                continue
            for l in adjacency[k]:
                if l == j or l == i:
                    continue
                cls = (classes[i], classes[j], classes[k], classes[l])
                t = find_torsion(torsion_params, cls)
                c_coeffs = [float(t.get(f"c{n}")) for n in range(6)]
                torsion_f.addTorsion(i, j, k, l, *c_coeffs)
    system.addForce(bond_f)
    system.addForce(angle_f)
    system.addForce(torsion_f)

    context = mm.Context(system, mm.VerletIntegrator(0.001), mm.Platform.getPlatformByName("Reference"))
    context.setPositions([mm.Vec3(*p) / 10.0 for p in positions])
    state = context.getState(getEnergy=True, getForces=True)
    energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    force = [[c / 10.0 for c in f] for f in state.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole / unit.nanometer)]

    record = {
        "positions_angstrom": [list(map(float, p)) for p in positions],
        "atoms": [
            {
                "type": t,
                "charge": float(atom_params[t].get("charge")),
                "sigma": float(atom_params[t].get("sigma")) * 10.0,
                "epsilon": float(atom_params[t].get("epsilon")),
                "mass": 12.011 if t != "opls_140" else 1.008,
            }
            for t in opls_types
        ],
        "bonds": [
            {"i": i, "j": j, "r0": float(bond_params[frozenset((classes[i], classes[j]))].get("length")) * 10.0,
             "k": float(bond_params[frozenset((classes[i], classes[j]))].get("k")) / 100.0}
            for i, j in bonds
        ],
        "scee": 2.0,
        "scnb": 2.0,
        "openmm_energy_kj_mol": energy,
        "openmm_forces_kj_mol_angstrom": force,
    }
    # Angles and torsions are read back from the OpenMM forces, so the Rust side
    # receives exactly the terms OpenMM evaluated.
    record["angles"] = []
    for n in range(angle_f.getNumAngles()):
        i, j, k, theta0, kk = angle_f.getAngleParameters(n)
        record["angles"].append({"i": i, "j": j, "k": k, "theta0": theta0.value_in_unit(unit.radian),
                                  "kk": kk.value_in_unit(unit.kilojoule_per_mole / unit.radian**2)})
    record["rb_torsions"] = []
    for n in range(torsion_f.getNumTorsions()):
        i, j, k, l, c0, c1, c2, c3, c4, c5 = torsion_f.getTorsionParameters(n)
        record["rb_torsions"].append({"i": i, "j": j, "k": k, "l": l,
                                       "c": [c0.value_in_unit(unit.kilojoule_per_mole), c1.value_in_unit(unit.kilojoule_per_mole),
                                             c2.value_in_unit(unit.kilojoule_per_mole), c3.value_in_unit(unit.kilojoule_per_mole),
                                             c4.value_in_unit(unit.kilojoule_per_mole), c5.value_in_unit(unit.kilojoule_per_mole)]})

    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    (out / "butane_opls.json").write_text(json.dumps(record, indent=1), encoding="utf-8")
    print("wrote butane_opls.json: energy", energy, "atoms", len(atoms), "bonds", len(bonds),
          "angles", len(record["angles"]), "rb torsions", len(record["rb_torsions"]))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
