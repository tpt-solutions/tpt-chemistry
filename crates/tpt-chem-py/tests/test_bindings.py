"""Smoke tests for the tpt_chem extension. Run after `maturin build` with the wheel installed."""

import unittest

import tpt_chem

WATER = "3\nwater\nO 0 0 0.117\nH 0 0.757 -0.469\nH 0 -0.757 -0.469\n"
# Same value as `cargo run -p tpt-chemistry --example xyz_energy` on this geometry.
WATER_REFERENCE = -74.96294664


class Bindings(unittest.TestCase):
    def test_water_energy_matches_native_solver(self):
        self.assertAlmostEqual(tpt_chem.rhf_energy_xyz(WATER), WATER_REFERENCE, places=7)

    def test_malformed_xyz_raises_value_error(self):
        with self.assertRaises(ValueError):
            tpt_chem.rhf_energy_xyz("not an xyz file")

    def test_formula_and_canonical_smiles(self):
        self.assertEqual(tpt_chem.formula_from_smiles("CC(=O)O"), "C2H4O2")
        self.assertEqual(tpt_chem.canonical_smiles("OC(C)=O"), "CC(=O)O")
        self.assertEqual(tpt_chem.canonical_smiles("O=C(O)C"), "CC(=O)O")

    def test_bad_smiles_raises_value_error(self):
        with self.assertRaises(ValueError):
            tpt_chem.formula_from_smiles("C(")


if __name__ == "__main__":
    unittest.main()
