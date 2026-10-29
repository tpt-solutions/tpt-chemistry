// Node smoke test for the wasm-pack `--target web` build.
// Run after `wasm-pack build crates/tpt-chem-wasm --target web`:
//   node crates/tpt-chem-wasm/tests/smoke.mjs
import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import { fileURLToPath } from "node:url";
import { initSync, rhf_energy_xyz, formula_from_smiles, canonical_smiles } from "../pkg/tpt_chem_wasm.js";

const wasmPath = fileURLToPath(new URL("../pkg/tpt_chem_wasm_bg.wasm", import.meta.url));
initSync({ module: readFileSync(wasmPath) });

const WATER = "3\nwater\nO 0 0 0.117\nH 0 0.757 -0.469\nH 0 -0.757 -0.469\n";
// Native reference: `cargo run -p tpt-chemistry --example xyz_energy` on the same geometry.
const NATIVE_WATER_ENERGY = -74.96294664;
const energy = rhf_energy_xyz(WATER);
assert.equal(typeof energy, "number");
assert.ok(Math.abs(energy - NATIVE_WATER_ENERGY) < 1e-7, `${energy} vs ${NATIVE_WATER_ENERGY}`);

assert.equal(formula_from_smiles("CC(=O)O"), "C2H4O2");
assert.equal(canonical_smiles("OC(C)=O"), "CC(=O)O");
assert.throws(() => rhf_energy_xyz("not an xyz file"));
assert.throws(() => formula_from_smiles("C("));

console.log(`ok: water RHF/STO-3G = ${energy.toFixed(8)} Eh`);
