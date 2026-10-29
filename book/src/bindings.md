# Python and the browser

## Python (`tpt-chem-py`)

The Python module is built with [maturin](https://www.maturin.rs/) from
`crates/tpt-chem-py`. It is not part of the Cargo workspace, because the
extension links against libpython and would break `cargo test` on Windows.

```sh
cd crates/tpt-chem-py
maturin build --release
pip install target/wheels/tpt_chem-*.whl
```

```python
import tpt_chem

tpt_chem.rhf_energy_xyz(open("water.xyz").read())   # Hartree
tpt_chem.canonical_smiles("OC(C)=O")                # "CC(=O)O"
tpt_chem.formula_from_smiles("CC(=O)O")             # "C2H4O2"
```

## Browser and Node (`tpt-chem-wasm`)

`crates/tpt-chem-wasm` compiles to WebAssembly with wasm-pack:

```sh
wasm-pack build crates/tpt-chem-wasm --target web --out-dir pkg --release
```

`crates/tpt-chem-wasm/www/index.html` runs RHF/STO-3G on a geometry typed
into a text box, with no server-side code. `crates/tpt-chem-wasm/tests/smoke.mjs`
checks the same functions under Node.
