# XTC fixtures

- `gas_5.xtc` and `gas_5.oracle.txt`: five atoms, two frames, uncompressed
  path (≤ 9 atoms). Written and decoded once by mdtraj 1.11.2. The oracle
  holds mdtraj's decoded coordinates (nm).
- `ours_gas_100.xtc`: 100 atoms, 3 frames, 20 Å box. Random positions from a
  64-bit LCG (see `gas_frame` in `tests/xtc_fixtures.rs`).
- `ours_water_3000.xtc`: 1000 water-like triplets, 2 frames, 30 Å box.

The `ours_*` files are written by `XtcWriter` (precision 1000). Their oracles
hold the formula positions rounded to 0.001 nm, not a decoder's output. The
files also decode exactly with molstar's MIT-licensed XTC parser
(`src/mol-io/reader/xtc/parser.ts`, Copyright (c) 2020 mol* contributors),
checked once by a line-for-line Python port of its decoding loop. That port
is not committed.

Regenerate the `ours_*` files with
`cargo test -p tpt-chem-io --test xtc_fixtures -- --ignored regenerate_our_fixtures`.
