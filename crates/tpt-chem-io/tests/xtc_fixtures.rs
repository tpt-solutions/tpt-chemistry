//! XTC fixtures in `tests/data/xtc`.
//!
//! - `gas_5`: written and decoded by mdtraj (1.11.2) once. Five atoms take the
//!   uncompressed path, and the oracle is mdtraj's decode.
//! - `ours_*`: written by `XtcWriter` from formula positions. The expected
//!   values are those positions rounded to 0.001 nm, so no decoder produced
//!   them. The files also decode with molstar's MIT parser (checked once by
//!   hand, see the README).
//!
//! To rewrite the `ours_*` files, run
//! `cargo test -p tpt-chem-io --test xtc_fixtures -- --ignored regenerate_our_fixtures`.

use std::fs;
use std::io::Cursor;
use std::path::PathBuf;

use tpt_chem_core::vec3::Vec3;
use tpt_chem_io::xtc::{XtcFrame, XtcReader, XtcWriter};

const TOLERANCE_ANGSTROM: f64 = 1e-5;

struct Oracle {
    frames: Vec<OracleFrame>,
}

struct OracleFrame {
    step: i32,
    time: f32,
    cell_nm: [f64; 9],
    positions_nm: Vec<[f64; 3]>,
}

fn data_path(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/xtc")
        .join(name)
}

fn parse_oracle(text: &str) -> Oracle {
    let mut lines = text.lines();
    let header: Vec<&str> = lines.next().expect("header").split_whitespace().collect();
    let n_frames: usize = header[1].parse().unwrap();
    let n_atoms: usize = header[3].parse().unwrap();
    let mut frames = Vec::with_capacity(n_frames);
    for _ in 0..n_frames {
        let t: Vec<&str> = lines
            .next()
            .expect("frame line")
            .split_whitespace()
            .collect();
        let step: i32 = t[1].parse().unwrap();
        let time: f32 = t[3].parse().unwrap();
        let mut cell_nm = [0.0; 9];
        for (k, v) in cell_nm.iter_mut().enumerate() {
            *v = t[5 + k].parse().unwrap();
        }
        let positions_nm = (0..n_atoms)
            .map(|_| {
                let p: Vec<f64> = lines
                    .next()
                    .expect("atom line")
                    .split_whitespace()
                    .map(|s| s.parse().unwrap())
                    .collect();
                [p[0], p[1], p[2]]
            })
            .collect();
        frames.push(OracleFrame {
            step,
            time,
            cell_nm,
            positions_nm,
        });
    }
    Oracle { frames }
}

fn check_against_oracle(name: &str) {
    let bytes = fs::read(data_path(&format!("{name}.xtc"))).expect("fixture");
    let oracle =
        parse_oracle(&fs::read_to_string(data_path(&format!("{name}.oracle.txt"))).unwrap());
    let frames = XtcReader::new(Cursor::new(bytes)).collect().unwrap();
    assert_eq!(frames.len(), oracle.frames.len(), "{name}: frame count");

    for (k, (got, want)) in frames.iter().zip(&oracle.frames).enumerate() {
        assert_eq!(got.step, want.step, "{name} frame {k}: step");
        assert_eq!(got.time, want.time, "{name} frame {k}: time");
        for r in 0..3 {
            for c in 0..3 {
                let expect = want.cell_nm[3 * r + c] * 10.0;
                assert!(
                    (got.cell[r][c] - expect).abs() < TOLERANCE_ANGSTROM,
                    "{name} frame {k}: cell[{r}][{c}] {} vs {expect}",
                    got.cell[r][c]
                );
            }
        }
        assert_eq!(got.positions.len(), want.positions_nm.len());
        for (i, (p, q)) in got.positions.iter().zip(&want.positions_nm).enumerate() {
            let err = [
                (p.x - q[0] * 10.0).abs(),
                (p.y - q[1] * 10.0).abs(),
                (p.z - q[2] * 10.0).abs(),
            ];
            assert!(
                err.iter().all(|&e| e < TOLERANCE_ANGSTROM),
                "{name} frame {k} atom {i}: {p:?} vs {q:?} (nm) err {err:?}"
            );
        }
    }
}

/// Deterministic pseudo-random number in `[0, 1)` (a 64-bit LCG).
fn unit(seed: u64) -> f64 {
    let x = seed
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    (x >> 11) as f64 / (1u64 << 53) as f64
}

/// Random gas in a cubic box of side `box_a` (Å).
fn gas_frame(n: usize, f: usize, box_a: f64) -> Vec<Vec3> {
    (0..n)
        .map(|i| {
            let s = (f * 100_003 + i * 1_009) as u64;
            Vec3::new(
                unit(s * 3) * box_a,
                unit(s * 3 + 1) * box_a,
                unit(s * 3 + 2) * box_a,
            )
        })
        .collect()
}

/// Water-like O–H–H triplets on a jittered lattice in a cubic box (Å).
fn water_frame(n_mol: usize, f: usize, box_a: f64) -> Vec<Vec3> {
    let side = (n_mol as f64).cbrt().ceil() as usize;
    let mut out = Vec::with_capacity(3 * n_mol);
    let step = box_a / side as f64;
    for m in 0..n_mol {
        let (ix, iy, iz) = (m % side, (m / side) % side, m / (side * side));
        let c = Vec3::new(
            (ix as f64 + 0.5) * step,
            (iy as f64 + 0.5) * step,
            (iz as f64 + 0.5) * step,
        );
        let jitter = |k: usize| (unit((f * 1_000_003 + m * 97 + k) as u64) - 0.5) * 0.05;
        out.push(c + Vec3::new(jitter(0), jitter(1), jitter(2)));
        out.push(c + Vec3::new(0.96 + jitter(3), jitter(4), jitter(5)));
        out.push(c + Vec3::new(-0.24 + jitter(6), 0.93 + jitter(7), jitter(8)));
    }
    out
}

/// The position the file stores: quantized to 0.001 nm, as the writer does.
fn quantized_nm(p: Vec3) -> [f64; 3] {
    [p.x, p.y, p.z].map(|v| f64::from(((v / 10.0) as f32 * 1000.0).round()) / 1000.0)
}

/// Writes `<name>.xtc` and `<name>.oracle.txt` for the given frames.
fn write_fixture(name: &str, frames: &[XtcFrame]) {
    let mut w = XtcWriter::new(Cursor::new(Vec::new()));
    for f in frames {
        w.write_frame(f).unwrap();
    }
    fs::write(
        data_path(&format!("{name}.xtc")),
        w.into_inner().into_inner(),
    )
    .unwrap();

    let n = frames[0].positions.len();
    let mut text = format!("frames {} atoms {}\n", frames.len(), n);
    for f in frames {
        let cell: Vec<String> = f
            .cell
            .iter()
            .flatten()
            .map(|v| format!("{:?}", v / 10.0))
            .collect();
        text.push_str(&format!(
            "step {} time {:?} box {}\n",
            f.step,
            f.time,
            cell.join(" ")
        ));
        for p in &f.positions {
            let q = quantized_nm(*p);
            text.push_str(&format!("{:?} {:?} {:?}\n", q[0], q[1], q[2]));
        }
    }
    fs::write(data_path(&format!("{name}.oracle.txt")), text).unwrap();
}

fn frames_of(
    n_atoms: usize,
    n_frames: usize,
    box_a: f64,
    positions: impl Fn(usize, usize) -> Vec<Vec3>,
) -> Vec<XtcFrame> {
    (0..n_frames)
        .map(|f| XtcFrame {
            step: 10 * f as i32,
            time: 0.02 * f as f32,
            cell: [[box_a, 0.0, 0.0], [0.0, box_a, 0.0], [0.0, 0.0, box_a]],
            positions: positions(n_atoms, f),
        })
        .collect()
}

#[test]
#[ignore = "rewrites the ours_* fixtures; run on purpose with --ignored"]
fn regenerate_our_fixtures() {
    write_fixture(
        "ours_gas_100",
        &frames_of(100, 3, 20.0, |n, f| gas_frame(n, f, 20.0)),
    );
    write_fixture(
        "ours_water_3000",
        &frames_of(1000, 2, 30.0, |n, f| water_frame(n, f, 30.0)),
    );
}

#[test]
fn five_atoms_uncompressed_matches_mdtraj() {
    check_against_oracle("gas_5");
}

#[test]
fn our_gas_fixture_matches_formula() {
    check_against_oracle("ours_gas_100");
}

#[test]
fn our_water_fixture_matches_formula() {
    check_against_oracle("ours_water_3000");
}

#[test]
fn our_writer_output_round_trips_at_default_precision() {
    let frames = frames_of(50, 3, 25.0, |n, f| gas_frame(n, f, 25.0));
    let mut w = XtcWriter::new(Cursor::new(Vec::new()));
    for f in &frames {
        w.write_frame(f).unwrap();
    }
    let bytes = w.into_inner().into_inner();
    let back = XtcReader::new(Cursor::new(bytes)).collect().unwrap();
    assert_eq!(back.len(), frames.len());
    for (a, b) in back.iter().zip(&frames) {
        assert_eq!(a.step, b.step);
        for (p, q) in a.positions.iter().zip(&b.positions) {
            // Half a 0.001 nm unit is 0.005 Å.
            assert!(
                (p.x - q.x).abs() <= 0.0051
                    && (p.y - q.y).abs() <= 0.0051
                    && (p.z - q.z).abs() <= 0.0051
            );
        }
    }
}
