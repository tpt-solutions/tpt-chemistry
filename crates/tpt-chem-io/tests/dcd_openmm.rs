//! DCD files written by OpenMM (MIT licence) must read back through
//! `DcdReader`, and our writer must produce the same frame records.
//!
//! The fixture `data/openmm_periodic_4atoms.dcd` is written by
//! `tools/openmm_dcd_fixture.py` with OpenMM's `DCDFile` from the frames below
//! (Å). The topology is periodic, so each frame carries a unit cell.

use std::io::Cursor;

use tpt_chem_core::vec3::Vec3;
use tpt_chem_io::dcd::{DcdReader, DcdWriter};

const FIXTURE: &[u8] = include_bytes!("data/openmm_periodic_4atoms.dcd");

const FRAMES: [[[f64; 3]; 4]; 3] = [
    [
        [0.0, 0.0, 0.0],
        [1.5, -2.25, 3.0],
        [0.125, 0.25, 0.5],
        [-4.0, 4.0, 2.5],
    ],
    [
        [0.5, 0.0, 0.0],
        [1.0, -2.0, 3.5],
        [0.0, 0.75, 0.5],
        [-3.5, 4.25, 2.0],
    ],
    [
        [1.0, 0.5, -0.5],
        [0.5, -1.75, 3.25],
        [0.25, 1.0, 0.75],
        [-3.0, 4.5, 1.5],
    ],
];

const CELL: [f64; 6] = [10.0, 0.0, 11.0, 0.0, 0.0, 12.0];

/// Bytes before the first frame: the 84-byte header, the title record and the
/// atom-count record. A title record is 4 + 4 + 80·ntitle + 4 bytes.
fn frames_offset(bytes: &[u8]) -> usize {
    let title_len = u32::from_le_bytes([bytes[92], bytes[93], bytes[94], bytes[95]]) as usize;
    92 + 4 + title_len + 4 + 12
}

#[test]
fn header_matches_the_openmm_fixture() {
    let reader = DcdReader::new(Cursor::new(FIXTURE)).unwrap();
    let header = reader.header();
    assert_eq!(header.n_atoms, 4);
    assert_eq!(header.n_frames, 3);
    assert!(header.has_unitcell);
    // OpenMM writes DELTA in CHARMM AKMA units: a 1 ps step is 20.4548 AKMA.
    assert!(
        (header.delta - 20.4548).abs() < 1e-3,
        "delta {}",
        header.delta
    );
}

#[test]
fn coordinates_and_cells_match_openmm() {
    let frames = DcdReader::new(Cursor::new(FIXTURE))
        .unwrap()
        .collect()
        .unwrap();
    assert_eq!(frames.len(), 3);
    for (frame, expected) in frames.iter().zip(&FRAMES) {
        for (got, want) in frame.positions.iter().zip(expected) {
            let want = Vec3::new(want[0], want[1], want[2]);
            assert!(got.dist(want) < 1e-5, "{got:?} vs {want:?}");
        }
        let cell = frame.unitcell.expect("periodic fixture has cells");
        for (a, b) in cell.iter().zip(CELL) {
            assert!((a - b).abs() < 1e-9, "cell {cell:?}");
        }
    }
}

#[test]
fn our_writer_produces_the_same_frame_records_as_openmm() {
    // Same frames and cells through DcdWriter. The title differs (OpenMM uses
    // two titles, we use three), so compare everything after the title.
    // OpenMM stores the third cell length as 12.000000000000002, one ulp above
    // 12; passing that exact value makes the records comparable byte for byte.
    let openmm_cell = [10.0, 0.0, 11.0, 0.0, 0.0, 12.000_000_000_000_002];
    let mut writer =
        DcdWriter::new(Cursor::new(Vec::new()), 4, "compare", 0, 1, 20.4548, true).unwrap();
    for frame in FRAMES {
        let positions: Vec<Vec3> = frame.iter().map(|p| Vec3::new(p[0], p[1], p[2])).collect();
        writer.write_frame(&positions, Some(openmm_cell)).unwrap();
    }
    let ours = writer.finish().unwrap().into_inner();

    let openmm_frames = frames_offset(FIXTURE);
    let our_frames = frames_offset(&ours);
    assert_eq!(
        &ours[our_frames..],
        &FIXTURE[openmm_frames..],
        "frame records differ from OpenMM's"
    );
}
