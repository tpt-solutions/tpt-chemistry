//! CHARMM/NAMD DCD binary trajectories.
//!
//! A DCD file is a sequence of Fortran unformatted records. Each record is a
//! 4-byte length, the payload, and the same 4-byte length again. The layout
//! follows the CHARMM/NAMD writer:
//!
//! - header record (84 bytes): `CORD`, NSET (frames), ISTART, NSAVC, six
//!   words of zeros, DELTA (f32 plus a unit-cell flag word in CHARMM files),
//!   eight words of zeros, and the CHARMM version (24)
//! - title record: a count of 80-byte titles and their text
//! - atom-count record: N
//! - per frame: optionally a unit-cell record (six f64), then X, Y and Z
//!   records of N f32 each
//!
//! Positions are in Å and the unit cell is stored as the six CHARMM values
//! `A, cos γ, B, cos β, cos α, C`, the order CHARMM writes. The reader accepts
//! either byte order, and the writer always writes little-endian.
//!
//! The layout was checked against OpenMM's `DCDFile` writer (MIT licence):
//! [`DcdReader`] reads its output (`tests/dcd_openmm.rs`), and [`DcdWriter`]
//! produces byte-identical frame records for the same frames. Time steps are
//! CHARMM AKMA units, as OpenMM writes them.
//!
//! # Examples
//!
//! ```
//! use std::io::Cursor;
//! use tpt_chem_core::vec3::Vec3;
//! use tpt_chem_io::dcd::{DcdReader, DcdWriter};
//!
//! let mut writer = DcdWriter::new(Cursor::new(Vec::new()), 2, "demo", 1, 1, 0.001, false).unwrap();
//! writer.write_frame(&[Vec3::ZERO, Vec3::new(1.0, 2.0, 3.0)], None).unwrap();
//! let bytes = writer.finish().unwrap().into_inner();
//!
//! let mut reader = DcdReader::new(Cursor::new(bytes)).unwrap();
//! assert_eq!(reader.header().n_atoms, 2);
//! let frame = reader.next_frame().unwrap().unwrap();
//! assert_eq!(frame.positions[1], Vec3::new(1.0, 2.0, 3.0));
//! ```

use std::io::{Read, Seek, SeekFrom, Write};

use tpt_chem_core::vec3::Vec3;

use crate::error::{io_err, IoError, Result};

const FORMAT: &str = "DCD";
/// Record marker value of the header and of the title record's length prefix.
const HEADER_RECORD: u32 = 84;
const TITLE_WIDTH: usize = 80;
const TITLES: usize = 3;
/// CHARMM version word written when DELTA is a float (CHARMM convention).
const CHARMM_VERSION: u32 = 24;

/// Header fields of a DCD file.
#[derive(Clone, Debug, PartialEq)]
pub struct DcdHeader {
    /// Atoms per frame.
    pub n_atoms: usize,
    /// Frames declared in the header (zero when the writer did not finish).
    pub n_frames: u32,
    /// First timestep.
    pub istart: u32,
    /// Timesteps between saved frames.
    pub nsavc: u32,
    /// Integration time step in CHARMM AKMA units (1 AKMA ≈ 48.888 fs).
    pub delta: f64,
    /// Whether each frame carries a unit-cell record.
    pub has_unitcell: bool,
    /// Title text, with trailing padding removed.
    pub title: String,
}

/// One trajectory frame.
#[derive(Clone, Debug, PartialEq)]
pub struct DcdFrame {
    /// Positions (Å).
    pub positions: Vec<Vec3>,
    /// CHARMM unit cell `[A, cos γ, B, cos β, cos α, C]`, when present.
    pub unitcell: Option<[f64; 6]>,
}

/// Streaming DCD writer. Frames are appended as they are written, and the
/// frame count in the header is updated by [`DcdWriter::finish`].
pub struct DcdWriter<W: Write + Seek> {
    inner: W,
    n_atoms: usize,
    with_unitcell: bool,
    frames: u32,
}

impl<W: Write + Seek> DcdWriter<W> {
    /// Write the header and title for an `n_atoms` trajectory.
    ///
    /// `istart` and `nsavc` are the first timestep and the spacing between
    /// saved frames, and `delta` is the integration time step in the file's
    /// time units. `with_unitcell` makes every frame carry a unit cell.
    ///
    /// # Errors
    /// Propagates `std::io` failures.
    pub fn new(
        mut inner: W,
        n_atoms: usize,
        title: &str,
        istart: u32,
        nsavc: u32,
        delta: f32,
        with_unitcell: bool,
    ) -> Result<Self> {
        let mut header = Vec::with_capacity(84);
        header.extend_from_slice(b"CORD");
        push_u32(&mut header, 0); // NSET, rewritten by finish()
        push_u32(&mut header, istart);
        push_u32(&mut header, nsavc);
        for _ in 0..6 {
            push_u32(&mut header, 0);
        }
        header.extend_from_slice(&delta.to_le_bytes());
        push_u32(&mut header, u32::from(with_unitcell));
        for _ in 0..8 {
            push_u32(&mut header, 0);
        }
        push_u32(&mut header, CHARMM_VERSION);
        write_record(&mut inner, &header)?;

        let mut titles = vec![b' '; TITLE_WIDTH * TITLES];
        let bytes = title.as_bytes();
        let len = bytes.len().min(titles.len());
        titles[..len].copy_from_slice(&bytes[..len]);
        let mut title_record = Vec::with_capacity(4 + titles.len());
        push_u32(&mut title_record, TITLES as u32);
        title_record.extend_from_slice(&titles);
        write_record(&mut inner, &title_record)?;

        write_record(&mut inner, &(n_atoms as u32).to_le_bytes())?;
        Ok(DcdWriter {
            inner,
            n_atoms,
            with_unitcell,
            frames: 0,
        })
    }

    /// Append one frame.
    ///
    /// # Errors
    /// [`IoError::Unsupported`] when the position count differs from the atom
    /// count, or when a unit cell is missing from a file that stores one.
    /// Otherwise [`IoError`] from the underlying writer.
    pub fn write_frame(&mut self, positions: &[Vec3], unitcell: Option<[f64; 6]>) -> Result<()> {
        if positions.len() != self.n_atoms {
            return Err(IoError::Unsupported {
                format: FORMAT,
                message: format!(
                    "frame has {} positions but the file has {} atoms",
                    positions.len(),
                    self.n_atoms
                ),
            });
        }
        if self.with_unitcell {
            let cell = unitcell.ok_or_else(|| IoError::Unsupported {
                format: FORMAT,
                message: "this file stores unit cells, so every frame needs one".into(),
            })?;
            let mut record = Vec::with_capacity(48);
            for v in cell {
                record.extend_from_slice(&v.to_le_bytes());
            }
            write_record(&mut self.inner, &record)?;
        }
        for axis in 0..3 {
            let mut record = Vec::with_capacity(4 * self.n_atoms);
            for p in positions {
                let v = [p.x, p.y, p.z][axis] as f32;
                record.extend_from_slice(&v.to_le_bytes());
            }
            write_record(&mut self.inner, &record)?;
        }
        self.frames += 1;
        Ok(())
    }

    /// Frames written so far.
    pub fn frames_written(&self) -> u32 {
        self.frames
    }

    /// Write the final frame count into the header and return the sink.
    ///
    /// # Errors
    /// Propagates `std::io` failures.
    pub fn finish(mut self) -> Result<W> {
        self.inner.seek(SeekFrom::Start(8)).map_err(io_err)?;
        self.inner
            .write_all(&self.frames.to_le_bytes())
            .map_err(io_err)?;
        self.inner.seek(SeekFrom::End(0)).map_err(io_err)?;
        self.inner.flush().map_err(io_err)?;
        Ok(self.inner)
    }
}

/// Streaming DCD reader over any byte source.
pub struct DcdReader<R: Read> {
    inner: R,
    header: DcdHeader,
    big_endian: bool,
}

impl<R: Read> DcdReader<R> {
    /// Read the header and title, and detect the byte order.
    ///
    /// # Errors
    /// [`IoError`] when the header is not a DCD header or the file ends early.
    pub fn new(mut inner: R) -> Result<Self> {
        let mut first = [0u8; 4];
        inner.read_exact(&mut first).map_err(io_err)?;
        let big_endian = if u32::from_le_bytes(first) == HEADER_RECORD {
            false
        } else if u32::from_be_bytes(first) == HEADER_RECORD {
            true
        } else {
            return Err(bad("first record is not the 84-byte DCD header"));
        };

        let mut block = [0u8; 84];
        inner.read_exact(&mut block).map_err(io_err)?;
        if &block[0..4] != b"CORD" {
            return Err(bad("header does not start with the CORD tag"));
        }
        let word = |i: usize| -> u32 {
            let b = [
                block[4 + 4 * i],
                block[5 + 4 * i],
                block[6 + 4 * i],
                block[7 + 4 * i],
            ];
            if big_endian {
                u32::from_be_bytes(b)
            } else {
                u32::from_le_bytes(b)
            }
        };
        let n_frames = word(0);
        let istart = word(1);
        let nsavc = word(2);
        // Words 9–10 hold DELTA: an f32 plus a unit-cell flag in CHARMM files,
        // or one f64 otherwise. The CHARMM version sits in word 19.
        let charmm = word(19) != 0;
        let (delta, has_unitcell) = if charmm {
            (f64::from(f32::from_bits(word(9))), word(10) != 0)
        } else {
            let mut raw = [0u8; 8];
            raw.copy_from_slice(&block[4 + 36..4 + 44]);
            let v = if big_endian {
                f64::from_be_bytes(raw)
            } else {
                f64::from_le_bytes(raw)
            };
            (v, false)
        };
        let end = read_u32_at(&mut inner, big_endian)?;
        if end != HEADER_RECORD {
            return Err(bad("header record end marker is wrong"));
        }

        // Title record: a count of 80-byte titles, then the text.
        let title_len = read_u32_at(&mut inner, big_endian)? as usize;
        if title_len < 4 {
            return Err(bad("title record is too short"));
        }
        let mut title_bytes = vec![0u8; title_len];
        inner.read_exact(&mut title_bytes).map_err(io_err)?;
        let title = String::from_utf8_lossy(&title_bytes[4..])
            .trim_end_matches([' ', '\0'])
            .to_string();
        let title_end = read_u32_at(&mut inner, big_endian)?;
        if title_end as usize != title_len {
            return Err(bad("title record end marker does not match its length"));
        }

        // Atom-count record.
        let len = read_u32_at(&mut inner, big_endian)?;
        if len != 4 {
            return Err(bad("atom-count record must be 4 bytes"));
        }
        let n_atoms = read_u32_at(&mut inner, big_endian)? as usize;
        let len_end = read_u32_at(&mut inner, big_endian)?;
        if len_end != 4 {
            return Err(bad("atom-count record end marker is wrong"));
        }

        Ok(DcdReader {
            inner,
            header: DcdHeader {
                n_atoms,
                n_frames,
                istart,
                nsavc,
                delta,
                has_unitcell,
                title,
            },
            big_endian,
        })
    }

    /// The parsed header.
    pub fn header(&self) -> &DcdHeader {
        &self.header
    }

    /// Read the next frame, or `None` at the end of the file.
    ///
    /// # Errors
    /// [`IoError`] on a truncated or malformed frame.
    pub fn next_frame(&mut self) -> Result<Option<DcdFrame>> {
        let n = self.header.n_atoms;
        let mut marker = [0u8; 4];
        match self.inner.read(&mut marker[..1]) {
            Ok(0) => return Ok(None),
            Ok(_) => {}
            Err(e) => return Err(io_err(e)),
        }
        self.inner.read_exact(&mut marker[1..]).map_err(io_err)?;
        let first = self.decode_u32(marker);

        let unitcell = if self.header.has_unitcell {
            if first != 48 {
                return Err(bad("unit-cell record must be 48 bytes"));
            }
            let mut cell = [0.0f64; 6];
            for c in &mut cell {
                let mut raw = [0u8; 8];
                self.inner.read_exact(&mut raw).map_err(io_err)?;
                *c = if self.big_endian {
                    f64::from_be_bytes(raw)
                } else {
                    f64::from_le_bytes(raw)
                };
            }
            let end = read_u32_at(&mut self.inner, self.big_endian)?;
            if end != 48 {
                return Err(bad("unit-cell record end marker is wrong"));
            }
            Some(cell)
        } else {
            // The first marker belongs to the X record.
            None
        };
        let mut axes = Vec::with_capacity(3);
        let mut pending = if self.header.has_unitcell {
            None
        } else {
            Some(first)
        };
        for _ in 0..3 {
            let len = match pending.take() {
                Some(len) => len,
                None => read_u32_at(&mut self.inner, self.big_endian)?,
            };
            if len as usize != 4 * n {
                return Err(bad(
                    "coordinate record length does not match the atom count",
                ));
            }
            let mut values = Vec::with_capacity(n);
            for _ in 0..n {
                let mut raw = [0u8; 4];
                self.inner.read_exact(&mut raw).map_err(io_err)?;
                let bits = if self.big_endian {
                    u32::from_be_bytes(raw)
                } else {
                    u32::from_le_bytes(raw)
                };
                values.push(f64::from(f32::from_bits(bits)));
            }
            let end = read_u32_at(&mut self.inner, self.big_endian)?;
            if end as usize != 4 * n {
                return Err(bad("coordinate record end marker is wrong"));
            }
            axes.push(values);
        }
        let positions = (0..n)
            .map(|i| Vec3::new(axes[0][i], axes[1][i], axes[2][i]))
            .collect();
        Ok(Some(DcdFrame {
            positions,
            unitcell,
        }))
    }

    /// Collect every remaining frame.
    ///
    /// # Errors
    /// [`IoError`] on a truncated or malformed frame.
    pub fn collect(mut self) -> Result<Vec<DcdFrame>> {
        let mut frames = Vec::new();
        while let Some(frame) = self.next_frame()? {
            frames.push(frame);
        }
        Ok(frames)
    }

    fn decode_u32(&self, raw: [u8; 4]) -> u32 {
        if self.big_endian {
            u32::from_be_bytes(raw)
        } else {
            u32::from_le_bytes(raw)
        }
    }
}

fn bad(message: &str) -> IoError {
    IoError::Parse {
        format: FORMAT,
        line: None,
        message: message.into(),
    }
}

fn push_u32(buf: &mut Vec<u8>, v: u32) {
    buf.extend_from_slice(&v.to_le_bytes());
}

/// Write a Fortran unformatted record: length, payload, length.
fn write_record<W: Write>(w: &mut W, payload: &[u8]) -> Result<()> {
    let len = (payload.len() as u32).to_le_bytes();
    w.write_all(&len).map_err(io_err)?;
    w.write_all(payload).map_err(io_err)?;
    w.write_all(&len).map_err(io_err)?;
    Ok(())
}

fn read_u32_at<R: Read>(r: &mut R, big_endian: bool) -> Result<u32> {
    let mut raw = [0u8; 4];
    r.read_exact(&mut raw).map_err(io_err)?;
    Ok(if big_endian {
        u32::from_be_bytes(raw)
    } else {
        u32::from_le_bytes(raw)
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    fn frames() -> Vec<Vec<Vec3>> {
        vec![
            vec![
                Vec3::ZERO,
                Vec3::new(1.5, -2.25, 3.0),
                Vec3::new(0.1, 0.2, 0.3),
            ],
            vec![
                Vec3::new(4.0, 5.0, 6.0),
                Vec3::new(-1.0, 0.5, 0.0),
                Vec3::new(9.0, 9.0, 9.0),
            ],
        ]
    }

    fn write(with_cell: bool) -> Vec<u8> {
        let mut w = DcdWriter::new(
            Cursor::new(Vec::new()),
            3,
            "test run",
            0,
            10,
            0.002,
            with_cell,
        )
        .unwrap();
        for (k, f) in frames().iter().enumerate() {
            let cell = with_cell.then_some([10.0 + k as f64, 0.0, 10.0, 0.0, 0.0, 10.0]);
            w.write_frame(f, cell).unwrap();
        }
        w.finish().unwrap().into_inner()
    }

    #[test]
    fn round_trip_without_unit_cell() {
        let bytes = write(false);
        let mut r = DcdReader::new(Cursor::new(bytes)).unwrap();
        assert_eq!(r.header().n_atoms, 3);
        assert_eq!(r.header().n_frames, 2);
        assert_eq!(r.header().nsavc, 10);
        assert_eq!(r.header().title, "test run");
        assert!(!r.header().has_unitcell);
        assert!((r.header().delta - 0.002).abs() < 1e-7);
        for expected in frames() {
            let got = r.next_frame().unwrap().unwrap();
            assert!(got.unitcell.is_none());
            for (a, b) in got.positions.iter().zip(&expected) {
                assert!(a.dist(*b) < 1e-6, "{a:?} vs {b:?}");
            }
        }
        assert!(r.next_frame().unwrap().is_none());
    }

    #[test]
    fn round_trip_with_unit_cell() {
        let bytes = write(true);
        let frames_read = DcdReader::new(Cursor::new(bytes))
            .unwrap()
            .collect()
            .unwrap();
        assert_eq!(frames_read.len(), 2);
        assert_eq!(
            frames_read[0].unitcell,
            Some([10.0, 0.0, 10.0, 0.0, 0.0, 10.0])
        );
        assert_eq!(frames_read[1].unitcell.unwrap()[0], 11.0);
        assert!(frames_read[1].positions[2].dist(Vec3::new(9.0, 9.0, 9.0)) < 1e-6);
    }

    #[test]
    fn header_bytes_match_the_charmm_layout() {
        let bytes = write(false);
        assert_eq!(&bytes[0..4], &84u32.to_le_bytes());
        assert_eq!(&bytes[4..8], b"CORD");
        assert_eq!(&bytes[8..12], &2u32.to_le_bytes()); // frames
        assert_eq!(&bytes[12..16], &0u32.to_le_bytes()); // istart
        assert_eq!(&bytes[16..20], &10u32.to_le_bytes()); // nsavc
        assert_eq!(&bytes[8 + 19 * 4..8 + 20 * 4], &24u32.to_le_bytes()); // version (word 19)
        assert_eq!(&bytes[88..92], &84u32.to_le_bytes()); // end of header record
    }

    #[test]
    fn big_endian_files_are_read() {
        // Build the same file with big-endian integers and floats by hand.
        let mut out = Vec::new();
        let be = |v: u32| v.to_be_bytes();
        let mut header = Vec::new();
        header.extend_from_slice(b"CORD");
        header.extend_from_slice(&be(1)); // nset
        header.extend_from_slice(&be(0)); // istart
        header.extend_from_slice(&be(1)); // nsavc
        for _ in 0..6 {
            header.extend_from_slice(&be(0));
        }
        header.extend_from_slice(&0.5f32.to_be_bytes());
        header.extend_from_slice(&be(0)); // no unit cell
        for _ in 0..8 {
            header.extend_from_slice(&be(0));
        }
        header.extend_from_slice(&be(24));
        for part in [be(84).to_vec(), header, be(84).to_vec()] {
            out.extend_from_slice(&part);
        }
        // Title record: count + 240 bytes.
        let mut title = be(3).to_vec();
        title.extend(std::iter::repeat_n(b' ', 240));
        out.extend_from_slice(&be(244));
        out.extend_from_slice(&title);
        out.extend_from_slice(&be(244));
        // One atom, one frame.
        out.extend_from_slice(&be(4));
        out.extend_from_slice(&be(1));
        out.extend_from_slice(&be(4));
        for v in [1.25f32, 2.5, -3.0] {
            out.extend_from_slice(&be(4));
            out.extend_from_slice(&v.to_be_bytes());
            out.extend_from_slice(&be(4));
        }
        let mut r = DcdReader::new(Cursor::new(out)).unwrap();
        assert_eq!(r.header().n_atoms, 1);
        let frame = r.next_frame().unwrap().unwrap();
        assert!(frame.positions[0].dist(Vec3::new(1.25, 2.5, -3.0)) < 1e-6);
    }

    #[test]
    fn wrong_atom_count_is_rejected() {
        let mut w = DcdWriter::new(Cursor::new(Vec::new()), 3, "t", 0, 1, 0.001, false).unwrap();
        assert!(w.write_frame(&[Vec3::ZERO], None).is_err());
    }

    #[test]
    fn missing_unit_cell_is_rejected_when_the_file_stores_one() {
        let mut w = DcdWriter::new(Cursor::new(Vec::new()), 1, "t", 0, 1, 0.001, true).unwrap();
        assert!(w.write_frame(&[Vec3::ZERO], None).is_err());
    }

    #[test]
    fn non_dcd_input_is_rejected() {
        assert!(DcdReader::new(Cursor::new(b"not a dcd file at all".to_vec())).is_err());
    }
}
