//! GROMACS XTC trajectories: XDR big-endian frames with compressed
//! coordinates.
//!
//! Each frame is `magic 1995`, the atom count, the step, the time (ps), the
//! 3×3 box (nm), and the atom count again. Small systems (≤ 9 atoms) store
//! `x y z` floats directly. Larger ones store the coordinates as integers
//! scaled by `precision` (1000 for 0.001 nm), with the integer ranges and the
//! bit-packed stream described below.
//!
//! Compressed stream: each atom is coded either absolutely or as a *run*. A
//! run starts with an absolute anchor atom, then `K` small deltas (each
//! coded in `smallidx` bits per component, offset by half the range). The
//! first delta is taken from the anchor and the output order puts that atom
//! first, which suits water-like triplets. The delta width adapts by ±1 per
//! run through `is_smaller`, as in the reference format.
//!
//! Units follow the crate convention: [`XtcFrame`] holds Ångström, and the
//! file stores nm.
//!
//! The format was implemented from its published layout. The compressed
//! coordinates are checked against committed fixtures decoded by mdtraj once
//! (see `tests/data/xtc/README.md` and `tests/xtc_mdtraj.rs`). Frames with a
//! coordinate range that needs more
//! than 24 bits per component are rejected.
//!
//! # Examples
//!
//! ```
//! use std::io::Cursor;
//! use tpt_chem_core::vec3::Vec3;
//! use tpt_chem_io::xtc::{XtcFrame, XtcReader, XtcWriter};
//!
//! let frame = XtcFrame {
//!     step: 10,
//!     time: 0.5,
//!     cell: [[20.0, 0.0, 0.0], [0.0, 20.0, 0.0], [0.0, 0.0, 20.0]],
//!     positions: vec![Vec3::new(1.0, 2.0, 3.0), Vec3::new(1.2, 2.1, 3.3)],
//! };
//! let mut writer = XtcWriter::new(Cursor::new(Vec::new()));
//! writer.write_frame(&frame).unwrap();
//! let bytes = writer.into_inner().into_inner();
//!
//! let mut reader = XtcReader::new(Cursor::new(bytes));
//! let back = reader.next_frame().unwrap().unwrap();
//! assert_eq!(back.step, 10);
//! assert!((back.positions[1].x - 1.2).abs() < 1e-3);
//! ```

use std::io::{Read, Write};

use tpt_chem_core::vec3::Vec3;

use crate::error::{io_err, IoError, Result};

const FORMAT: &str = "XTC";
const MAGIC: i32 = 1995;
/// Frames with at most this many atoms store uncompressed floats.
const UNCOMPRESSED_MAX_ATOMS: usize = 9;
/// Index of the smallest delta width in [`MAGICINTS`].
const FIRSTIDX: usize = 9;
/// Index of the largest delta width in [`MAGICINTS`].
const LASTIDX: usize = MAGICINTS.len() - 1;
/// Largest per-component integer range the stream can code (24 bits).
const MAX_RANGE: u64 = 0x00ff_ffff;
/// Deltas per run, bounded by the 5-bit run length (3 components each).
const MAX_RUN_DELTAS: usize = 10;
/// Default `precision`: coordinates are kept to 0.001 nm.
const DEFAULT_PRECISION: f32 = 1000.0;
const NM_TO_ANG: f64 = 10.0;

/// Delta widths used by the coordinate coder, indexed by `smallidx`.
#[rustfmt::skip]
const MAGICINTS: [u32; 73] = [
    0, 0, 0, 0, 0, 0, 0, 0, 0,
    8, 10, 12, 16, 20, 25, 32, 40, 50, 64, 80, 101, 128, 161, 203, 256,
    322, 406, 512, 645, 812, 1024, 1290, 1625, 2048, 2580, 3250, 4096, 5160,
    6501, 8192, 10321, 13003, 16384, 20642, 26007, 32768, 41285, 52015,
    65536, 82570, 104031, 131072, 165140, 208063, 262144, 330280, 416127,
    524288, 660561, 832255, 1048576, 1321122, 1664510, 2097152, 2642245,
    3329021, 4194304, 5284491, 6658042, 8388608, 10568983, 13316085, 16777216,
];

/// One XTC frame. Positions and cell are in Ångström, the time in ps.
#[derive(Clone, Debug, PartialEq)]
pub struct XtcFrame {
    /// Simulation step.
    pub step: i32,
    /// Simulation time (ps).
    pub time: f32,
    /// Box vectors as rows (Ångström). Zero rows mean no periodic cell.
    pub cell: [[f64; 3]; 3],
    /// Atom positions (Ångström).
    pub positions: Vec<Vec3>,
}

/// Bits used for a mixed-radix triple with ranges `sizes`: the bit length of
/// the product. Run deltas use `smallidx` as their width instead, as in the
/// reference decoder (`small_bits` below).
fn sizeofints(sizes: [u32; 3]) -> u32 {
    let product: u128 = sizes.iter().map(|&s| u128::from(s)).product();
    128 - product.leading_zeros()
}

/// MSB-first bit reader over the compressed payload.
struct BitReader<'a> {
    data: &'a [u8],
    pos: usize,
}

impl BitReader<'_> {
    fn bits(&mut self, n: u32) -> Result<u32> {
        let mut v = 0u32;
        for _ in 0..n {
            let byte = *self
                .data
                .get(self.pos >> 3)
                .ok_or(IoError::UnexpectedEof { format: FORMAT })?;
            v = (v << 1) | u32::from((byte >> (7 - (self.pos & 7))) & 1);
            self.pos += 1;
        }
        Ok(v)
    }

    /// Reads a mixed-radix triple coded in `bits` bits (least significant
    /// byte first, then the top partial byte).
    fn ints(&mut self, bits: u32, sizes: [u32; 3]) -> Result<[u32; 3]> {
        let mut v: u128 = 0;
        let mut remaining = bits;
        let mut shift = 0;
        while remaining > 8 {
            v |= u128::from(self.bits(8)?) << shift;
            shift += 8;
            remaining -= 8;
        }
        if remaining > 0 {
            v |= u128::from(self.bits(remaining)?) << shift;
        }
        let s = sizes.map(u128::from);
        let n2 = v % s[2];
        v /= s[2];
        let n1 = v % s[1];
        v /= s[1];
        Ok([v as u32, n1 as u32, n2 as u32])
    }
}

/// MSB-first bit writer, the counterpart of [`BitReader`].
#[derive(Default)]
struct BitWriter {
    bytes: Vec<u8>,
    nbits: usize,
}

impl BitWriter {
    fn bits(&mut self, n: u32, v: u32) {
        for k in (0..n).rev() {
            if self.nbits & 7 == 0 {
                self.bytes.push(0);
            }
            let bit = ((v >> k) & 1) as u8;
            let last = self.bytes.last_mut().expect("byte pushed above");
            *last |= bit << (7 - (self.nbits & 7));
            self.nbits += 1;
        }
    }

    /// Writes a mixed-radix triple `n0 < s0`, `n1 < s1`, `n2 < s2` in `bits`
    /// bits, the inverse of [`BitReader::ints`].
    fn ints(&mut self, bits: u32, sizes: [u32; 3], nums: [u32; 3]) {
        let s = sizes.map(u128::from);
        let mut v = u128::from(nums[0]);
        v = v * s[1] + u128::from(nums[1]);
        v = v * s[2] + u128::from(nums[2]);
        let mut remaining = bits;
        let mut shift = 0;
        while remaining > 8 {
            self.bits(8, ((v >> shift) & 0xff) as u32);
            shift += 8;
            remaining -= 8;
        }
        if remaining > 0 {
            self.bits(
                remaining,
                ((v >> shift) & ((1u128 << remaining) - 1)) as u32,
            );
        }
    }
}

/// Decodes the compressed integer coordinates of `n_atoms` atoms.
fn decompress(
    data: &[u8],
    n_atoms: usize,
    minint: [i32; 3],
    maxint: [i32; 3],
    smallidx0: u32,
) -> Result<Vec<[i32; 3]>> {
    let sizeint = range_sizes(minint, maxint)?;
    let bitsize = sizeofints(sizeint);
    let mut smallidx = smallidx0 as usize;
    if !(FIRSTIDX..=LASTIDX).contains(&smallidx) {
        return Err(parse_err("delta width index out of range"));
    }
    let mut r = BitReader { data, pos: 0 };
    let mut out: Vec<[i32; 3]> = Vec::with_capacity(n_atoms);
    while out.len() < n_atoms {
        let abs = r.ints(bitsize, sizeint)?;
        let anchor = offset(abs, minint);

        let flag = r.bits(1)?;
        let mut is_smaller: i32 = 0;
        let mut run = 0usize;
        if flag == 1 {
            let run5 = r.bits(5)? as usize;
            is_smaller = (run5 % 3) as i32;
            run = run5 - is_smaller as usize;
            is_smaller -= 1;
        }

        if run > 0 {
            let width = MAGICINTS[smallidx];
            let smallnum = i64::from(width / 2);
            let widths = [width; 3];
            let small_bits = smallidx as u32;
            let mut prev = anchor;
            for k in (0..run).step_by(3) {
                let raw = r.ints(small_bits, widths)?;
                let cur = [0, 1, 2].map(|c| i64::from(raw[c]) + i64::from(prev[c]) - smallnum);
                let cur = cur.map(|v| v as i32);
                if k == 0 {
                    // The first delta is relative to the anchor, and it is
                    // written before the anchor.
                    out.push(cur);
                    out.push(anchor);
                } else {
                    out.push(cur);
                }
                prev = cur;
            }
            if out.len() > n_atoms {
                return Err(parse_err("run extends past the atom count"));
            }
        } else {
            out.push(anchor);
        }

        if is_smaller < 0 {
            smallidx = smallidx
                .checked_sub(1)
                .filter(|&i| i >= FIRSTIDX)
                .ok_or_else(|| parse_err("delta width underflow"))?;
        } else if is_smaller > 0 {
            smallidx += 1;
            if smallidx > LASTIDX {
                return Err(parse_err("delta width overflow"));
            }
        }
    }
    Ok(out)
}

/// Per-component integer range `max − min + 1`, checked against the 24-bit
/// limit of the coder.
fn range_sizes(minint: [i32; 3], maxint: [i32; 3]) -> Result<[u32; 3]> {
    let mut sizes = [0u32; 3];
    for c in 0..3 {
        let size = i64::from(maxint[c]) - i64::from(minint[c]) + 1;
        if size < 1 || size as u64 > MAX_RANGE {
            return Err(parse_err("coordinate range exceeds the 24-bit XTC limit"));
        }
        sizes[c] = size as u32;
    }
    Ok(sizes)
}

fn offset(v: [u32; 3], base: [i32; 3]) -> [i32; 3] {
    [0, 1, 2].map(|c| (i64::from(v[c]) + i64::from(base[c])) as i32)
}

fn parse_err(message: &str) -> IoError {
    IoError::Parse {
        format: FORMAT,
        line: None,
        message: message.to_string(),
    }
}

/// Output of [`compress`]: the header fields and the bit-packed payload.
struct Compressed {
    minint: [i32; 3],
    maxint: [i32; 3],
    smallidx: u32,
    payload: Vec<u8>,
}

/// Encodes integer coordinates.
fn compress(q: &[[i32; 3]]) -> Result<Compressed> {
    let n = q.len();
    let mut minint = [i32::MAX; 3];
    let mut maxint = [i32::MIN; 3];
    for p in q {
        for c in 0..3 {
            minint[c] = minint[c].min(p[c]);
            maxint[c] = maxint[c].max(p[c]);
        }
    }
    if n == 0 {
        minint = [0; 3];
        maxint = [0; 3];
    }
    let sizeint = range_sizes(minint, maxint)?;
    let bitsize = sizeofints(sizeint);

    // Choose the narrowest delta width that covers every delta the encoder
    // may use: consecutive atoms, and atoms two apart (the run chain).
    let mut worst: i64 = 0;
    for j in 0..n {
        for gap in [1, 2] {
            if j + gap < n {
                for (a, b) in q[j + gap].iter().zip(&q[j]) {
                    worst = worst.max((i64::from(*a) - i64::from(*b)).abs());
                }
            }
        }
    }
    let mut smallidx = FIRSTIDX;
    while smallidx < LASTIDX && i64::from(MAGICINTS[smallidx] / 2) <= worst {
        smallidx += 1;
    }
    let width = MAGICINTS[smallidx];
    let smallnum = i64::from(width / 2);
    let widths = [width; 3];
    let small_bits = smallidx as u32;
    let fits = |d: [i64; 3]| {
        d.iter()
            .all(|&v| v >= -smallnum && v < i64::from(width) - smallnum)
    };
    let sub = |a: [i32; 3], b: [i32; 3]| [0, 1, 2].map(|c| i64::from(a[c]) - i64::from(b[c]));

    let mut w = BitWriter::default();
    let mut i = 0;
    while i < n {
        if i + 1 < n {
            let anchor = q[i + 1];
            let first = sub(q[i], anchor);
            if fits(first) {
                let mut deltas = vec![first];
                let mut cur = q[i];
                let mut j = i + 2;
                while j < n && deltas.len() < MAX_RUN_DELTAS {
                    let d = sub(q[j], cur);
                    if !fits(d) {
                        break;
                    }
                    deltas.push(d);
                    cur = q[j];
                    j += 1;
                }
                // Anchor, run flag, run length (`is_smaller` = 0 is coded as 1).
                w.ints(bitsize, sizeint, sub_u(anchor, minint));
                w.bits(1, 1);
                w.bits(5, (3 * deltas.len() + 1) as u32);
                for d in &deltas {
                    let raw = [0, 1, 2].map(|c| (d[c] + smallnum) as u32);
                    w.ints(small_bits, widths, raw);
                }
                i += 1 + deltas.len();
                continue;
            }
        }
        w.ints(bitsize, sizeint, sub_u(q[i], minint));
        w.bits(1, 0);
        i += 1;
    }
    Ok(Compressed {
        minint,
        maxint,
        smallidx: smallidx as u32,
        payload: w.bytes,
    })
}

fn sub_u(a: [i32; 3], b: [i32; 3]) -> [u32; 3] {
    [0, 1, 2].map(|c| (i64::from(a[c]) - i64::from(b[c])) as u32)
}

/// Reads XTC frames from a byte stream.
pub struct XtcReader<R: Read> {
    inner: R,
}

impl<R: Read> XtcReader<R> {
    /// Wraps `inner`. No bytes are read until [`XtcReader::next_frame`].
    pub fn new(inner: R) -> Self {
        XtcReader { inner }
    }

    /// Reads the next frame, or `None` at the end of the stream.
    ///
    /// # Errors
    /// [`IoError::Parse`] for a bad magic number or a corrupt payload, and
    /// [`IoError::UnexpectedEof`] for a truncated frame.
    pub fn next_frame(&mut self) -> Result<Option<XtcFrame>> {
        let mut first = [0u8; 4];
        let got = self.inner.read(&mut first).map_err(io_err)?;
        if got == 0 {
            return Ok(None);
        }
        self.read_exact_or_eof(&mut first[got..])?;
        let magic = i32::from_be_bytes(first);
        if magic != MAGIC {
            return Err(parse_err("bad XTC magic number"));
        }
        let n_atoms = self.read_i32()?;
        let n_atoms = usize::try_from(n_atoms).map_err(|_| parse_err("negative atom count"))?;
        let step = self.read_i32()?;
        let time = self.read_f32()?;
        let mut cell = [[0.0f64; 3]; 3];
        for row in cell.iter_mut() {
            for v in row.iter_mut() {
                *v = f64::from(self.read_f32()?) * NM_TO_ANG;
            }
        }
        let n_again = self.read_i32()?;
        if usize::try_from(n_again).ok() != Some(n_atoms) {
            return Err(parse_err("atom count repeated inconsistently"));
        }

        let positions = if n_atoms <= UNCOMPRESSED_MAX_ATOMS {
            let mut out = Vec::with_capacity(n_atoms);
            for _ in 0..n_atoms {
                let v = [self.read_f32()?, self.read_f32()?, self.read_f32()?];
                out.push(nm_to_angstrom(v));
            }
            out
        } else {
            let precision = self.read_f32()?;
            if precision <= 0.0 {
                return Err(parse_err("non-positive precision"));
            }
            let minint = [self.read_i32()?, self.read_i32()?, self.read_i32()?];
            let maxint = [self.read_i32()?, self.read_i32()?, self.read_i32()?];
            let smallidx = self.read_i32()?;
            let nbytes = self.read_i32()?;
            let nbytes = usize::try_from(nbytes).map_err(|_| parse_err("negative byte count"))?;
            let mut payload = vec![0u8; nbytes];
            self.read_exact_or_eof(&mut payload)?;
            let pad = (4 - nbytes % 4) % 4;
            let mut skip = vec![0u8; pad];
            self.read_exact_or_eof(&mut skip)?;

            let smallidx = u32::try_from(smallidx).map_err(|_| parse_err("negative smallidx"))?;
            let q = decompress(&payload, n_atoms, minint, maxint, smallidx)?;
            q.into_iter()
                .map(|p| {
                    let nm = [0, 1, 2].map(|c| (p[c] as f32) / precision);
                    nm_to_angstrom(nm)
                })
                .collect()
        };
        Ok(Some(XtcFrame {
            step,
            time,
            cell,
            positions,
        }))
    }

    /// Reads every remaining frame.
    ///
    /// # Errors
    /// See [`XtcReader::next_frame`].
    pub fn collect(mut self) -> Result<Vec<XtcFrame>> {
        let mut frames = Vec::new();
        while let Some(f) = self.next_frame()? {
            frames.push(f);
        }
        Ok(frames)
    }

    fn read_exact_or_eof(&mut self, buf: &mut [u8]) -> Result<()> {
        self.inner.read_exact(buf).map_err(|e| {
            if e.kind() == std::io::ErrorKind::UnexpectedEof {
                IoError::UnexpectedEof { format: FORMAT }
            } else {
                io_err(e)
            }
        })
    }

    fn read_i32(&mut self) -> Result<i32> {
        let mut b = [0u8; 4];
        self.read_exact_or_eof(&mut b)?;
        Ok(i32::from_be_bytes(b))
    }

    fn read_f32(&mut self) -> Result<f32> {
        let mut b = [0u8; 4];
        self.read_exact_or_eof(&mut b)?;
        Ok(f32::from_be_bytes(b))
    }
}

fn nm_to_angstrom(v: [f32; 3]) -> Vec3 {
    Vec3::new(
        f64::from(v[0]) * NM_TO_ANG,
        f64::from(v[1]) * NM_TO_ANG,
        f64::from(v[2]) * NM_TO_ANG,
    )
}

/// Writes XTC frames to a byte stream.
pub struct XtcWriter<W: Write> {
    inner: W,
    precision: f32,
}

impl<W: Write> XtcWriter<W> {
    /// Writes with the default precision of 0.001 nm.
    pub fn new(inner: W) -> Self {
        XtcWriter {
            inner,
            precision: DEFAULT_PRECISION,
        }
    }

    /// Writes with `precision` integer units per nm (for example 1000 for
    /// 0.001 nm). Coordinates are rounded to the nearest unit.
    ///
    /// # Errors
    /// [`IoError::Parse`] if `precision` is not positive and finite.
    pub fn with_precision(inner: W, precision: f32) -> Result<Self> {
        if !(precision > 0.0 && precision.is_finite()) {
            return Err(parse_err("precision must be positive and finite"));
        }
        Ok(XtcWriter { inner, precision })
    }

    /// Writes one frame.
    ///
    /// # Errors
    /// [`IoError::Parse`] if the coordinate range does not fit the coder, or
    /// [`IoError::Io`] on a write failure.
    pub fn write_frame(&mut self, frame: &XtcFrame) -> Result<()> {
        let n = frame.positions.len();
        let n_i32 = i32::try_from(n).map_err(|_| parse_err("too many atoms for XTC"))?;
        let mut buf: Vec<u8> = Vec::new();
        buf.extend_from_slice(&MAGIC.to_be_bytes());
        buf.extend_from_slice(&n_i32.to_be_bytes());
        buf.extend_from_slice(&frame.step.to_be_bytes());
        buf.extend_from_slice(&frame.time.to_be_bytes());
        for row in &frame.cell {
            for &v in row {
                buf.extend_from_slice(&((v / NM_TO_ANG) as f32).to_be_bytes());
            }
        }
        buf.extend_from_slice(&n_i32.to_be_bytes());

        let nm: Vec<[f32; 3]> = frame
            .positions
            .iter()
            .map(|p| [p.x, p.y, p.z].map(|v| (v / NM_TO_ANG) as f32))
            .collect();
        if n <= UNCOMPRESSED_MAX_ATOMS {
            for p in &nm {
                for &v in p {
                    buf.extend_from_slice(&v.to_be_bytes());
                }
            }
        } else {
            let mut q = Vec::with_capacity(n);
            for p in &nm {
                let scaled = [0, 1, 2].map(|c| (p[c] * self.precision).round());
                if scaled.iter().any(|v| v.abs() >= i32::MAX as f32) {
                    return Err(parse_err("coordinate too large for XTC integers"));
                }
                q.push(scaled.map(|v| v as i32));
            }
            let c = compress(&q)?;
            buf.extend_from_slice(&self.precision.to_be_bytes());
            for v in c.minint.iter().chain(&c.maxint) {
                buf.extend_from_slice(&v.to_be_bytes());
            }
            buf.extend_from_slice(&(c.smallidx as i32).to_be_bytes());
            buf.extend_from_slice(&(c.payload.len() as i32).to_be_bytes());
            buf.extend_from_slice(&c.payload);
            buf.resize(buf.len() + (4 - c.payload.len() % 4) % 4, 0);
        }
        self.inner.write_all(&buf).map_err(io_err)
    }

    /// Returns the underlying writer.
    pub fn into_inner(self) -> W {
        self.inner
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn magic_table_has_expected_shape() {
        assert_eq!(MAGICINTS.len(), 73);
        assert_eq!(MAGICINTS[FIRSTIDX], 8);
        assert_eq!(MAGICINTS[LASTIDX], 16_777_216);
        assert!(MAGICINTS[FIRSTIDX..].windows(2).all(|w| w[0] < w[1]));
    }

    #[test]
    fn mixed_radix_round_trip() {
        let sizes = [1000u32, 77, 65_536];
        let bits = sizeofints(sizes);
        let mut w = BitWriter::default();
        let nums = [999u32, 76, 40_000];
        w.ints(bits, sizes, nums);
        let mut r = BitReader {
            data: &w.bytes,
            pos: 0,
        };
        assert_eq!(r.ints(bits, sizes).unwrap(), nums);
    }

    #[test]
    fn small_frame_is_stored_uncompressed() {
        let frame = XtcFrame {
            step: 3,
            time: 1.5,
            cell: [[0.0; 3]; 3],
            positions: vec![Vec3::new(0.25, -1.5, 2.0)],
        };
        let mut w = XtcWriter::new(Cursor::new(Vec::new()));
        w.write_frame(&frame).unwrap();
        let bytes = w.into_inner().into_inner();
        // magic, natoms, step, time, 9 box floats, natoms, 3 coordinates.
        assert_eq!(bytes.len(), 4 * (4 + 9 + 1 + 3));
        let back = XtcReader::new(Cursor::new(bytes)).collect().unwrap();
        assert_eq!(back.len(), 1);
        assert!((back[0].positions[0].x - 0.25).abs() < 1e-6);
    }

    #[test]
    fn compressed_round_trip_within_precision() {
        // Water-like triplets so runs and the adaptive width are used.
        let mut positions = Vec::new();
        for m in 0..40 {
            let c = f64::from(m) * 0.31;
            positions.push(Vec3::new(c, 0.1 * c, 0.2));
            positions.push(Vec3::new(c + 0.96, 0.1 * c, 0.2));
            positions.push(Vec3::new(c - 0.24, 0.1 * c + 0.93, 0.2));
        }
        let frame = XtcFrame {
            step: 42,
            time: 0.084,
            cell: [[30.0, 0.0, 0.0], [0.0, 30.0, 0.0], [0.0, 0.0, 30.0]],
            positions: positions.clone(),
        };
        let mut w = XtcWriter::new(Cursor::new(Vec::new()));
        w.write_frame(&frame).unwrap();
        let bytes = w.into_inner().into_inner();
        let back = XtcReader::new(Cursor::new(bytes)).collect().unwrap();
        assert_eq!(back[0].step, 42);
        for (a, b) in back[0].positions.iter().zip(&positions) {
            // Half a unit of 0.001 nm is 0.005 Å.
            assert!((a.x - b.x).abs() <= 0.0051, "{} vs {}", a.x, b.x);
            assert!((a.y - b.y).abs() <= 0.0051);
            assert!((a.z - b.z).abs() <= 0.0051);
        }
    }

    #[test]
    fn bad_magic_is_rejected() {
        let bytes = vec![0u8; 64];
        let err = XtcReader::new(Cursor::new(bytes)).next_frame().unwrap_err();
        assert!(matches!(err, IoError::Parse { .. }));
    }

    #[test]
    fn truncated_frame_is_reported() {
        let frame = XtcFrame {
            step: 0,
            time: 0.0,
            cell: [[0.0; 3]; 3],
            positions: (0..20).map(|i| Vec3::new(f64::from(i), 0.0, 0.0)).collect(),
        };
        let mut w = XtcWriter::new(Cursor::new(Vec::new()));
        w.write_frame(&frame).unwrap();
        let mut bytes = w.into_inner().into_inner();
        bytes.truncate(bytes.len() - 8);
        let err = XtcReader::new(Cursor::new(bytes)).next_frame().unwrap_err();
        assert!(matches!(err, IoError::UnexpectedEof { .. }));
    }

    #[test]
    fn empty_stream_yields_no_frames() {
        let frames = XtcReader::new(Cursor::new(Vec::new())).collect().unwrap();
        assert!(frames.is_empty());
    }
}
