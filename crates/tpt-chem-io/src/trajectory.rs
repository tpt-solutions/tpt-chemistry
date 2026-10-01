//! Streaming XYZ trajectory I/O for large MD outputs.
//!
//! Frames are read and written one at a time, so a multi-gigabyte
//! trajectory never needs to be held in memory — spec.txt §4's
//! "Trajectory streaming for large MD outputs".

use std::io::{BufReader, BufWriter, Read, Write};

use tpt_chem_core::molecule::Molecule;

use crate::error::Result;
use crate::xyz;

/// Incremental reader over a (possibly huge) multi-frame XYZ stream.
pub struct XyzTrajectoryReader<R: Read> {
    inner: BufReader<R>,
}

impl<R: Read> XyzTrajectoryReader<R> {
    /// Wrap any byte source.
    pub fn new(inner: R) -> Self {
        XyzTrajectoryReader {
            inner: BufReader::new(inner),
        }
    }

    /// Read the next frame, or `None` at the end of the stream.
    ///
    /// # Errors
    /// [`crate::IoError`] on malformed frames.
    pub fn next_frame(&mut self) -> Result<Option<Molecule>> {
        xyz::read_frame(&mut self.inner)
    }

    /// Collect every remaining frame (convenience for small files).
    ///
    /// # Errors
    /// [`crate::IoError`] on malformed frames.
    pub fn collect(mut self) -> Result<Vec<Molecule>> {
        let mut frames = Vec::new();
        while let Some(mol) = self.next_frame()? {
            frames.push(mol);
        }
        Ok(frames)
    }
}

/// Incremental writer producing a multi-frame XYZ stream.
pub struct XyzTrajectoryWriter<W: Write> {
    inner: BufWriter<W>,
    frames: usize,
}

impl<W: Write> XyzTrajectoryWriter<W> {
    /// Wrap any byte sink.
    pub fn new(inner: W) -> Self {
        XyzTrajectoryWriter {
            inner: BufWriter::new(inner),
            frames: 0,
        }
    }

    /// Append one frame with a comment line.
    ///
    /// # Errors
    /// Propagates `std::io` failures.
    pub fn write_frame(&mut self, mol: &Molecule, comment: &str) -> std::io::Result<()> {
        xyz::write_frame_to(&mut self.inner, mol, comment)?;
        self.frames += 1;
        Ok(())
    }

    /// Frames written so far.
    pub fn frames_written(&self) -> usize {
        self.frames
    }

    /// Flush and return the underlying writer.
    ///
    /// # Errors
    /// Propagates `std::io` failures from the final flush.
    pub fn finish(mut self) -> std::io::Result<W> {
        self.inner.flush()?;
        Ok(self.inner.into_inner()?)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    const WATER: &str = "3\nw\nO 0 0 0.117\nH 0 0.757 -0.469\nH 0 -0.757 -0.469\n";

    #[test]
    fn stream_roundtrip() {
        let mut source = String::new();
        source.push_str(WATER);
        source.push_str(WATER);

        let frames = XyzTrajectoryReader::new(Cursor::new(source.clone()))
            .collect()
            .unwrap();
        assert_eq!(frames.len(), 2);

        // Re-serialize through the writer and compare parse results.
        let mut buf = Vec::new();
        {
            let mut w = XyzTrajectoryWriter::new(&mut buf);
            for mol in &frames {
                w.write_frame(mol, "streamed").unwrap();
            }
            assert_eq!(w.frames_written(), 2);
            w.finish().unwrap();
        }
        let text = String::from_utf8(buf).unwrap();
        let again = XyzTrajectoryReader::new(Cursor::new(text))
            .collect()
            .unwrap();
        assert_eq!(again.len(), 2);
        for (a, b) in frames.iter().zip(again.iter()) {
            assert_eq!(a.len(), b.len());
            for i in 0..a.len() {
                assert_eq!(a.atom(i.into()).pos, b.atom(i.into()).pos);
            }
        }
    }

    #[test]
    fn empty_stream_is_none() {
        let mut r = XyzTrajectoryReader::new(Cursor::new(String::new()));
        assert!(r.next_frame().unwrap().is_none());
    }
}
