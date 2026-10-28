//! Run manifests: a small JSON record written alongside every output so a
//! result can be traced back to the code, inputs, seed and parameters that
//! produced it.
//!
//! The manifest has no timestamp on purpose. The same run then writes the
//! same bytes, which makes manifests diffable and lets tests compare them.
//! Add a timestamp in the caller if you need one.
//!
//! Output files are identified by their byte length and an FNV-1a 64-bit
//! hash. That is a fingerprint for detecting changes, not a cryptographic
//! digest.
//!
//! # Examples
//!
//! ```
//! use tpt_chemistry::manifest::RunManifest;
//!
//! let m = RunManifest::new("md-cli", "0.1.0")
//!     .seed(42)
//!     .param("temperature_k", 300.0)
//!     .param("steps", 1000);
//! let json = m.to_json();
//! assert!(json.contains("\"seed\": 42"));
//! assert!(json.contains("\"temperature_k\": \"300\""));
//! ```

use std::fmt::Display;
use std::fs;
use std::io;
use std::path::{Path, PathBuf};

/// Version of the umbrella crate that produced the run.
pub const UMBRELLA_VERSION: &str = env!("CARGO_PKG_VERSION");

/// Feature flags this build of the umbrella crate was compiled with.
///
/// Only features that gate a constituent crate are listed.
pub fn enabled_features() -> Vec<&'static str> {
    let mut out = Vec::new();
    if cfg!(feature = "core") {
        out.push("core");
    }
    if cfg!(feature = "io") {
        out.push("io");
    }
    if cfg!(feature = "md") {
        out.push("md");
    }
    if cfg!(feature = "quantum") {
        out.push("quantum");
    }
    if cfg!(feature = "kinetics") {
        out.push("kinetics");
    }
    if cfg!(feature = "aqueous") {
        out.push("aqueous");
    }
    if cfg!(feature = "crystal") {
        out.push("crystal");
    }
    if cfg!(feature = "verify") {
        out.push("verify");
    }
    out
}

/// FNV-1a 64-bit hash of `bytes`.
pub fn fnv1a64(bytes: &[u8]) -> u64 {
    const OFFSET: u64 = 0xcbf2_9ce4_8422_2325;
    const PRIME: u64 = 0x0000_0100_0000_01b3;
    bytes
        .iter()
        .fold(OFFSET, |h, &b| (h ^ u64::from(b)).wrapping_mul(PRIME))
}

/// A recorded output file: its path, byte length and FNV-1a hash.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OutputRecord {
    /// Path as given to [`RunManifest::output_file`].
    pub path: String,
    /// Length in bytes.
    pub bytes: u64,
    /// FNV-1a 64-bit hash of the contents.
    pub fnv1a64: u64,
}

/// Provenance record for one run.
#[derive(Debug, Clone, PartialEq)]
pub struct RunManifest {
    tool: String,
    tool_version: String,
    seed: Option<u64>,
    params: Vec<(String, String)>,
    inputs: Vec<String>,
    outputs: Vec<OutputRecord>,
}

impl RunManifest {
    /// Starts a manifest for `tool` at `tool_version`.
    pub fn new(tool: &str, tool_version: &str) -> Self {
        Self {
            tool: tool.to_string(),
            tool_version: tool_version.to_string(),
            seed: None,
            params: Vec::new(),
            inputs: Vec::new(),
            outputs: Vec::new(),
        }
    }

    /// Records the random seed used by the run.
    #[must_use]
    pub fn seed(mut self, seed: u64) -> Self {
        self.seed = Some(seed);
        self
    }

    /// Records a named parameter. The value is stored as its `Display` text.
    #[must_use]
    pub fn param(mut self, key: &str, value: impl Display) -> Self {
        self.params.push((key.to_string(), value.to_string()));
        self
    }

    /// Records an input path (as given, not resolved).
    #[must_use]
    pub fn input(mut self, path: &str) -> Self {
        self.inputs.push(path.to_string());
        self
    }

    /// Records an output file by reading it and hashing its contents.
    ///
    /// # Errors
    /// Returns any I/O error from reading `path`.
    pub fn output_file(mut self, path: &Path) -> io::Result<Self> {
        let bytes = fs::read(path)?;
        self.outputs.push(OutputRecord {
            path: path.display().to_string(),
            bytes: bytes.len() as u64,
            fnv1a64: fnv1a64(&bytes),
        });
        Ok(self)
    }

    /// Records an output that was produced in memory and has not been
    /// written yet, from its contents.
    #[must_use]
    pub fn output_bytes(mut self, path: &str, contents: &[u8]) -> Self {
        self.outputs.push(OutputRecord {
            path: path.to_string(),
            bytes: contents.len() as u64,
            fnv1a64: fnv1a64(contents),
        });
        self
    }

    /// Serializes to pretty-printed JSON with a fixed key order.
    ///
    /// Parameter values are written as strings so that integers and floats
    /// keep the exact text the caller gave them.
    pub fn to_json(&self) -> String {
        let mut s = String::from("{\n");
        s.push_str(&kv("tool", &quote(&self.tool)));
        s.push_str(&kv("tool_version", &quote(&self.tool_version)));
        s.push_str(&kv("umbrella_version", &quote(UMBRELLA_VERSION)));
        let features: Vec<String> = enabled_features().iter().map(|f| quote(f)).collect();
        s.push_str(&kv("features", &array(&features)));
        s.push_str(&kv(
            "seed",
            &self.seed.map_or("null".to_string(), |v| v.to_string()),
        ));
        let params: Vec<String> = self
            .params
            .iter()
            .map(|(k, v)| format!("{}: {}", quote(k), quote(v)))
            .collect();
        s.push_str(&kv("params", &object(&params)));
        let inputs: Vec<String> = self.inputs.iter().map(|p| quote(p)).collect();
        s.push_str(&kv("inputs", &array(&inputs)));
        let outputs: Vec<String> = self
            .outputs
            .iter()
            .map(|o| {
                format!(
                    "{{\"path\": {}, \"bytes\": {}, \"fnv1a64\": \"{:016x}\"}}",
                    quote(&o.path),
                    o.bytes,
                    o.fnv1a64
                )
            })
            .collect();
        // The last key carries no trailing comma.
        s.push_str(&format!("  \"outputs\": {}\n", array(&outputs)));
        s.push_str("}\n");
        s
    }

    /// Writes the manifest to `<output>.manifest.json` next to `output`.
    ///
    /// Returns the path written.
    ///
    /// # Errors
    /// Returns any I/O error from writing the file.
    pub fn write_alongside(&self, output: &Path) -> io::Result<PathBuf> {
        let mut name = output.as_os_str().to_owned();
        name.push(".manifest.json");
        let path = PathBuf::from(name);
        fs::write(&path, self.to_json())?;
        Ok(path)
    }
}

fn kv(key: &str, value: &str) -> String {
    format!("  {}: {},\n", quote(key), value)
}

fn array(items: &[String]) -> String {
    if items.is_empty() {
        "[]".to_string()
    } else {
        format!("[{}]", items.join(", "))
    }
}

fn object(items: &[String]) -> String {
    if items.is_empty() {
        "{}".to_string()
    } else {
        format!("{{{}}}", items.join(", "))
    }
}

/// JSON string literal with the required escapes.
fn quote(s: &str) -> String {
    let mut out = String::with_capacity(s.len() + 2);
    out.push('"');
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if (c as u32) < 0x20 => out.push_str(&format!("\\u{:04x}", c as u32)),
            c => out.push(c),
        }
    }
    out.push('"');
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fnv_matches_reference_vectors() {
        assert_eq!(fnv1a64(b""), 0xcbf2_9ce4_8422_2325);
        assert_eq!(fnv1a64(b"a"), 0xaf63_dc4c_8601_ec8c);
    }

    #[test]
    fn json_is_stable_and_escaped() {
        let m = RunManifest::new("t\"ool", "1.0")
            .seed(7)
            .param("note", "a\\b\nc")
            .input("in.xyz")
            .output_bytes("out.xyz", b"abc");
        let a = m.to_json();
        let b = m.to_json();
        assert_eq!(a, b, "output must be deterministic");
        assert!(a.contains("\"tool\": \"t\\\"ool\""));
        assert!(a.contains("\"note\": \"a\\\\b\\nc\""));
        assert!(a.contains("\"seed\": 7"));
        assert!(a.contains("\"inputs\": [\"in.xyz\"]"));
        assert!(a.contains(&format!("\"fnv1a64\": \"{:016x}\"", fnv1a64(b"abc"))));
    }

    #[test]
    fn missing_seed_and_empty_lists_are_valid_json_shapes() {
        let json = RunManifest::new("x", "0").to_json();
        assert!(json.contains("\"seed\": null"));
        assert!(json.contains("\"params\": {}"));
        assert!(json.contains("\"outputs\": []"));
        assert!(json.ends_with("}\n"));
    }

    #[test]
    fn write_alongside_appends_suffix() {
        let dir = std::env::temp_dir().join(format!("tpt-manifest-{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        let out = dir.join("run.xyz");
        fs::write(&out, "3\n").unwrap();
        let m = RunManifest::new("t", "0").output_file(&out).unwrap();
        let written = m.write_alongside(&out).unwrap();
        assert_eq!(written, dir.join("run.xyz.manifest.json"));
        let text = fs::read_to_string(&written).unwrap();
        assert!(text.contains("\"bytes\": 2"));
        fs::remove_dir_all(&dir).unwrap();
    }
}
