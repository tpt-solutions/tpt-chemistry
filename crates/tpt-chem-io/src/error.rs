//! Error types shared by the format parsers.

use core::fmt;

/// Convenience alias for IO results.
pub type Result<T> = core::result::Result<T, IoError>;

/// Errors from parsing or writing chemical file formats.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum IoError {
    /// The input ended before a complete record was read.
    UnexpectedEof {
        /// Which format was being parsed.
        format: &'static str,
    },
    /// A record could not be parsed.
    Parse {
        /// Which format was being parsed.
        format: &'static str,
        /// 1-based line number, when known.
        line: Option<usize>,
        /// What went wrong.
        message: String,
    },
    /// A section or record required by the format is absent.
    MissingSection {
        /// Which format was being parsed.
        format: &'static str,
        /// The missing section keyword.
        section: String,
    },
    /// A number field could not be parsed.
    BadNumber {
        /// Which format was being parsed.
        format: &'static str,
        /// 1-based line number, when known.
        line: Option<usize>,
        /// The offending token.
        token: String,
    },
    /// An element symbol is unknown.
    UnknownElement {
        /// Which format was being parsed.
        format: &'static str,
        /// The offending symbol.
        symbol: String,
        /// 1-based line number, when known.
        line: Option<usize>,
    },
    /// Wraps a `std::io` failure from file operations.
    Io(String),
}

impl fmt::Display for IoError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            IoError::UnexpectedEof { format } => {
                write!(f, "{format}: unexpected end of input")
            }
            IoError::Parse {
                format,
                line,
                message,
            } => match line {
                Some(n) => write!(f, "{format} line {n}: {message}"),
                None => write!(f, "{format}: {message}"),
            },
            IoError::MissingSection { format, section } => {
                write!(f, "{format}: missing section {section:?}")
            }
            IoError::BadNumber {
                format,
                line,
                token,
            } => match line {
                Some(n) => write!(f, "{format} line {n}: bad number {token:?}"),
                None => write!(f, "{format}: bad number {token:?}"),
            },
            IoError::UnknownElement {
                format,
                symbol,
                line,
            } => match line {
                Some(n) => write!(f, "{format} line {n}: unknown element {symbol:?}"),
                None => write!(f, "{format}: unknown element {symbol:?}"),
            },
            IoError::Io(msg) => write!(f, "io error: {msg}"),
        }
    }
}

impl std::error::Error for IoError {}

/// Build a `Parse` error with the given line.
pub(crate) fn parse_err(format: &'static str, line: usize, message: impl Into<String>) -> IoError {
    IoError::Parse {
        format,
        line: Some(line),
        message: message.into(),
    }
}

/// Build a `BadNumber` error with the given line and token.
pub(crate) fn bad_number(format: &'static str, line: usize, token: &str) -> IoError {
    IoError::BadNumber {
        format,
        line: Some(line),
        token: token.to_string(),
    }
}

/// Wrap a `std::io::Error`.
pub(crate) fn io_err(e: std::io::Error) -> IoError {
    IoError::Io(e.to_string())
}

/// Parse an `f64` field, mapping failure to [`IoError::BadNumber`].
pub(crate) fn parse_f64(format: &'static str, line: usize, token: &str) -> Result<f64> {
    token
        .trim()
        .parse::<f64>()
        .map_err(|_| bad_number(format, line, token))
}
