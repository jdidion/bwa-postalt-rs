extern crate structopt;

use std::fmt;
use std::path::PathBuf;
use std::str::FromStr;

use rust_htslib::bam::Format;
pub use structopt::StructOpt;

// wrap htslib Format type because it doesn't implement FromStr
#[derive(Debug)]
pub enum HtsFormat {
    Sam,
    Bam,
    Cram,
}

impl HtsFormat {
    pub fn from_string_ignore_case(s: &str) -> HtsFormat {
        match s.to_lowercase().as_str() {
            "sam" => HtsFormat::Sam,
            "bam" => HtsFormat::Bam,
            "cram" => HtsFormat::Cram,
            _ => panic!("invalid format {}", s),
        }
    }

    pub fn convert(fmt: HtsFormat) -> Format {
        match fmt {
            HtsFormat::Sam => Format::Sam,
            HtsFormat::Bam => Format::Bam,
            HtsFormat::Cram => Format::Cram,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ParseHtsFormatError {
    msg: String,
}

impl fmt::Display for ParseHtsFormatError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        "provided string was not one of {sam,bam,cram}".fmt(f)
    }
}

impl FromStr for HtsFormat {
    type Err = ParseHtsFormatError;

    fn from_str(s: &str) -> Result<Self, ParseHtsFormatError> {
        match s.to_lowercase().as_str() {
            "sam" => Ok(HtsFormat::Sam),
            "bam" => Ok(HtsFormat::Bam),
            "cram" => Ok(HtsFormat::Cram),
            other => Err(ParseHtsFormatError {
                msg: format!("invalid format {}", other),
            }),
        }
    }
}

#[derive(Debug, StructOpt)]
pub struct Opt {
    /// Prefix of output files containting sequences matching HLA genes
    #[structopt(short = "p")]
    pub prefix: Option<String>,

    /// Reduce mapQ to 0 if not overlapping lifted best and pa < RATIO
    #[structopt(name = "RATIO", short = "r", default_value = "1")]
    pub min_pa_ratio: f32,

    /// ALT-to-REF alignment SAM/BAM file
    #[structopt(parse(from_os_str))]
    pub alt: PathBuf,

    /// SAM/BAM file to process; stdin if not present
    #[structopt(parse(from_os_str))]
    pub sam: Option<PathBuf>,

    /// Output SAM/BAM file; stdout if not present
    #[structopt(parse(from_os_str))]
    pub output: Option<PathBuf>,

    /// Output file format; detected from output or sam
    /// filename if not present
    #[structopt(long = "format")]
    pub output_format: Option<HtsFormat>,
}
