use clap::Parser;
use log::info;
use metrics_processor::OutputFormat;  // Import the enum

pub mod calculations;
pub mod extract_from_bam;
pub mod feather;
pub mod file_info;
pub mod histograms;
pub mod metrics;
pub mod metrics_processor;
pub mod phased;
pub mod splicing;
pub mod text_output;
pub mod tsv_output;
pub mod utils;

// The arguments end up in the Cli struct
#[derive(Parser, Debug)]
#[clap(author, version, about="Tool to extract QC metrics from cram or bam", long_about = None)]
pub struct Cli {
    /// cram or bam file to check
    #[clap(value_parser, default_value = "-")]
    input: String,

    /// Number of parallel decompression threads to use
    #[clap(short, long, value_parser, default_value_t = 4)]
    threads: usize,

    /// reference for decompressing cram
    #[clap(long, value_parser)]
    reference: Option<String>,

    /// Minimal length of read to be considered
    #[clap(short, long, value_parser, default_value_t = 0)]
    min_read_len: usize,

    /// If histograms have to be generated (optionally specify output file)
    #[clap(long, value_parser, value_name = "FILE", num_args = 0..=1)]
    hist: Option<Option<String>>,

    /// Write data to an arrow format file
    #[clap(long, value_parser)]
    arrow: Option<String>,

    /// Provide normalized number of reads per chromosome
    #[clap(long, value_parser)]
    karyotype: bool,

    /// Calculate metrics for phased reads
    #[clap(long, value_parser)]
    phased: bool,

    /// Provide metrics for spliced data
    #[clap(long, value_parser)]
    spliced: bool,

    /// Provide metrics for unaligned reads
    #[clap(long, value_parser)]
    ubam: bool,

    /// Output format (text, json, or tsv)
    #[clap(long, value_parser, default_value_t = OutputFormat::Text)]
    format: OutputFormat,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    env_logger::init();
    let mut args = Cli::parse();
    utils::is_file(&args.input).unwrap_or_else(|_| panic!("Path to input file {} is invalid", args.input));
    if args.ubam {
        args.karyotype = false;
        args.phased = false;
        args.spliced = false;
    };
    info!("Collected arguments");
    let (metrics, header) = extract_from_bam::extract(&args);
    info!("Extracted metrics");
    metrics_processor::process_metrics(metrics, &args, header)?;
    info!("Finished");
    Ok(())
}


#[cfg(test)]
#[ctor::ctor]
fn init() {
    env_logger::init();
}

#[test]
fn verify_app() {
    use clap::CommandFactory;
    Cli::command().debug_assert()
}

#[test]
fn extract() {
    let args = Cli {
        input: "test-data/small-test-phased.bam".to_string(),
        threads: 8,
        reference: None,
        min_read_len: 0,
        hist: Some(None),
        arrow: Some("test.feather".to_string()),
        karyotype: true,
        phased: true,
        spliced: false,
        ubam: false,
        format: OutputFormat::Text,
    };
    let (metrics, header) = extract_from_bam::extract(&args);
    assert!(metrics_processor::process_metrics(metrics, &args, header).is_ok())
}

// this test is ignored because it uses a local reference file
#[ignore]
#[test]
fn extract_cram() {
    let args = Cli {
        input: "test-data/small-test-phased.cram".to_string(),
        threads: 8,
        reference: Some("/home/wdecoster/reference/GRCh38.fa".to_string()),
        min_read_len: 0,
        hist: None,
        arrow: None,
        karyotype: false,
        phased: false,
        spliced: false,
        ubam: false,
        format: OutputFormat::Text,
    };
    let (metrics, header) = extract_from_bam::extract(&args);
    assert!(metrics_processor::process_metrics(metrics, &args, header).is_ok())
}

#[test]
fn extract_ubam() {
    let args = Cli {
        input: "test-data/small-test-ubam.bam".to_string(),
        threads: 8,
        reference: None,
        min_read_len: 0,
        hist: Some(None),
        arrow: Some("test.feather".to_string()),
        karyotype: false,
        phased: false,
        spliced: false,
        ubam: true,
        format: OutputFormat::Text,
    };
    let (metrics, header) = extract_from_bam::extract(&args);
    assert!(metrics_processor::process_metrics(metrics, &args, header).is_ok())
}

// this test is ignored because it uses a local reference file and takes a very long time
#[ignore]
#[test]
fn extract_url() {
    let args = Cli {
        input: "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1KG_ONT_VIENNA/hg38/HG00096.hg38.cram".to_string(),
        threads: 8,
        reference: Some("/home/wdecoster/local/1KG_ONT_VIENNA_hg38.fa.gz".to_string()),
        min_read_len: 0,
        hist: Some(None),
        arrow: None,
        karyotype: false,
        phased: false,
        spliced: false,
        ubam: false,
        format: OutputFormat::Text,
    };
    let (metrics, header) = extract_from_bam::extract(&args);
    assert!(metrics_processor::process_metrics(metrics, &args, header).is_ok())
}

#[test]
fn extract_json() {
    let args = Cli {
        input: "test-data/small-test-phased.bam".to_string(),
        threads: 8,
        reference: None,
        min_read_len: 0,
        hist: Some(None),
        arrow: None,
        karyotype: true,
        phased: true,
        spliced: false,
        ubam: false,
        format: OutputFormat::Json,
    };
    let (metrics, header) = extract_from_bam::extract(&args);
    assert!(metrics_processor::process_metrics(metrics, &args, header).is_ok())
}

#[test]
fn extract_tsv() {
    let args = Cli {
        input: "test-data/small-test-phased.bam".to_string(),
        threads: 8,
        reference: None,
        min_read_len: 0,
        hist: Some(Some("hist.txt".to_string())),
        arrow: None,
        karyotype: true,
        phased: true,
        spliced: false,
        ubam: false,
        format: OutputFormat::Tsv,
    };
    let (metrics, header) = extract_from_bam::extract(&args);
    assert!(metrics_processor::process_metrics(metrics, &args, header).is_ok())
}

#[test]
fn extract_with_high_min_length() {
    // Use a minimum read length higher than any read in the test file
    let args = Cli {
        input: "test-data/small-test-phased.bam".to_string(),
        threads: 8,
        reference: None,
        min_read_len: 1_000_000, // Set very high to ensure no reads match
        hist: None,
        arrow: None,
        karyotype: false,
        phased: false,
        spliced: false,
        ubam: false,
        format: OutputFormat::Text,
    };
    
    // The test should still run without panicking
    let (metrics, header) = extract_from_bam::extract(&args);
    assert!(metrics.lengths.as_ref().unwrap().is_empty());
    assert!(metrics_processor::process_metrics(metrics, &args, header).is_ok());
}

#[test]
fn extract_json_with_high_min_length() {
    let args = Cli {
        input: "test-data/small-test-phased.bam".to_string(),
        threads: 8,
        reference: None,
        min_read_len: 1_000_000, // Set very high to ensure no reads match
        hist: None,
        arrow: None,
        karyotype: false,
        phased: false,
        spliced: false,
        ubam: false,
        format: OutputFormat::Json,
    };
    
    let (metrics, header) = extract_from_bam::extract(&args);
    assert!(metrics.lengths.as_ref().unwrap().is_empty());
    assert!(metrics_processor::process_metrics(metrics, &args, header).is_ok());
}

#[test]
fn extract_tsv_with_high_min_length() {
    let args = Cli {
        input: "test-data/small-test-phased.bam".to_string(),
        threads: 8,
        reference: None,
        min_read_len: 1_000_000, // Set very high to ensure no reads match
        hist: None,
        arrow: None,
        karyotype: false,
        phased: false,
        spliced: false,
        ubam: false,
        format: OutputFormat::Tsv,
    };
    
    let (metrics, header) = extract_from_bam::extract(&args);
    assert!(metrics.lengths.as_ref().unwrap().is_empty());
    assert!(metrics_processor::process_metrics(metrics, &args, header).is_ok());
}