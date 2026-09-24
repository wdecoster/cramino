use clap::Parser;
use log::info;
use metrics_processor::OutputFormat; // Import the enum

pub mod calculations;
pub mod extract_from_bam;
pub mod feather;
pub mod file_info;
pub mod histograms;
pub mod metrics;
pub mod metrics_processor;
pub mod phased;
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
    #[clap(short, long, value_parser = parse_threads, default_value_t = 4)]
    threads: usize,

    /// reference for decompressing cram
    #[clap(long, value_parser)]
    reference: Option<String>,

    /// Minimal aligned length of a read to be considered (full length for unmapped reads)
    #[clap(short, long, value_parser, default_value_t = 0)]
    min_read_len: usize,

    /// If histograms have to be generated (optionally specify output file as --hist=FILE)
    #[clap(long, value_parser, value_name = "FILE", num_args = 0..=1, require_equals = true, conflicts_with = "hist_count")]
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

    /// Also include unmapped reads, and estimate the identity from base qualities (disables --karyotype, --phased and --spliced)
    #[clap(long, value_parser)]
    ubam: bool,

    /// Output format (text, json, or tsv)
    #[clap(long, value_parser, default_value_t = OutputFormat::Text)]
    format: OutputFormat,

    /// Scale histogram bins by total basepairs in each bin (not just read count)
    #[clap(long, value_parser)]
    pub scaled: bool,

    /// Output read length histogram bin counts in TSV format (optionally specify output file as --hist-count=FILE), cannot be combined with --hist
    #[clap(long, value_parser, value_name = "FILE", num_args = 0..=1, require_equals = true, conflicts_with = "hist")]
    pub hist_count: Option<Option<String>>,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Rust ignores SIGPIPE, so that writing to a closed pipe (e.g. `cramino file.bam | head`)
    // panics. Restoring the default makes cramino exit quietly instead, like other tools.
    // SAFETY: called before any threads are started
    #[cfg(unix)]
    unsafe {
        libc::signal(libc::SIGPIPE, libc::SIG_DFL);
    }
    env_logger::init();
    let mut args = Cli::parse();
    utils::is_file(&args.input).unwrap_or_else(|err| utils::exit_with_error(&err));
    check_stdin_input(&args.input);
    // tsv output is a single table, which histograms on stdout would break
    if matches!(args.format, OutputFormat::Tsv)
        && (matches!(args.hist, Some(None)) || matches!(args.hist_count, Some(None)))
    {
        use clap::CommandFactory;
        Cli::command()
            .error(
                clap::error::ErrorKind::ArgumentConflict,
                "--hist and --hist-count require a file with --format tsv, e.g. --hist=histograms.txt",
            )
            .exit();
    }
    if utils::is_remote(&args.input) {
        // SAFETY: nothing has spawned threads yet, so the environment can be modified
        unsafe { extract_from_bam::setup_ssl_certificates() };
    }
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

fn parse_threads(value: &str) -> Result<usize, String> {
    match value.parse::<usize>() {
        Ok(0) => Err("the number of threads has to be at least 1".to_string()),
        Ok(threads) => Ok(threads),
        Err(_) => Err(format!("'{value}' is not a positive number")),
    }
}

fn check_stdin_input(input: &str) {
    if input == "-" {
        eprintln!(
            "Reading from stdin. If this is unexpected, make sure your input file is correctly specified."
        );
        // Check if stdin is connected to a terminal (interactive) using std library
        if std::io::IsTerminal::is_terminal(&std::io::stdin()) {
            eprintln!(
                "Warning: stdin appears to be a terminal. Did you mean to specify an input file?"
            );
        }
    }
}

#[cfg(test)]
#[ctor::ctor(unsafe)]
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
        scaled: false,
        hist_count: None,
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
        scaled: false,
        hist_count: None,
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
        scaled: false,
        hist_count: None,
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
        scaled: false,
        hist_count: None,
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
        scaled: false,
        hist_count: None,
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
        scaled: false,
        hist_count: None,
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
        scaled: false,
        hist_count: None,
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
        scaled: false,
        hist_count: None,
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
        scaled: false,
        hist_count: None,
    };

    let (metrics, header) = extract_from_bam::extract(&args);
    assert!(metrics.lengths.as_ref().unwrap().is_empty());
    assert!(metrics_processor::process_metrics(metrics, &args, header).is_ok());
}

#[test]
fn extract_hist_scaled() {
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
        scaled: true, // Set scaled to true for this test
        hist_count: None,
    };
    let (metrics, header) = extract_from_bam::extract(&args);
    assert!(metrics_processor::process_metrics(metrics, &args, header).is_ok())
}

#[test]
fn extract_hist_count() {
    let args = Cli {
        input: "test-data/small-test-phased.bam".to_string(),
        threads: 8,
        reference: None,
        min_read_len: 0,
        hist: None,
        arrow: None,
        karyotype: false,
        phased: false,
        spliced: false,
        ubam: false,
        format: OutputFormat::Text,
        scaled: false,
        hist_count: Some(None),
    };
    let (metrics, header) = extract_from_bam::extract(&args);
    assert!(metrics_processor::process_metrics(metrics, &args, header).is_ok())
}

#[test]
fn extract_hist_count_with_high_min_length() {
    // Test that --hist-count works with empty results
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
        scaled: false,
        hist_count: Some(None),
    };

    let (metrics, header) = extract_from_bam::extract(&args);
    assert!(metrics.lengths.as_ref().unwrap().is_empty());
    assert!(metrics_processor::process_metrics(metrics, &args, header).is_ok());
}

/// Pins which records enter the N50 calculations.
/// The N50 (and N75) is over the full length of each read, i.e. the primary alignments including
/// clipped bases, so that it is the same for a BAM and the ubam of the same reads.
/// The N50 aligned (like the yield, median length aligned and mean length aligned) is over the primary and supplementary
/// alignments without clipped bases. Secondary alignments are never counted, and unmapped reads
/// only with --ubam.
#[test]
fn n50_of_reads_and_of_alignments() {
    let bam_args = Cli {
        input: "test-data/small-test-phased.bam".to_string(),
        threads: 8,
        reference: None,
        min_read_len: 0,
        hist: None,
        arrow: None,
        karyotype: false,
        phased: false,
        spliced: false,
        ubam: false,
        format: OutputFormat::Text,
        scaled: false,
        hist_count: None,
    };
    let (bam, _) = extract_from_bam::extract(&bam_args);
    // 6176 primary alignments, the 1240 supplementary and 689 secondary alignments are not reads
    assert_eq!(bam.read_lengths.len(), 6176);
    assert_eq!(bam.all_reads, 6176);
    assert_eq!(bam.num_reads, 6176);
    let read_total = bam.read_lengths.iter().map(|l| *l as u128).sum::<u128>();
    assert_eq!(read_total, 76748714);
    assert_eq!(
        calculations::get_n(&bam.read_lengths, read_total, 0.50),
        21885
    );

    // 6176 primary + 1240 supplementary alignments
    let lengths = bam.lengths.unwrap();
    assert_eq!(lengths.len(), 7416);
    // the aligned bases exceed the read bases, as supplementary alignments can overlap
    let aligned_total = lengths.iter().map(|l| *l as u128).sum::<u128>();
    assert_eq!(aligned_total, 76888339);
    assert_eq!(calculations::get_n(&lengths, aligned_total, 0.50), 19294);

    // the ubam of the same reads has the same read length N50
    let ubam_args = Cli {
        input: "test-data/small-test-ubam.bam".to_string(),
        ubam: true,
        ..bam_args
    };
    let (ubam, _) = extract_from_bam::extract(&ubam_args);
    let read_total = ubam.read_lengths.iter().map(|l| *l as u128).sum::<u128>();
    assert_eq!(ubam.read_lengths.len(), 6176);
    assert_eq!(read_total, 76748714);
    assert_eq!(
        calculations::get_n(&ubam.read_lengths, read_total, 0.50),
        21885
    );
    assert_eq!(ubam.lengths.unwrap().len(), 6176);
}

/// The minimum length applies to the same aligned length that is reported, and is inclusive
#[test]
fn min_read_len_filters_on_reported_length() {
    let args = Cli {
        input: "test-data/small-test-phased.bam".to_string(),
        threads: 8,
        reference: None,
        min_read_len: 5000,
        hist: None,
        arrow: None,
        karyotype: false,
        phased: false,
        spliced: false,
        ubam: false,
        format: OutputFormat::Text,
        scaled: false,
        hist_count: None,
    };
    let (data, _) = extract_from_bam::extract(&args);
    let lengths = data.lengths.unwrap();
    assert!(!lengths.is_empty());
    assert!(lengths.iter().all(|l| *l >= 5000));
}
