use clap::Parser;
use extract_from_bam::Data;
use log::{error, info};
use std::path::PathBuf;

pub mod calculations;
pub mod extract_from_bam;
pub mod feather;
pub mod file_info;
pub mod histograms;
pub mod karyotype;
pub mod phased;
pub mod splicing;
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

    /// If histograms have to be generated
    #[clap(long, value_parser)]
    hist: bool,

    /// If a checksum has to be calculated [DEPRECATED]
    #[clap(long, value_parser)]
    checksum: bool,

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
}

pub fn is_file(pathname: &str) -> Result<(), String> {
    if pathname == "-"
        || pathname.starts_with("http")
        || pathname.starts_with("ftp")
        || pathname.starts_with("s3")
    {
        return Ok(());
    }
    let path = PathBuf::from(pathname);
    if path.is_file() {
        Ok(())
    } else {
        Err(format!("Input file {} is invalid", path.display()))
    }
}

fn main() -> Result<(), rust_htslib::errors::Error> {
    env_logger::init();
    let mut args = Cli::parse();
    is_file(&args.input).unwrap_or_else(|_| panic!("Path to input file {} is invalid", args.input));
    if args.ubam {
        args.karyotype = false;
        args.phased = false;
        args.spliced = false;
    };
    if args.checksum {
        eprintln!(
            "Calculating a checksum is deprecated, let me you if you disagree and want it back"
        );
    }
    info!("Collected arguments");
    let (metrics, header) = extract_from_bam::extract(&args);

    metrics_from_bam(metrics, args, header)?;
    info!("Finished");
    Ok(())
}

fn metrics_from_bam(
    metrics: Data,
    args: Cli,
    header: rust_htslib::bam::Header,
) -> Result<(), rust_htslib::errors::Error> {
    let bam = file_info::BamFile { path: args.input };
    println!("File name\t{}", bam.file_name());

    let genome_size = utils::get_genome_size(&header)?;
    generate_main_output(
        metrics.lengths.as_ref().unwrap(),
        metrics.identities.as_ref(),
        genome_size,
        metrics.all_counts,
    );

    println!("Path\t{}", bam);
    println!("Creation time\t{}", bam.file_time());

    let phaseblocks = if args.phased {
        Some(phased::phase_metrics(
            metrics.tids.as_ref().unwrap(),
            metrics.starts.unwrap(),
            metrics.ends.unwrap(),
            metrics.phasesets.as_ref().unwrap(),
        ))
    } else {
        None
    };
    if args.karyotype {
        karyotype::make_karyotype(metrics.tids.as_ref().unwrap(), header);
    }
    let exon_counts = if let Some(mut exon_counts) = metrics.exons {
        exon_counts.sort_unstable();
        splicing::splice_metrics(&exon_counts);
        exon_counts
    } else {
        Vec::new()
    };
    if args.hist {
        histograms::make_histogram_lengths(metrics.lengths.as_ref().unwrap());
        if !args.ubam {
            histograms::make_histogram_identities(metrics.identities.as_ref().unwrap());
        }
        if args.phased {
            histograms::make_histogram_phaseblocks(&phaseblocks.unwrap())
        }
        if args.spliced {
            histograms::make_histogram_exons(&exon_counts);
        }
    }
    Ok(())
}

fn generate_main_output(
    lengths: &Vec<u128>,
    identities: Option<&Vec<f64>>,
    genome_size: u64,
    all_reads: usize,
) {
    let num_reads = lengths.len();
    if num_reads < 2 {
        error!("Not enough reads to calculate metrics!");
        panic!();
    }
    let data_yield: u128 = lengths.iter().sum::<u128>();
    println!("Number of alignments\t{num_reads}");
    println!(
        "% from total reads\t{:.2}",
        (num_reads as f64) / (all_reads as f64) * 100.0
    );
    println!("Yield [Gb]\t{:.2}", data_yield as f64 / 1e9);
    println!(
        "Mean coverage\t{:.2}",
        data_yield as f64 / genome_size as f64
    );
    let data_yield_long = lengths.iter().filter(|l| l > &&25000).sum::<u128>();
    println!("Yield [Gb] (>25kb)\t{:.2}", data_yield_long as f64 / 1e9);
    println!("N50\t{}", calculations::get_n(lengths, data_yield, 0.50));
    println!("N75\t{}", calculations::get_n(lengths, data_yield, 0.75));
    println!("Median length\t{:.2}", calculations::median_length(lengths));
    println!("Mean length\t{:.2}", data_yield / num_reads as u128);
    if let Some(identities) = identities {
        println!("Median identity\t{:.2}", calculations::median(identities));
        println!(
            "Mean identity\t{:.2}",
            identities.iter().sum::<f64>() / (num_reads as f64)
        );
        // modal accuracy has lower precision because it gets inflated and divided by 10, losing everything after the first decimal
        println!(
            "Modal identity\t{:.1}",
            calculations::modal_accuracy(identities)
        );
    }
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
        hist: true,
        checksum: false,
        arrow: Some("test.feather".to_string()),
        karyotype: true,
        phased: true,
        spliced: false,
        ubam: false,
    };
    let (metrics, header) = extract_from_bam::extract(&args);
    assert!(metrics_from_bam(metrics, args, header).is_ok())
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
        hist: false,
        checksum: false,
        arrow: None,
        karyotype: false,
        phased: false,
        spliced: false,
        ubam: false,
    };
    let (metrics, header) = extract_from_bam::extract(&args);
    assert!(metrics_from_bam(metrics, args, header).is_ok())
}

#[test]
fn extract_ubam() {
    let args = Cli {
        input: "test-data/small-test-ubam.bam".to_string(),
        threads: 8,
        reference: None,
        min_read_len: 0,
        hist: true,
        checksum: false,
        arrow: Some("test.feather".to_string()),
        karyotype: false,
        phased: false,
        spliced: false,
        ubam: true,
    };
    let (metrics, header) = extract_from_bam::extract(&args);
    assert!(metrics_from_bam(metrics, args, header).is_ok())
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
        hist: true,
        checksum: false,
        arrow: None,
        karyotype: false,
        phased: false,
        spliced: false,
        ubam: false,
    };
    let (metrics, header) = extract_from_bam::extract(&args);
    assert!(metrics_from_bam(metrics, args, header).is_ok())
}
