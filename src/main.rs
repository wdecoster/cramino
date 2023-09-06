use clap::AppSettings::DeriveDisplayOrder;
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
#[structopt(global_settings=&[DeriveDisplayOrder])]
#[clap(author, version, about="Tool to extract QC metrics from cram or bam", long_about = None)]
pub struct Cli {
    /// cram or bam file to check
    #[clap(value_parser, validator=is_file)]
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

    /// If a checksum has to be calculated
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
    if pathname == "-" {
        return Ok(());
    }
    let path = PathBuf::from(pathname);
    if path.is_file() {
        Ok(())
    } else {
        Err(format!("Input file {} is invalid", path.display()))
    }
}

fn main() {
    env_logger::init();
    let mut args = Cli::parse();
    if args.ubam {
        args.karyotype = false;
        args.phased = false;
        args.spliced = false;
    };
    info!("Collected arguments");
    let metrics = extract_from_bam::extract(&args);

    metrics_from_bam(metrics, args);
    info!("Finished");
}

fn metrics_from_bam(metrics: Data, args: Cli) {
    let bam = file_info::BamFile { path: args.input };
    println!("File name\t{}", bam.file_name());

    generate_main_output(
        metrics.lengths.as_ref().unwrap(),
        metrics.identities.as_ref(),
        utils::get_genome_size(&bam.path),
    );

    println!("Path\t{}", bam);
    println!("Creation time\t{}", bam.file_time());
    if args.checksum {
        println!("Checksum\t{}", bam.checksum());
    }

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
        karyotype::make_karyotype(metrics.tids.as_ref().unwrap(), bam.to_string());
    }
    if args.spliced {
        splicing::splice_metrics(metrics.exons.unwrap());
    }
    if args.hist {
        histograms::make_histogram_lengths(metrics.lengths.as_ref().unwrap());
        if !args.ubam {
            histograms::make_histogram_identities(metrics.identities.as_ref().unwrap());
        }
        if args.phased {
            histograms::make_histogram_phaseblocks(&phaseblocks.unwrap())
        }
    }
}

fn generate_main_output(lengths: &Vec<u64>, identities: Option<&Vec<f64>>, genome_size: u64) {
    let num_reads = lengths.len();
    if num_reads < 2 {
        error!("Not enough reads to calculate metrics!");
        panic!();
    }
    let data_yield: u64 = lengths.iter().sum::<u64>();
    println!("Number of reads\t{num_reads}");
    println!("Yield [Gb]\t{:.2}", data_yield as f64 / 1e9);
    println!(
        "Mean coverage\t{:.2}",
        data_yield as f64 / genome_size as f64
    );
    let data_yield_long = lengths.iter().filter(|l| l > &&25000).sum::<u64>();
    println!("Yield [Gb] (>25kb)\t{:.2}", data_yield_long as f64 / 1e9);
    println!("N50\t{}", calculations::get_n50(lengths, data_yield));
    println!("N75\t{}", calculations::get_n75(lengths, data_yield));
    println!("Median length\t{:.2}", calculations::median_length(lengths));
    println!("Mean length\t{:.2}", data_yield / num_reads as u64);
    if let Some(identities) = identities {
        println!("Median identity\t{:.2}", calculations::median(identities));
        println!(
            "Mean identity\t{:.2}",
            identities.iter().sum::<f64>() / (num_reads as f64)
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
        checksum: true,
        arrow: Some("test.feather".to_string()),
        karyotype: true,
        phased: true,
        spliced: false,
        ubam: false,
    };
    let metrics = extract_from_bam::extract(&args);
    metrics_from_bam(metrics, args)
}

#[test]
fn extract_ubam() {
    let args = Cli {
        input: "test-data/small-test-phased.bam".to_string(),
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
    let metrics = extract_from_bam::extract(&args);
    metrics_from_bam(metrics, args)
}
