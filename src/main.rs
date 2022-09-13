use clap::AppSettings::DeriveDisplayOrder;
use clap::Parser;
use histo_fp::Histogram;
use log::{error, info};
use std::{cmp::Ordering, path::PathBuf};

pub mod calculations;
pub mod extract_from_bam;
pub mod file_info;

// The arguments end up in the Cli struct
#[derive(Parser, Debug)]
#[structopt(global_settings=&[DeriveDisplayOrder])]
#[clap(author, version, about="Tool to extract QC metrics from bam", long_about = None)]
struct Cli {
    /// bam file to check
    #[clap(value_parser, validator=is_file)]
    bam: String,

    /// Number of parallel threads to use
    #[clap(short, long, value_parser, default_value_t = 4)]
    threads: usize,

    /// If histograms have to be generated
    #[clap(short, long, value_parser)]
    hist: bool,
}

fn is_file(pathname: &str) -> Result<(), String> {
    let path = PathBuf::from(pathname);
    if path.is_file() {
        Ok(())
    } else {
        Err(format!("Input file {} is invalid", path.display()))
    }
}

fn main() {
    env_logger::init();
    let args = Cli::parse();
    info!("Collected arguments");
    metrics_from_bam(args.bam, args.threads, args.hist);
    info!("Finished");
}

fn metrics_from_bam(bam: String, threads: usize, hist: bool) {
    let (mut lengths, mut identities): (Vec<u64>, Vec<f32>) =
        extract_from_bam::extract(&bam, threads);
    let num_reads = lengths.len();
    if num_reads < 2 {
        error!("Not enough reads to calculate metrics!");
        panic!();
    }
    lengths.sort_unstable();
    identities.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());
    let data_yield: u64 = lengths.iter().sum::<u64>();
    let bam = file_info::BamFile { path: bam };
    println!("File name\t{}", bam.file_name());
    println!("Number of reads\t{num_reads}");
    println!("Yield [Gb]\t{:.2}", data_yield as f64 / 1e9);
    println!("N50\t{}", calculations::get_n50(&lengths, data_yield));
    println!(
        "Median length\t{:.2}",
        calculations::median_length(&lengths)
    );
    println!("Mean length\t{:.2}", data_yield / num_reads as u64);
    println!("Median identity\t{:.2}", calculations::median(&identities));
    println!(
        "Mean identity\t{:.2}",
        identities.iter().sum::<f32>() / num_reads as f32
    );
    println!("Path\t{}", bam);
    println!("Checksum\t{}", bam.checksum());
    println!("Creation time\t{}", bam.file_time());
    if hist {
        println!(
            "\n\nHistogram for lengths\n{}",
            make_histogram_lengths(lengths)
        );
        println!(
            "\n\nHistogram for identities\n{}",
            make_histogram_identities(identities)
        );
    }
}

fn make_histogram_lengths(array: Vec<u64>) -> Histogram {
    let mut histogram = Histogram::with_buckets(100, Some(2));
    for value in array.into_iter() {
        histogram.add(value as f64);
    }

    histogram
}

fn make_histogram_identities(array: Vec<f32>) -> Histogram {
    let bins = 100.0
        - array
            .iter()
            .min_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Less))
            .unwrap();
    let mut histogram = Histogram::with_buckets(bins as u64, None);
    for value in array.into_iter() {
        histogram.add(value as f64);
    }

    histogram
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
    metrics_from_bam("/home/wdecoster/test-data/test.bam".to_string(), 8, true)
}
