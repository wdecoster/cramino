use clap::AppSettings::DeriveDisplayOrder;
use clap::Parser;
use log::{error, info};
pub mod calculations;
pub mod extract_from_bam;
pub mod file_info;

// The arguments end up in the Cli struct
#[derive(Parser, Debug)]
#[structopt(global_settings=&[DeriveDisplayOrder])]
#[clap(author, version, about="Tool to extract QC metrics from bam", long_about = None)]
struct Cli {
    /// bam file to check
    #[clap(value_parser)]
    bam: String,

    /// Number of parallel threads to use
    #[clap(short, long, value_parser, default_value_t = 8)]
    threads: usize,
}

fn main() {
    env_logger::init();
    let args = Cli::parse();
    info!("Collected arguments");
    metrics_from_bam(args.bam, args.threads);
    info!("Finished");
}

fn metrics_from_bam(bam: String, threads: usize) {
    let (mut lengths, mut identities): (Vec<u32>, Vec<f32>) =
        extract_from_bam::extract(&bam, threads);
    if lengths.len() < 2 {
        error!("Not enough reads to calculate metrics!");
        panic!();
    }
    lengths.sort_unstable();
    identities.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());
    let data_yield: u32 = lengths.iter().sum::<u32>();
    let num_reads = lengths.len();
    let n50 = calculations::get_n50(&lengths, data_yield);
    let median_length = calculations::median(&lengths);
    let median_pid = calculations::median(&identities);
    // get histograms for pid and length
    let bam = file_info::BamFile { path: bam };
    let bam_name = bam.file_name();
    let checksum = bam.checksum();
    let file_time = bam.file_time();
    println!("bam\tnumber_of_reads\tyield\tn50\tmedian_length\tmedian_pid\thistogram_length\thistogram_pid\tbam_path\tbam_checksum\tbam_created");
    println!(
        "{bam_name}\t{num_reads}\t{data_yield}\t{n50}\t{median_length}\t{median_pid}\t{bam}\t{checksum}\t{file_time}"
    )
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
    metrics_from_bam("/home/wdecoster/test-data/test.bam".to_string(), 8)
}
