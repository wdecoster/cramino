use clap::AppSettings::DeriveDisplayOrder;
use clap::{Parser, Subcommand};
use log::{error, info};
pub mod calculations;
pub mod extract_from_bam;
pub mod file_info;

// The arguments end up in the Cli struct
#[derive(Parser, Debug)]
#[structopt(global_settings=&[DeriveDisplayOrder])]
#[clap(author, version, about="Tool to extract QC metrics from bam by sampling", long_about = None)]
struct Cli {
    #[clap(subcommand)]
    command: Commands,
}
// Every subcommand is a variation of the Commands Enum, and has its arguments defined below
#[derive(Debug, Subcommand)]
enum Commands {
    /// extract metrics
    #[clap(arg_required_else_help = true)]
    Extract {
        /// bam file to check
        #[clap(value_parser)]
        bam: String,

        /// fraction of reads to use
        #[clap(short, long, value_parser, default_value_t = 0.01)]
        fraction: f32,

        /// Number of parallel threads to use
        #[clap(short, long, value_parser, default_value_t = 8)]
        threads: usize,
    },
    Combine {},
}

fn main() {
    env_logger::init();
    let args = Cli::parse();
    info!("Collected arguments");
    match args.command {
        Commands::Extract {
            bam,
            fraction,
            threads,
        } => extract_metrics(bam, fraction, threads),
        Commands::Combine {} => {
            unimplemented!();
        }
    }
}

///extract metrics from a bam file by subsampling reads
fn extract_metrics(bamp: String, fraction: f32, threads: usize) {
    // the fraction determines how many reads to skip every time by inverting it
    // every "inverse"th read is used
    let inverse = (1.0 / fraction) as usize;
    let (mut lengths, mut pids): (Vec<u32>, Vec<f32>) =
        extract_from_bam::extract(&bamp, inverse, threads);
    // to get the yield, sum all lengths of the selected reads and multiply that again by the inverse of the fraction
    if lengths.len() < 2 {
        error!("Not enough reads to calculate metrics!");
        panic!();
    }
    lengths.sort_unstable();
    pids.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());
    let data_yield: u32 = lengths.iter().sum::<u32>() * inverse as u32;
    let numreads = lengths.len() * inverse;
    let n50 = calculations::get_n50(&lengths, data_yield);
    let median_length = calculations::median(&lengths);
    let median_pid = calculations::median(&pids);
    // get histograms for pid and length
    let bamname = file_info::file_name(&bamp);
    let checksum = file_info::checksum(&bamp);
    let file_time = file_info::file_time(&bamp);
    println!("bam\tnumber_of_reads\tyield\tn50\tmedian_length\tmedian_pid\thistogram_length\thistogram_pid\tbam_path\tbam_checksum\tbam_created");
    println!(
        "{bamname}\t{numreads}\t{data_yield}\t{n50}\t{median_length}\t{median_pid}\t{bamp}\t{checksum}\t{file_time}"
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
    extract_metrics("~/test-data/test.bam".to_string(), 0.01, 8)
}
