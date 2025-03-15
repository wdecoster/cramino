use clap::Parser;
use extract_from_bam::Data;
use log::info;
use rust_htslib::bam;
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
pub mod metrics;

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

    /// Output results in JSON format
    #[clap(long, value_parser)]
    json: bool,
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
    metrics_data: Data,
    args: Cli,
    header: rust_htslib::bam::Header,
) -> Result<(), rust_htslib::errors::Error> {
    let bam = file_info::BamFile { path: args.input.clone() };
    
    // Create a metrics object
    let mut metrics = metrics::Metrics::new(metrics::FileInfo {
        name: bam.file_name(),
        path: bam.to_string(),
        creation_time: bam.file_time(),
    });
    
    let lengths = metrics_data.lengths.as_ref().unwrap();
    let num_alignments = lengths.len();
    let num_reads = metrics_data.num_reads;
    let all_alignments = metrics_data.all_counts;
    
    // Fill in alignment stats
    metrics.alignment_stats = metrics::AlignmentStats {
        num_alignments,
        percent_from_total: (num_reads as f64) / (all_alignments as f64) * 100.0,
        num_reads,
    };
    
    // Calculate and fill read stats
    let genome_size = utils::get_genome_size(&header)?;
    let (data_yield, data_yield_long) = calculate_data_yield(lengths);
    
    metrics.read_stats = metrics::ReadStats {
        yield_gb: data_yield as f64 / 1e9,
        mean_coverage: data_yield as f64 / genome_size as f64,
        yield_gb_long: data_yield_long as f64 / 1e9,
        n50: calculations::get_n(lengths, data_yield, 0.50),
        n75: calculations::get_n(lengths, data_yield, 0.75),
        median_length: calculations::median_length(lengths),
        mean_length: data_yield as f64 / lengths.len() as f64,
    };
    
    // Add identity metrics if available
    if let Some(identities) = metrics_data.identities.as_ref() {
        metrics.identity_stats = Some(metrics::IdentityStats {
            median_identity: calculations::median(identities),
            mean_identity: identities.iter().sum::<f64>() / (identities.len() as f64),
            modal_identity: calculations::modal_accuracy(identities),
        });
    }
    
    // Add phase metrics if available
    if args.phased {
        let phaseblocks = phased::phase_metrics(
            metrics_data.tids.as_ref().unwrap(),
            metrics_data.starts.clone().unwrap(),
            metrics_data.ends.clone().unwrap(),
            metrics_data.phasesets.as_ref().unwrap(),
        );
        
        if !phaseblocks.is_empty() {
            let phased_bases = phaseblocks.iter().sum::<i64>();
            let phased_reads = metrics_data.phasesets.as_ref().unwrap().iter()
                .filter(|p| p.is_some())
                .count();
                
            metrics.phase_stats = Some(metrics::PhaseStats {
                fraction_phased: (phased_reads as f32) / (num_reads as f32),
                num_phaseblocks: phaseblocks.len(),
                total_bases_phased_gb: phased_bases as f64 / 1e9,
                median_phaseblock_length: phased::median(&phaseblocks),
                n50_phaseblock_length: phased::get_n50(&phaseblocks, phased_bases),
            });
        }
    }
    
    // Add karyotype data if requested
    if args.karyotype {
        let head_view = bam::HeaderView::from_header(&header);
        let mut tidcount = std::collections::HashMap::new();
        
        for tid in metrics_data.tids.as_ref().unwrap() {
            *tidcount.entry(*tid).or_default() += 1;
        }
        
        let mut karyotype_data = Vec::new();
        for (tid, count) in tidcount.iter() {
            if *tid >= 0 {  // Skip unmapped reads (tid = -1)
                let chrom = std::str::from_utf8(head_view.tid2name((*tid).try_into().unwrap())).unwrap();
                let chrom_length = head_view.target_len((*tid).try_into().unwrap()).unwrap();
                let norm_count = (*count as f32) / (chrom_length as f32);
                
                karyotype_data.push(metrics::ChromosomeData {
                    chromosome: chrom.to_string(),
                    count: *count,
                    normalized_count: norm_count,
                });
            }
        }
        
        metrics.karyotype_stats = Some(karyotype_data);
    }
    
    // Add splicing metrics if requested
    if args.spliced && metrics_data.exons.is_some() {
        let exon_counts = metrics_data.exons.as_ref().unwrap();
        let num_reads = exon_counts.len();
        let num_single_exon = exon_counts.iter().filter(|&&x| x == 1).count();
        
        metrics.splice_stats = Some(metrics::SpliceStats {
            median_exons: calculations::median_splice(exon_counts),
            mean_exons: (exon_counts.iter().sum::<usize>() as f32) / (num_reads as f32),
            fraction_unspliced: (num_single_exon as f32) / (num_reads as f32),
        });
    }
    
    // Output either as JSON or text format
    if args.json {
        println!("{}", serde_json::to_string_pretty(&metrics).unwrap());
    } else {
        // Print text output using the collected metrics
        print_text_output(&metrics);
        
        // Generate histograms if requested
        if args.hist {
            histograms::make_histogram_lengths(lengths);
            
            if !args.ubam && metrics_data.identities.is_some() {
                histograms::make_histogram_identities(metrics_data.identities.as_ref().unwrap());
            }
            
            if args.phased && metrics.phase_stats.is_some() {
                let phaseblocks = phased::phase_metrics(
                    metrics_data.tids.as_ref().unwrap(),
                    metrics_data.starts.unwrap(),
                    metrics_data.ends.unwrap(),
                    metrics_data.phasesets.as_ref().unwrap(),
                );
                histograms::make_histogram_phaseblocks(&phaseblocks);
            }
            
            if args.spliced && metrics_data.exons.is_some() {
                histograms::make_histogram_exons(metrics_data.exons.as_ref().unwrap());
            }
        }
    }
    
    Ok(())
}

fn print_text_output(metrics: &metrics::Metrics) {
    // Print file info
    println!("File name\t{}", metrics.file_info.name);
    
    // Print alignment stats
    println!("Number of alignments\t{}", metrics.alignment_stats.num_alignments);
    println!(
        "% from total alignments\t{:.2}",
        metrics.alignment_stats.percent_from_total
    );
    println!("Number of reads\t{}", metrics.alignment_stats.num_reads);
    
    // Print read stats
    println!("Yield [Gb]\t{:.2}", metrics.read_stats.yield_gb);
    println!("Mean coverage\t{:.2}", metrics.read_stats.mean_coverage);
    println!("Yield [Gb] (>25kb)\t{:.2}", metrics.read_stats.yield_gb_long);
    println!("N50\t{}", metrics.read_stats.n50);
    println!("N75\t{}", metrics.read_stats.n75);
    println!("Median length\t{:.2}", metrics.read_stats.median_length);
    println!("Mean length\t{:.2}", metrics.read_stats.mean_length);
    
    // Print file path and creation time
    println!("Path\t{}", metrics.file_info.path);
    println!("Creation time\t{}", metrics.file_info.creation_time);
    
    // Print identity stats if available
    if let Some(identity_stats) = &metrics.identity_stats {
        println!("Median identity\t{:.2}", identity_stats.median_identity);
        println!("Mean identity\t{:.2}", identity_stats.mean_identity);
        println!("Modal identity\t{:.1}", identity_stats.modal_identity);
    }
    
    // Print phase stats if available
    if let Some(phase_stats) = &metrics.phase_stats {
        println!("Fraction reads phased\t{:.2}", phase_stats.fraction_phased);
        println!("Number of phaseblocks\t{}", phase_stats.num_phaseblocks);
        println!("Total bases phased [Gb]\t{:.2}", phase_stats.total_bases_phased_gb);
        println!("Median phaseblock length\t{:.2}", phase_stats.median_phaseblock_length);
        println!("N50 phaseblock length\t{}", phase_stats.n50_phaseblock_length);
    }
    
    // Print karyotype stats if available
    if let Some(karyotype_stats) = &metrics.karyotype_stats {
        if !karyotype_stats.is_empty() {
            // Calculate median for normalization
            let counts: Vec<f32> = karyotype_stats.iter().map(|c| c.normalized_count).collect();
            let median_count = if !counts.is_empty() {
                let mut counts_clone = counts.clone();
                counts_clone.sort_by(|a, b| a.partial_cmp(b).unwrap());
                counts_clone[counts_clone.len() / 2]
            } else {
                1.0  // Default if no data
            };
            
            println!("\n\n# Normalized read count per chromosome\n");
            let mut sorted_stats = karyotype_stats.clone();
            sorted_stats.sort_by(|a, b| a.chromosome.cmp(&b.chromosome));
            
            for chrom_data in sorted_stats {
                println!(
                    "{}\t{:.2}", 
                    chrom_data.chromosome,
                    chrom_data.normalized_count / median_count
                );
            }
        } else {
            println!("\n\n# Warning - no contigs found in BAM file!\n");
        }
    }
    
    // Print splice stats if available
    if let Some(splice_stats) = &metrics.splice_stats {
        println!("Median number of exons\t{}", splice_stats.median_exons);
        println!("Mean number of exons\t{:.2}", splice_stats.mean_exons);
        println!("Fraction unspliced reads\t{:.2}", splice_stats.fraction_unspliced);
    }
}

// Helper function to calculate data yield
fn calculate_data_yield(lengths: &[u128]) -> (u128, u128) {
    lengths.iter().fold((0u128, 0u128), |(total, long), &len| {
        let long_increment = if len > 25000 { len } else { 0 };
        (total + len, long + long_increment)
    })
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
        json: false,
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
        json: false,
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
        json: false,
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
        json: false,
    };
    let (metrics, header) = extract_from_bam::extract(&args);
    assert!(metrics_from_bam(metrics, args, header).is_ok())
}

#[test]
fn extract_json() {
    let args = Cli {
        input: "test-data/small-test-phased.bam".to_string(),
        threads: 8,
        reference: None,
        min_read_len: 0,
        hist: false,
        checksum: false,
        arrow: None,
        karyotype: true,
        phased: true,
        spliced: false,
        ubam: false,
        json: true,
    };
    let (metrics, header) = extract_from_bam::extract(&args);
    assert!(metrics_from_bam(metrics, args, header).is_ok())
}
