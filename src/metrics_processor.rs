use crate::{
    Cli, calculations, extract_from_bam::Data, file_info, histograms, metrics, phased, utils,
};
use clap::builder::{TypedValueParser, ValueParserFactory};
use rust_htslib::bam;
use std::collections::HashMap;
use std::fmt;
use std::str::FromStr;

#[derive(Debug, PartialEq, Clone, Copy)]
pub enum OutputFormat {
    Text,
    Json,
    Tsv,
}

// Implement Display for pretty printing
impl fmt::Display for OutputFormat {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            OutputFormat::Text => write!(f, "text"),
            OutputFormat::Json => write!(f, "json"),
            OutputFormat::Tsv => write!(f, "tsv"),
        }
    }
}

// Implement FromStr for parsing from command line
impl FromStr for OutputFormat {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "text" => Ok(OutputFormat::Text),
            "json" => Ok(OutputFormat::Json),
            "tsv" => Ok(OutputFormat::Tsv),
            _ => Err(format!("Unknown output format: {}", s)),
        }
    }
}

// Implement ValueParserFactory for clap integration
impl ValueParserFactory for OutputFormat {
    type Parser = OutputFormatValueParser;

    fn value_parser() -> Self::Parser {
        OutputFormatValueParser
    }
}

#[derive(Clone, Debug)]
pub struct OutputFormatValueParser;

impl TypedValueParser for OutputFormatValueParser {
    type Value = OutputFormat;

    fn parse_ref(
        &self,
        cmd: &clap::Command,
        arg: Option<&clap::Arg>,
        value: &std::ffi::OsStr,
    ) -> Result<Self::Value, clap::Error> {
        let value_str = value
            .to_str()
            .ok_or_else(|| clap::Error::new(clap::error::ErrorKind::InvalidUtf8).with_cmd(cmd))?;

        OutputFormat::from_str(value_str).map_err(|_err| {
            let mut err = clap::Error::new(clap::error::ErrorKind::InvalidValue).with_cmd(cmd);
            if let Some(arg) = arg {
                err.insert(
                    clap::error::ContextKind::InvalidArg,
                    clap::error::ContextValue::String(arg.to_string()),
                );
            }
            err.insert(
                clap::error::ContextKind::InvalidValue,
                clap::error::ContextValue::String(value_str.to_string()),
            );
            err.insert(
                clap::error::ContextKind::SuggestedValue,
                clap::error::ContextValue::Strings(vec![
                    "text".to_string(),
                    "json".to_string(),
                    "tsv".to_string(),
                ]),
            );
            err
        })
    }
}

pub fn process_metrics(
    metrics_data: Data,
    args: &Cli,
    header: rust_htslib::bam::Header,
) -> Result<(), Box<dyn std::error::Error>> {
    let bam = file_info::BamFile {
        path: args.input.clone(),
    };

    // Create a metrics object
    let mut metrics_obj = metrics::Metrics::new(metrics::FileInfo {
        name: bam.file_name(),
        path: bam.to_string(),
        creation_time: bam.file_time(),
    });

    let lengths = metrics_data.lengths.as_ref().unwrap();
    let hist_requested = args.hist.is_some() || args.hist_count.is_some();

    // Check if no reads passed the filters
    if lengths.is_empty() {
        eprintln!("Warning: No reads pass your filtering criteria");

        // Set minimal metrics with zeros
        metrics_obj.alignment_stats = metrics::AlignmentStats {
            num_alignments: 0,
            percent_from_total: 0.0,
            num_reads: 0,
        };

        metrics_obj.read_stats = metrics::ReadStats {
            yield_gb: 0.0,
            mean_coverage: 0.0,
            yield_gb_long: 0.0,
            n50: 0,
            n75: 0,
            median_length: 0.0,
            mean_length: 0.0,
        };

        // Output based on selected format
        match args.format {
            OutputFormat::Text => {
                crate::text_output::print_text_output(&metrics_obj);
                // Handle --hist-count flag (output empty histogram counts after metrics)
                if let Some(hist_count_file) = &args.hist_count {
                    let value_label = if args.scaled { "bases" } else { "count" };
                    if let Some(file) = hist_count_file {
                        std::fs::write(file, format!("\nbin_start\tbin_end\t{}\n", value_label))?;
                    } else {
                        println!("\nbin_start\tbin_end\t{}", value_label);
                    }
                }
            }
            OutputFormat::Json => {
                if hist_requested {
                    metrics_obj.histograms = Some(histograms::build_histograms(&metrics_data));
                }
                println!("{}", serde_json::to_string_pretty(&metrics_obj).unwrap());
                // Handle --hist-count flag (output empty histogram counts after metrics)
                if let Some(hist_count_file) = &args.hist_count {
                    let value_label = if args.scaled { "bases" } else { "count" };
                    if let Some(file) = hist_count_file {
                        std::fs::write(file, format!("\nbin_start\tbin_end\t{}\n", value_label))?;
                    } else {
                        println!("\nbin_start\tbin_end\t{}", value_label);
                    }
                }
            }
            OutputFormat::Tsv => {
                crate::tsv_output::print_tsv_output(&metrics_obj);
                // Handle --hist-count flag (output empty histogram counts after metrics)
                if let Some(hist_count_file) = &args.hist_count {
                    let value_label = if args.scaled { "bases" } else { "count" };
                    if let Some(file) = hist_count_file {
                        std::fs::write(file, format!("\nbin_start\tbin_end\t{}\n", value_label))?;
                    } else {
                        println!("\nbin_start\tbin_end\t{}", value_label);
                    }
                }
            }
        }

        return Ok(());
    }

    // Continue with normal processing if we have reads
    let num_alignments = lengths.len();
    let num_reads = metrics_data.num_reads;
    let all_alignments = metrics_data.all_counts;

    // Fill in alignment stats
    metrics_obj.alignment_stats = metrics::AlignmentStats {
        num_alignments,
        percent_from_total: (num_reads as f64) / (all_alignments as f64) * 100.0,
        num_reads,
    };

    // Calculate and fill read stats
    let genome_size = utils::get_genome_size(&header)?;
    let (data_yield, data_yield_long) = utils::calculate_data_yield(lengths);

    metrics_obj.read_stats = metrics::ReadStats {
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
        metrics_obj.identity_stats = Some(metrics::IdentityStats {
            median_identity: calculations::median(identities),
            mean_identity: identities.iter().sum::<f64>() / (identities.len() as f64),
            modal_identity: calculations::modal_accuracy(identities),
        });
    }

    // Add phase metrics if available
    let phaseblocks = if args.phased {
        let phaseblocks = phased::phase_metrics(
            metrics_data.tids.as_ref().expect("TIDs data is missing"),
            metrics_data.starts.clone().expect("Starts data is missing"),
            metrics_data.ends.clone().expect("Ends data is missing"),
            metrics_data
                .phasesets
                .as_ref()
                .expect("Phase sets data is missing"),
        );

        if !phaseblocks.is_empty() {
            let phased_bases = phaseblocks.iter().sum::<i64>();
            let phased_reads = metrics_data
                .phasesets
                .as_ref()
                .expect("Phase sets data is missing")
                .iter()
                .filter(|p| p.is_some())
                .count();

            metrics_obj.phase_stats = Some(metrics::PhaseStats {
                fraction_phased: (phased_reads as f32) / (num_reads as f32),
                num_phaseblocks: phaseblocks.len(),
                total_bases_phased_gb: phased_bases as f64 / 1e9,
                median_phaseblock_length: phased::median(&phaseblocks),
                n50_phaseblock_length: phased::get_n50(&phaseblocks, phased_bases),
            });
        }
        Some(phaseblocks)
    } else {
        None
    };

    // Add karyotype data if requested
    if args.karyotype {
        let head_view = bam::HeaderView::from_header(&header);
        let mut tidcount = HashMap::new();

        for tid in metrics_data.tids.as_ref().expect("TIDs data is missing") {
            *tidcount.entry(*tid).or_default() += 1;
        }

        let mut karyotype_data = Vec::new();
        for (tid, count) in tidcount.iter() {
            if *tid >= 0 {
                // Skip unmapped reads (tid = -1)
                let chrom = std::str::from_utf8(
                    head_view.tid2name((*tid).try_into().expect("Failed to convert TID to usize")),
                )
                .unwrap();
                let chrom_length = head_view
                    .target_len((*tid).try_into().expect("Failed to convert TID to usize"))
                    .unwrap();
                let norm_count = (*count as f32) / (chrom_length as f32);

                karyotype_data.push(metrics::ChromosomeData {
                    chromosome: chrom.to_string(),
                    count: *count,
                    normalized_count: norm_count,
                });
            }
        }

        metrics_obj.karyotype_stats = Some(karyotype_data);
    }

    // Add splicing metrics if requested
    if args.spliced && metrics_data.exons.is_some() {
        let exon_counts = metrics_data.exons.as_ref().unwrap();
        let num_reads = exon_counts.len();
        let num_single_exon = exon_counts.iter().filter(|&&x| x == 1).count();

        metrics_obj.splice_stats = Some(metrics::SpliceStats {
            median_exons: calculations::median_splice(exon_counts),
            mean_exons: (exon_counts.iter().sum::<usize>() as f32) / (num_reads as f32),
            fraction_unspliced: (num_single_exon as f32) / (num_reads as f32),
        });
    }

    // Output based on selected format
    match args.format {
        OutputFormat::Text => {
            // Print text output using the collected metrics
            crate::text_output::print_text_output(&metrics_obj);
            if let Some(hist_file) = &args.hist {
                histograms::create_histograms(&metrics_data, hist_file, phaseblocks, args.scaled)?;
            }
            // Handle --hist-count flag (output histogram counts after metrics)
            if let Some(hist_count_file) = &args.hist_count {
                histograms::output_histogram_counts(&metrics_data, hist_count_file, args.scaled)?;
            }
        }
        OutputFormat::Json => {
            if hist_requested {
                metrics_obj.histograms = Some(histograms::build_histograms(&metrics_data));
            }
            println!("{}", serde_json::to_string_pretty(&metrics_obj).unwrap());
            if let Some(hist_file) = &args.hist {
                histograms::create_histograms(&metrics_data, hist_file, phaseblocks, args.scaled)?;
            }
            // Handle --hist-count flag (output histogram counts after metrics)
            if let Some(hist_count_file) = &args.hist_count {
                histograms::output_histogram_counts(&metrics_data, hist_count_file, args.scaled)?;
            }
        }
        OutputFormat::Tsv => {
            crate::tsv_output::print_tsv_output(&metrics_obj);
            if let Some(hist_file) = &args.hist {
                histograms::create_histograms(&metrics_data, hist_file, phaseblocks, args.scaled)?;
            }
            // Handle --hist-count flag (output histogram counts after metrics)
            if let Some(hist_count_file) = &args.hist_count {
                histograms::output_histogram_counts(&metrics_data, hist_count_file, args.scaled)?;
            }
        }
    }

    Ok(())
}
