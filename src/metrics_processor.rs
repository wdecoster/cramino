use crate::{
    Cli, calculations, extract_from_bam::Data, file_info, histograms, metrics, phased, utils,
};
use clap::builder::{TypedValueParser, ValueParserFactory};
use rust_htslib::bam;
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
    // With json output, histograms without a file are only included in the json,
    // as anything else on stdout would make it invalid json
    let json = matches!(args.format, OutputFormat::Json);
    let hist = args.hist.as_ref().filter(|file| !json || file.is_some());
    let hist_count = args
        .hist_count
        .as_ref()
        .filter(|file| !json || file.is_some());

    let genome_size = utils::get_genome_size(&header)?;

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
            // the coverage is 0 rather than NA if it can be calculated, i.e. with reference sequences
            mean_coverage: (genome_size > 0).then_some(0.0),
            ..Default::default()
        };

        // Output based on selected format
        match args.format {
            OutputFormat::Text => {
                crate::text_output::print_text_output(&metrics_obj);
                // Handle --hist-count flag (output empty histogram counts after metrics)
                if let Some(hist_count_file) = hist_count {
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
                if let Some(hist_count_file) = hist_count {
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
                if let Some(hist_count_file) = hist_count {
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
    let all_reads = metrics_data.all_reads;

    // Fill in alignment stats
    metrics_obj.alignment_stats = metrics::AlignmentStats {
        num_alignments,
        percent_from_total: (num_reads as f64) / (all_reads as f64) * 100.0,
        num_reads,
    };

    // Calculate and fill read stats
    let (data_yield, data_yield_long) = utils::calculate_data_yield(lengths);
    let read_yield = metrics_data
        .read_lengths
        .iter()
        .map(|l| *l as u128)
        .sum::<u128>();

    let read_lengths = &metrics_data.read_lengths;

    metrics_obj.read_stats = metrics::ReadStats {
        yield_gb: data_yield as f64 / 1e9,
        mean_coverage: (genome_size > 0).then(|| data_yield as f64 / genome_size as f64),
        yield_gb_long: data_yield_long as f64 / 1e9,
        n50: calculations::get_n(read_lengths, read_yield, 0.50),
        n75: calculations::get_n(read_lengths, read_yield, 0.75),
        // there can be no reads left when only supplementary alignments pass --min-read-len
        median_length: calculations::median(read_lengths),
        mean_length: if read_lengths.is_empty() {
            0.0
        } else {
            read_yield as f64 / read_lengths.len() as f64
        },
        n50_aligned: calculations::get_n(lengths, data_yield, 0.50),
        n75_aligned: calculations::get_n(lengths, data_yield, 0.75),
        median_length_aligned: calculations::median(lengths),
        mean_length_aligned: data_yield as f64 / lengths.len() as f64,
    };

    // Add identity metrics if available, which are estimated from the base qualities with --ubam.
    // They are absent if no read has an identity, e.g. without base qualities
    if let Some(identities) = metrics_data.identities.as_ref().filter(|i| !i.is_empty()) {
        let median = calculations::median(identities);
        let mean = identities.iter().sum::<f64>() / (identities.len() as f64);
        let modal = calculations::modal_accuracy(identities);
        if metrics_data.is_ubam {
            metrics_obj.estimated_identity_stats = Some(metrics::EstimatedIdentityStats {
                median_estimated_identity: median,
                mean_estimated_identity: mean,
                modal_estimated_identity: modal,
            });
        } else {
            metrics_obj.identity_stats = Some(metrics::IdentityStats {
                median_identity: median,
                mean_identity: mean,
                modal_identity: modal,
            });
        }
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
            // the phaseblocks are collected in genomic order, while both the median and
            // the N50 need them sorted by length. The original vector is left untouched
            // for the histogram, which bins the lengths.
            let mut phaseblocks_by_length = phaseblocks.clone();
            phaseblocks_by_length.sort_unstable_by(|a, b| b.cmp(a));
            metrics_obj.phase_stats = Some(metrics::PhaseStats {
                // supplementary alignments are not counted, so that this is a fraction of reads.
                // There are no reads left if only supplementary alignments pass --min-read-len
                fraction_phased: if num_reads == 0 {
                    0.0
                } else {
                    (metrics_data.num_phased_reads as f32) / (num_reads as f32)
                },
                num_phaseblocks: phaseblocks.len(),
                total_bases_phased_gb: phased_bases as f64 / 1e9,
                median_phaseblock_length: phased::median(&phaseblocks_by_length),
                n50_phaseblock_length: phased::get_n50(&phaseblocks_by_length, phased_bases),
            });
        }
        Some(phaseblocks)
    } else {
        None
    };

    // Add karyotype data if requested
    if args.karyotype {
        metrics_obj.karyotype_stats = Some(karyotype(
            metrics_data
                .karyotype_tids
                .as_ref()
                .expect("Karyotype TIDs data is missing"),
            &header,
        ));
    }

    // Add splicing metrics if requested
    if args.spliced && metrics_data.exons.as_ref().is_some_and(|e| !e.is_empty()) {
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
            if let Some(hist_file) = hist {
                histograms::create_histograms(&metrics_data, hist_file, phaseblocks, args.scaled)?;
            }
            // Handle --hist-count flag (output histogram counts after metrics)
            if let Some(hist_count_file) = hist_count {
                histograms::output_histogram_counts(&metrics_data, hist_count_file, args.scaled)?;
            }
        }
        OutputFormat::Json => {
            if hist_requested {
                metrics_obj.histograms = Some(histograms::build_histograms(&metrics_data));
            }
            println!("{}", serde_json::to_string_pretty(&metrics_obj).unwrap());
            if let Some(hist_file) = hist {
                histograms::create_histograms(&metrics_data, hist_file, phaseblocks, args.scaled)?;
            }
            // Handle --hist-count flag (output histogram counts after metrics)
            if let Some(hist_count_file) = hist_count {
                histograms::output_histogram_counts(&metrics_data, hist_count_file, args.scaled)?;
            }
        }
        OutputFormat::Tsv => {
            crate::tsv_output::print_tsv_output(&metrics_obj);
            if let Some(hist_file) = hist {
                histograms::create_histograms(&metrics_data, hist_file, phaseblocks, args.scaled)?;
            }
            // Handle --hist-count flag (output histogram counts after metrics)
            if let Some(hist_count_file) = hist_count {
                histograms::output_histogram_counts(&metrics_data, hist_count_file, args.scaled)?;
            }
        }
    }

    Ok(())
}

/// Number of reads (primary alignments) per chromosome, for all chromosomes in the order of
/// the header. The normalized count is the number of reads per bp of the chromosome, divided
/// by the median of that over the chromosomes with reads, so that a chromosome with the
/// typical coverage has a normalized count of 1, and one without reads (e.g. chrY) of 0.
fn karyotype(tids: &[i32], header: &bam::Header) -> Vec<metrics::ChromosomeData> {
    let head_view = bam::HeaderView::from_header(header);
    let mut counts = vec![0usize; head_view.target_count() as usize];
    // unmapped reads (tid -1) are skipped
    for tid in tids.iter().filter_map(|tid| usize::try_from(*tid).ok()) {
        counts[tid] += 1;
    }
    let per_bp: Vec<(usize, usize, f32)> = counts
        .into_iter()
        .enumerate()
        .map(|(tid, count)| {
            let length = head_view.target_len(tid as u32).unwrap();
            (tid, count, count as f32 / length as f32)
        })
        .collect();
    // chromosomes without reads are left out of the median, as there can be many
    // (e.g. alt contigs), which would otherwise make the median 0
    let mut sorted_per_bp: Vec<f32> = per_bp
        .iter()
        .filter(|(_, count, _)| *count > 0)
        .map(|(_, _, c)| *c)
        .collect();
    sorted_per_bp.sort_unstable_by(|a, b| a.total_cmp(b));
    if sorted_per_bp.is_empty() {
        return vec![];
    }
    let median_per_bp = calculations::median(&sorted_per_bp) as f32;
    per_bp
        .into_iter()
        .map(|(tid, count, per_bp)| metrics::ChromosomeData {
            chromosome: String::from_utf8_lossy(head_view.tid2name(tid as u32)).into_owned(),
            count,
            normalized_count: per_bp / median_per_bp,
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::header::HeaderRecord;

    #[test]
    fn test_karyotype_normalized_to_median_in_header_order() {
        let mut header = bam::Header::new();
        for (name, length) in [
            ("chr2", "1000"),
            ("chr10", "2000"),
            ("chrX", "1000"),
            ("chrY", "100"),
        ] {
            header.push_record(
                HeaderRecord::new(b"SQ")
                    .push_tag(b"SN", name)
                    .push_tag(b"LN", length),
            );
        }
        // chr2: 10 per kb, chr10: 10 per kb, chrX: 5 per kb, chrY: none. An unmapped read (-1) is skipped
        let tids: Vec<i32> = [vec![0; 10], vec![1; 20], vec![2; 5], vec![-1]].concat();
        let karyotype = karyotype(&tids, &header);
        let names: Vec<&str> = karyotype.iter().map(|c| c.chromosome.as_str()).collect();
        assert_eq!(names, ["chr2", "chr10", "chrX", "chrY"]);
        let normalized: Vec<f32> = karyotype.iter().map(|c| c.normalized_count).collect();
        assert_eq!(normalized, [1.0, 1.0, 0.5, 0.0]);
        assert_eq!(karyotype[1].count, 20);
    }
}
