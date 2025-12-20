use itertools::Itertools;
use std::cmp::max;
use std::fs::File;
use std::io::{self, Write};

use crate::{extract_from_bam, metrics, utils};

struct LengthHistogramData {
    step: u64,
    max_value: u64,
    counts: Vec<u64>,
    bases: Vec<u128>,
    overflow_count: u64,
    overflow_bases: u128,
}

fn compute_length_histogram_data(array: &[u128]) -> Option<LengthHistogramData> {
    if array.is_empty() {
        return None;
    }

    // dynamically set the maximum value based on the maximum read length, capped at 60k
    let max_read_length = array
        .iter()
        .copied()
        .max()
        .expect("Array is empty, cannot find max");
    let max_value = std::cmp::min(
        60_000u64,
        ((max_read_length as f64 / 10_000.0).ceil() as u64) * 10_000,
    );
    let step = 2000u64;
    let step_count = (max_value / step) as usize;
    let mut counts = vec![0u64; step_count];
    let mut bases = vec![0u128; step_count];
    let mut overflow_count = 0u64;
    let mut overflow_bases = 0u128;

    for &value in array {
        if value >= max_value as u128 {
            overflow_count += 1;
            overflow_bases += value;
        } else {
            let index = (value / step as u128) as usize;
            counts[index] += 1;
            bases[index] += value;
        }
    }

    Some(LengthHistogramData {
        step,
        max_value,
        counts,
        bases,
        overflow_count,
        overflow_bases,
    })
}

fn build_length_histogram(array: &[u128]) -> metrics::Histogram {
    let Some(hist) = compute_length_histogram_data(array) else {
        return metrics::Histogram {
            step: 2000,
            max_value: 0,
            bins: Vec::new(),
        };
    };

    let mut bins = Vec::with_capacity(hist.counts.len() + 1);
    for (index, (count, bases)) in hist.counts.iter().zip(hist.bases.iter()).enumerate() {
        let start = index as u64 * hist.step;
        let end = Some((index as u64 + 1) * hist.step);
        bins.push(metrics::HistogramBin {
            start,
            end,
            count: *count,
            bases: *bases,
        });
    }
    bins.push(metrics::HistogramBin {
        start: hist.max_value,
        end: None,
        count: hist.overflow_count,
        bases: hist.overflow_bases,
    });

    metrics::Histogram {
        step: hist.step,
        max_value: hist.max_value,
        bins,
    }
}

// Q-score bins track read counts and total bases per bin.
fn build_qscore_histogram(hist: &extract_from_bam::QScoreHistogramData) -> metrics::Histogram {
    let step = 1u64;
    let max_value = 40u64;
    let step_count = (max_value / step) as usize;

    let mut bins = Vec::with_capacity(step_count + 1);
    for index in 0..step_count {
        let count = *hist.counts.get(index).unwrap_or(&0);
        let bases = *hist.bases.get(index).unwrap_or(&0);
        bins.push(metrics::HistogramBin {
            start: index as u64 * step,
            end: Some((index as u64 + 1) * step),
            count,
            bases,
        });
    }
    let count = *hist.counts.get(step_count).unwrap_or(&0);
    let bases = *hist.bases.get(step_count).unwrap_or(&0);
    bins.push(metrics::HistogramBin {
        start: max_value,
        end: None,
        count,
        bases,
    });

    metrics::Histogram {
        step,
        max_value,
        bins,
    }
}

pub fn build_histograms(metrics_data: &extract_from_bam::Data) -> metrics::Histograms {
    let read_length = metrics_data
        .lengths
        .as_ref()
        .map(|lengths| build_length_histogram(lengths))
        .unwrap_or(metrics::Histogram {
            step: 2000,
            max_value: 0,
            bins: Vec::new(),
        });
    let q_score = metrics_data
        .q_score_hist
        .as_ref()
        .map(build_qscore_histogram);

    metrics::Histograms {
        read_length,
        q_score,
    }
}

fn output_histogram_counts_tsv<W: Write>(array: &[u128], writer: &mut W, scaled: bool) {
    let Some(hist) = compute_length_histogram_data(array) else {
        return;
    };

    let value_label = if scaled { "bases" } else { "count" };
    // Write TSV header with leading newline for formatting
    writeln!(writer, "\nbin_start\tbin_end\t{}", value_label)
        .expect("Unable to write histogram counts header");

    // Write each bin
    for (index, (count, bases)) in hist.counts.iter().zip(hist.bases.iter()).enumerate() {
        let value = if scaled { *bases } else { (*count).into() };
        writeln!(
            writer,
            "{}\t{}\t{}",
            index as u64 * hist.step,
            (index + 1) as u64 * hist.step,
            value
        )
        .expect("Unable to write histogram counts");
    }

    // Write overflow bin if it has any reads
    if hist.overflow_count > 0 {
        let value = if scaled {
            hist.overflow_bases
        } else {
            hist.overflow_count.into()
        };
        writeln!(writer, "{}+\tNA\t{}", hist.max_value, value)
            .expect("Unable to write overflow bin");
    }
}

// the histograms below are fully defined by the step size and the maximum value
// the step size is the size of each bin
// the maximum value is the value of the last bin
// the last bin (overflow), is printed separately
// the dotsize variable determines how many reads are represented by a single dot
// the functions are duplicated for flexibility in printing the name and formatting the histogram labels
// as well as for future customizations
// in principle it would be possible to enable the user to change the step size or max value, but I don't want to add too many options to the CLI

fn make_histogram_lengths<W: Write>(array: &[u128], writer: &mut W, scaled: bool) {
    let Some(hist) = compute_length_histogram_data(array) else {
        return;
    };
    let overflow = hist.overflow_count as usize;
    let overflow_bp = hist.overflow_bases;
    let dotsize = if scaled {
        let total_bp: u128 = hist.bases.iter().sum::<u128>() + overflow_bp;
        std::cmp::max((total_bp / 500) as usize, 1)
    } else {
        std::cmp::max(array.len() / 500, 1)
    };
    writeln!(
        writer,
        "\n\n# Histogram for read lengths:{}",
        if scaled {
            " (scaled by total basepairs)"
        } else {
            ""
        }
    )
    .expect("Unable to write histogram");
    for (index, (entry, bp)) in hist.counts.iter().zip(hist.bases.iter()).enumerate() {
        let bar = if scaled {
            "∎".repeat(((*bp as usize) / dotsize).max(0))
        } else {
            "∎".repeat(((*entry as usize) / dotsize).max(0))
        };
        writeln!(
            writer,
            "{: >11} {}",
            format!(
                "{}-{}",
                index as u64 * hist.step,
                (index + 1) as u64 * hist.step
            ),
            bar
        )
        .expect("Unable to write histogram");
    }
    if overflow > 0 {
        let bar = if scaled {
            "∎".repeat(((overflow_bp as usize) / dotsize).max(0))
        } else {
            "∎".repeat((overflow / dotsize).max(0))
        };
        writeln!(writer, "{: >11} {}", format!("{}+", hist.max_value), bar)
            .expect("Unable to write histogram");
    }
}

fn make_histogram_identities<W: Write>(array: &[f64], writer: &mut W) {
    let stepsize: u64 = 1;
    let max_value = 40;
    let step_count = max_value / stepsize as usize;
    let mut counts = vec![0; step_count + 1];
    for value in array.iter().map(|x| utils::accuracy_to_phred(*x)) {
        if value < max_value {
            counts[value] += 1;
        }
    }
    // the last bin is for all values above the last step
    counts[step_count] = array.len() - counts.iter().sum::<usize>();
    // the dotsize variable determines how many reads are represented by a single dot
    // I either have to set this dynamically or experiment with it further
    let dotsize = max(array.len() / 500, 1);
    writeln!(writer, "\n\n# Histogram for Phred-scaled accuracies:")
        .expect("Unable to write histogram");
    // print every entry in the vector, except the last one which is done separately
    for (index, entry) in counts.iter().dropping_back(1).enumerate() {
        writeln!(
            writer,
            "{: >6} {}",
            format!(
                "Q{}-{}",
                index as u64 * stepsize,
                (index + 1) * stepsize as usize
            ),
            "∎".repeat(entry / dotsize)
        )
        .expect("Unable to write histogram");
    }
    writeln!(
        writer,
        "{: >6} {}",
        format!("Q{}+", (counts.len() - 1) * stepsize as usize),
        "∎".repeat(counts.last().unwrap() / dotsize)
    )
    .expect("Unable to write histogram");
}

pub fn make_histogram_phaseblocks<W: Write>(array: &[i64], writer: &mut W) {
    // this is a tricky one, as the scale of the length of phaseblocks is hard to predict
    // I may have to increase its max value in the future
    // this configuration seemed sufficient for a randomly picked test file phased with LongShot
    // but presumably other tools can do better, especially with longer reads, and therefore I left some room for longer phase blocks
    let stepsize: i64 = 10000;
    let max_value = 1_000_000;
    let step_count = max_value / stepsize as usize;
    let mut counts = vec![0; step_count + 1];
    for value in array.iter() {
        let index = (value / stepsize) as usize;
        if index < counts.len() {
            counts[index] += 1;
        }
    }
    // the last bin is for all values above the last step
    counts[step_count] = array.len() - counts.iter().sum::<usize>();
    // the dotsize variable determines how many reads are represented by a single dot
    // I either have to set this dynamically or experiment with it further
    let dotsize = max(array.len() / 500, 1);
    writeln!(writer, "\n\n# Histogram for phaseblock lengths:").expect("Unable to write histogram");
    // print every entry in the vector, except the last one which is done separately
    for (index, entry) in counts.iter().dropping_back(1).enumerate() {
        writeln!(
            writer,
            "{: >14} {}",
            format!(
                "{}-{}",
                index as i64 * stepsize,
                (index + 1) * stepsize as usize
            ),
            "∎".repeat(entry / dotsize)
        )
        .expect("Unable to write histogram");
    }
    writeln!(
        writer,
        "{: >14} {}",
        format!("{}+", (counts.len() - 1) * stepsize as usize),
        "∎".repeat(counts.last().unwrap() / dotsize)
    )
    .expect("Unable to write histogram");
}

fn make_histogram_exons<W: Write>(array: &[usize], writer: &mut W) {
    let max_value = 15;
    let step_count = max_value;
    let mut counts = vec![0; step_count + 1];
    for value in array.iter() {
        if *value < max_value {
            counts[*value] += 1;
        }
    }
    // the last bin is for all values above the last step
    counts[step_count] = array.len() - counts.iter().sum::<usize>();
    // the dotsize variable determines how many reads are represented by a single dot
    // I either have to set this dynamically or experiment with it further
    let dotsize = max(array.len() / 500, 1);
    writeln!(writer, "\n\n# Histogram for number of exons:").expect("Unable to write histogram");
    // print the second entry in the vector. The first entry is 0 exons, which is not used (empty)
    // 1 exon is renamed to unspliced
    writeln!(
        writer,
        "{: >9} {}",
        format!("unspliced"),
        "∎".repeat(counts[1] / dotsize)
    )
    .expect("Unable to write histogram");
    // print every entry in the vector, except the first two (first one empty, and second one already done) and last one (done later)
    for (index, entry) in counts.iter().skip(2).dropping_back(1).enumerate() {
        writeln!(
            writer,
            "{: >9} {}",
            format!("{} exons", (index + 2)),
            "∎".repeat(entry / dotsize)
        )
        .expect("Unable to write histogram");
    }
    writeln!(
        writer,
        "{: >9} {}",
        format!("{}+ exons", (counts.len() - 1)),
        "∎".repeat(counts.last().unwrap() / dotsize)
    )
    .expect("Unable to write histogram");
}

pub fn create_histograms(
    metrics_data: &extract_from_bam::Data,
    hist_file: &Option<String>,
    phaseblocks: Option<Vec<i64>>,
    scaled: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut writer: Box<dyn Write> = if let Some(file) = hist_file {
        Box::new(File::create(file)?)
    } else {
        Box::new(io::stdout())
    };
    if let Some(lengths) = &metrics_data.lengths {
        make_histogram_lengths(lengths, &mut writer, scaled);
    }
    if let Some(identities) = &metrics_data.identities {
        make_histogram_identities(identities, &mut writer);
    }
    if let Some(phaseblocks) = phaseblocks {
        make_histogram_phaseblocks(&phaseblocks, &mut writer);
    }
    if let Some(exons) = &metrics_data.exons {
        make_histogram_exons(exons, &mut writer);
    }
    Ok(())
}

pub fn output_histogram_counts(
    metrics_data: &extract_from_bam::Data,
    hist_count_file: &Option<String>,
    scaled: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut writer: Box<dyn Write> = if let Some(file) = hist_count_file {
        Box::new(File::create(file)?)
    } else {
        Box::new(io::stdout())
    };
    if let Some(lengths) = &metrics_data.lengths {
        output_histogram_counts_tsv(lengths, &mut writer, scaled);
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::metrics;

    #[test]
    fn json_histograms_include_bins_with_scaled_and_unscaled_values() {
        let lengths = vec![1000u128, 3000, 5000];
        let identities = vec![90.0, 99.0, 90.0];
        let q10 = utils::accuracy_to_phred(90.0);
        let q99 = utils::accuracy_to_phred(99.0);
        let mut q_score_counts = vec![0u64; 41];
        let mut q_score_bases = vec![0u128; 41];
        q_score_counts[q10] = 2;
        q_score_bases[q10] = 6000;
        q_score_counts[q99] = 1;
        q_score_bases[q99] = 3000;

        let data = extract_from_bam::Data {
            lengths: Some(lengths),
            num_reads: 3,
            all_counts: 3,
            identities: Some(identities),
            q_score_hist: Some(extract_from_bam::QScoreHistogramData {
                counts: q_score_counts,
                bases: q_score_bases,
            }),
            tids: None,
            starts: None,
            ends: None,
            phasesets: None,
            exons: None,
        };

        let histograms = build_histograms(&data);
        assert_eq!(histograms.read_length.bins.len(), 6);
        assert_eq!(histograms.read_length.bins[0].count, 1);
        assert_eq!(histograms.read_length.bins[0].bases, 1000);
        assert_eq!(histograms.read_length.bins[1].count, 1);
        assert_eq!(histograms.read_length.bins[1].bases, 3000);
        assert_eq!(histograms.read_length.bins[2].count, 1);
        assert_eq!(histograms.read_length.bins[2].bases, 5000);

        let q_score = histograms.q_score.as_ref().expect("Missing Q-score histogram");
        assert_eq!(q_score.bins.len(), 41);
        assert_eq!(q_score.bins[q10].count, 2);
        assert_eq!(q_score.bins[q10].bases, 6000);
        assert_eq!(q_score.bins[q99].count, 1);
        assert_eq!(q_score.bins[q99].bases, 3000);

        let mut metrics_obj = metrics::Metrics::new(metrics::FileInfo {
            name: "test".to_string(),
            path: "test".to_string(),
            creation_time: "now".to_string(),
        });
        metrics_obj.histograms = Some(histograms);
        let json_value = serde_json::to_value(&metrics_obj).expect("Serialize metrics to JSON");
        assert!(json_value.get("histograms").is_some());
        assert!(json_value["histograms"]["read_length"]["bins"].is_array());
        assert!(json_value["histograms"]["q_score"]["bins"].is_array());
    }

    #[test]
    fn histogram_counts_tsv_scaled_uses_bases() {
        let lengths = vec![1000u128, 3000, 5000];
        let mut output = Vec::new();
        output_histogram_counts_tsv(&lengths, &mut output, true);

        let output = String::from_utf8(output).expect("TSV output is valid UTF-8");
        let lines: Vec<&str> = output.trim().lines().collect();
        assert_eq!(lines[0], "bin_start\tbin_end\tbases");
        assert_eq!(lines[1], "0\t2000\t1000");
        assert_eq!(lines[2], "2000\t4000\t3000");
        assert_eq!(lines[3], "4000\t6000\t5000");
    }
}
