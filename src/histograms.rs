use itertools::Itertools;
use std::cmp::max;
use std::fs::File;
use std::io::{self, Write};

use crate::extract_from_bam;

fn output_histogram_counts_tsv<W: Write>(array: &[u128], writer: &mut W) {
    // Handle empty array case
    if array.is_empty() {
        return;
    }

    // dynamically set the maximum value based on the maximum read length, capped at 60k
    let max_read_length = array
        .iter()
        .copied()
        .max()
        .expect("Array is empty, cannot find max");
    let max_value = std::cmp::min(
        60_000,
        (((max_read_length as f64) / 10_000.0).ceil() as usize) * 10_000,
    );
    let stepsize: u128 = 2000;
    let step_count = (max_value as u128 / stepsize) as usize;
    let mut counts = vec![0; step_count];
    let mut overflow = 0; // Track overflow reads

    for &value in array {
        if value >= max_value as u128 {
            overflow += 1;
        } else {
            let index = (value / stepsize) as usize;
            counts[index] += 1;
        }
    }

    // Write TSV header with leading newline for formatting
    writeln!(writer, "\nbin_start\tbin_end\tcount")
        .expect("Unable to write histogram counts header");

    // Write each bin
    for (index, count) in counts.iter().enumerate() {
        writeln!(
            writer,
            "{}\t{}\t{}",
            index as u128 * stepsize,
            (index + 1) as u128 * stepsize,
            count
        )
        .expect("Unable to write histogram counts");
    }

    // Write overflow bin if it has any reads
    if overflow > 0 {
        writeln!(writer, "{}+\tNA\t{}", max_value, overflow).expect("Unable to write overflow bin");
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
    // dynamically set the maximum value based on the maximum read length, capped at 60k
    let max_read_length = array
        .iter()
        .copied()
        .max()
        .expect("Array is empty, cannot find max");
    let max_value = std::cmp::min(
        60_000,
        (((max_read_length as f64) / 10_000.0).ceil() as usize) * 10_000,
    );
    let stepsize: u128 = 2000;
    let step_count = max_value / stepsize as usize;
    let mut counts = vec![0; step_count];
    let mut basepairs = vec![0u128; step_count];
    let mut overflow = 0; // Track overflow reads
    let mut overflow_bp = 0u128;
    for &value in array {
        if value >= max_value as u128 {
            overflow += 1; // Count overflow reads directly
            overflow_bp += value;
        } else {
            let index = (value / stepsize) as usize;
            counts[index] += 1;
            basepairs[index] += value;
        }
    }
    let dotsize = if scaled {
        let total_bp: u128 = basepairs.iter().sum::<u128>() + overflow_bp;
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
    for (index, (entry, bp)) in counts.iter().zip(basepairs.iter()).enumerate() {
        let bar = if scaled {
            "∎".repeat(((*bp as usize) / dotsize).max(0))
        } else {
            "∎".repeat((*entry / dotsize).max(0))
        };
        writeln!(
            writer,
            "{: >11} {}",
            format!(
                "{}-{}",
                index as u128 * stepsize,
                (index + 1) * stepsize as usize
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
        writeln!(writer, "{: >11} {}", format!("{}+", max_value), bar)
            .expect("Unable to write histogram");
    }
}

fn make_histogram_identities<W: Write>(array: &[f64], writer: &mut W) {
    let stepsize: u64 = 1;
    let max_value = 40;
    let step_count = max_value / stepsize as usize;
    let mut counts = vec![0; step_count + 1];
    for value in array.iter().map(|x| crate::utils::accuracy_to_phred(*x)) {
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
) -> Result<(), Box<dyn std::error::Error>> {
    let mut writer: Box<dyn Write> = if let Some(file) = hist_count_file {
        Box::new(File::create(file)?)
    } else {
        Box::new(io::stdout())
    };
    if let Some(lengths) = &metrics_data.lengths {
        output_histogram_counts_tsv(lengths, &mut writer);
    }
    Ok(())
}
