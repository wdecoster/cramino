use histo_fp::Histogram;
use itertools::Itertools;

pub fn make_histogram_lengths(array: &[u64]) {
    let stepsize: u64 = 2000;
    let max_value = 60_000;
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
    let dotsize = array.len() / 500;
    // sort the keys as numbers up to the '-', not strings
    println!("\n\n# Histogram for lengths:");
    for (index, entry) in counts.iter().dropping_back(1).enumerate() {
        println!(
            "{: >11} {}",
            format!(
                "{}-{}",
                index as u64 * stepsize,
                (index + 1) * stepsize as usize
            ),
            "∎".repeat(entry / dotsize)
        );
    }
    println!(
        "{: >11} {}",
        format!("{}+", (counts.len() - 1) * stepsize as usize),
        "∎".repeat(counts.last().unwrap() / dotsize)
    );
}

pub fn make_histogram_identities(array: &[f64]) {
    let stepsize: u64 = 1;
    let max_value = 40;
    let step_count = max_value / stepsize as usize;
    let mut counts = vec![0; step_count + 1];
    for value in array.iter().map(|x| accuracy_to_phred(*x)) {
        if value < max_value {
            counts[value] += 1;
        }
    }
    // the last bin is for all values above the last step
    counts[step_count] = array.len() - counts.iter().sum::<usize>();
    // the dotsize variable determines how many reads are represented by a single dot
    // I either have to set this dynamically or experiment with it further
    let dotsize = array.len() / 500;
    // sort the keys as numbers up to the '-', not strings
    println!("\n\n# Histogram for Phred-scaled accuracies:");
    for (index, entry) in counts.iter().dropping_back(1).enumerate() {
        println!(
            "{: >6} {}",
            format!(
                "Q{}-{}",
                index as u64 * stepsize,
                (index + 1) * stepsize as usize
            ),
            "∎".repeat(entry / dotsize)
        );
    }
    println!(
        "{: >6} {}",
        format!("Q{}+", (counts.len() - 1) * stepsize as usize),
        "∎".repeat(counts.last().unwrap() / dotsize)
    );
}

fn accuracy_to_phred(identity: f64) -> usize {
    // convert identity to phred scale
    // but return as usize (as that will be used for the histogram)
    // this is therefore not accurate for other applications
    (-10.0 * (1.0 - identity).log10()) as usize
}

pub fn make_histogram_phaseblocks(array: &[i64]) {
    let mut histogram = Histogram::with_buckets(100, Some(2));
    for value in array.iter() {
        histogram.add(*value as f64);
    }
    println!("\n\n# Histogram for phaseblocks\n{}", histogram);
}
