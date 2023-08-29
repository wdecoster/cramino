use histo_fp::Histogram;
use itertools::Itertools;
use std::cmp::Ordering;

pub fn make_histogram_lengths(array: &[u64]) {
    let mut counts = vec![0; 31];
    let stepsize = 2000;
    for value in array.iter() {
        let index = (value / stepsize) as usize;
        if index < counts.len() {
            counts[index] += 1;
        } else {
            let last_index = counts.len() - 1;
            counts[last_index] += 1;
        }
    }
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
    let bins = 100.0
        - array
            .iter()
            .min_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Less))
            .unwrap();
    let mut histogram = Histogram::with_buckets(bins as u64, None);
    for value in array.iter() {
        histogram.add(*value);
    }
    println!("\n\n# Histogram for identities\n{}", histogram);
}

pub fn make_histogram_phaseblocks(array: &[i64]) {
    let mut histogram = Histogram::with_buckets(100, Some(2));
    for value in array.iter() {
        histogram.add(*value as f64);
    }
    println!("\n\n# Histogram for phaseblocks\n{}", histogram);
}
