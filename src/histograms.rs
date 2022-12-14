use histo_fp::Histogram;
use std::cmp::Ordering;

pub fn make_histogram_lengths(array: &[u64]) {
    let mut histogram = Histogram::with_buckets(100, Some(2));
    for value in array.iter() {
        histogram.add(*value as f64);
    }
    println!("\n\n# Histogram for lengths\n{}", histogram);
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
