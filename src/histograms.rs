use histo_fp::Histogram;
use std::{cmp::Ordering, collections::HashMap};

pub fn make_histogram_lengths(array: &[u64]) {
    new_make_histogram_lengths(array);
    let mut histogram = Histogram::with_buckets(100, Some(2));
    for value in array.iter() {
        histogram.add(*value as f64);
    }
    println!("\n\n# Histogram for lengths\n{}", histogram);
}

// macro_rules! generate_match_arms {
//     ($start:expr, $end:expr, $count:expr) => {{
//         let mut arms = String::new();
//         let mut start = $start;
//         let step = ($end - $start + 1) / $count;
//         for i in 1..=$count {
//             let end = start + step - 1;
//             arms.push_str(&format!(
//                 "{}..={} => *counts.entry({}).or_insert(0) += 1,\n",
//                 start, end, i
//             ));
//             start = end + 1;
//         }
//         arms.push_str(&format!(
//             "_ => *counts.entry({}).or_insert(0) += 1,\n",
//             $count + 1
//         ));
//         arms
//     }};
// }

pub fn new_make_histogram_lengths(array: &[u64]) {
    let mut counts = HashMap::new();
    for value in array.iter() {
        match value {
            // should maybe write a macro to do this?
            1..=2000 => *counts.entry("1-2000").or_insert(0) += 1,
            2001..=4000 => *counts.entry("2001-4000").or_insert(0) += 1,
            4001..=6000 => *counts.entry("4001-6000").or_insert(0) += 1,
            6001..=8000 => *counts.entry("6001-8000").or_insert(0) += 1,
            8001..=10000 => *counts.entry("8001-10000").or_insert(0) += 1,
            10001..=12000 => *counts.entry("10001-12000").or_insert(0) += 1,
            12001..=14000 => *counts.entry("12001-14000").or_insert(0) += 1,
            14001..=16000 => *counts.entry("14001-16000").or_insert(0) += 1,
            16001..=18000 => *counts.entry("16001-18000").or_insert(0) += 1,
            18001..=20000 => *counts.entry("18001-20000").or_insert(0) += 1,
            20001..=22000 => *counts.entry("20001-22000").or_insert(0) += 1,
            22001..=24000 => *counts.entry("22001-24000").or_insert(0) += 1,
            24001..=26000 => *counts.entry("24001-26000").or_insert(0) += 1,
            26001..=28000 => *counts.entry("26001-28000").or_insert(0) += 1,
            28001..=30000 => *counts.entry("28001-30000").or_insert(0) += 1,
            30001..=32000 => *counts.entry("30001-32000").or_insert(0) += 1,
            32001..=34000 => *counts.entry("32001-34000").or_insert(0) += 1,
            34001..=36000 => *counts.entry("34001-36000").or_insert(0) += 1,
            36001..=38000 => *counts.entry("36001-38000").or_insert(0) += 1,
            38001..=40000 => *counts.entry("38001-40000").or_insert(0) += 1,
            40001..=42000 => *counts.entry("40001-42000").or_insert(0) += 1,
            42001..=44000 => *counts.entry("42001-44000").or_insert(0) += 1,
            44001..=46000 => *counts.entry("44001-46000").or_insert(0) += 1,
            46001..=48000 => *counts.entry("46001-48000").or_insert(0) += 1,
            48001..=50000 => *counts.entry("48001-50000").or_insert(0) += 1,
            50001..=52000 => *counts.entry("50001-52000").or_insert(0) += 1,
            52001..=54000 => *counts.entry("52001-54000").or_insert(0) += 1,
            54001..=56000 => *counts.entry("54001-56000").or_insert(0) += 1,
            56001..=58000 => *counts.entry("56001-58000").or_insert(0) += 1,
            58001..=60000 => *counts.entry("58001-60000").or_insert(0) += 1,
            _ => *counts.entry("60001-inf").or_insert(0) += 1,
        }
    }
    // the dotsize variable determines how many reads are represented by a single dot
    // I either have to set this dynamically or experiment with it further
    let dotsize = array.len() / 500;
    // sort the keys as numbers up to the '-', not strings
    let mut keys: Vec<&str> = counts.keys().copied().collect();
    keys.sort_by_key(|k| k.split('-').next().unwrap().parse::<usize>().unwrap());

    for key in keys {
        println!(
            "{: >11} {}",
            key,
            "∎".repeat(counts.get(key).unwrap() / dotsize)
        );
    }
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
