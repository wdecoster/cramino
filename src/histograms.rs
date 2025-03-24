use itertools::Itertools;
use std::cmp::max;

// the histograms below are fully defined by the step size and the maximum value
// the step size is the size of each bin
// the maximum value is the value of the last bin
// the last bin (overflow), is printed separately
// the dotsize variable determines how many reads are represented by a single dot
// the functions are duplicated for flexibility in printing the name and formatting the histogram labels
// as well as for future customizations
// in principle it would be possible to enable the user to change the step size or max value, but I don't want to add too many options to the CLI

pub fn make_histogram_lengths(array: &[u128]) {
    // dynamically set the maximum value based on the maximum read length, capped at 60k
    let max_read_length = array.iter().copied().max().expect("Array is empty, cannot find max");
    let max_value = std::cmp::min(
        60_000,
        (((max_read_length as f64) / 10_000.0).ceil() as usize) * 10_000
    );
    let stepsize: u128 = 2000;
    let step_count = max_value / stepsize as usize;
    let mut counts = vec![0; step_count];
    let mut overflow = 0; // Track overflow reads
    
    for &value in array {
        if value >= max_value as u128 {
            overflow += 1; // Count overflow reads directly
        } else {
            let index = (value / stepsize) as usize;
            counts[index] += 1;
        }
    }
    
    // the dotsize variable determines how many reads are represented by a single dot
    let dotsize = max(array.len() / 500, 1);
    
    println!("\n\n# Histogram for read lengths:");
    // print every entry in the vector
    for (index, entry) in counts.iter().dropping_back(1).enumerate() {
        println!(
            "{: >11} {}",
            format!(
                "{}-{}",
                index as u128 * stepsize,
                (index + 1) * stepsize as usize
            ),
            "∎".repeat(entry / dotsize)
        );
    }
    
    // Only print the overflow bin if there are actually overflow reads
    if overflow > 0 {
        println!(
            "{: >11} {}",
            format!("{}+", max_value),
            "∎".repeat(overflow / dotsize)
        );
    }
}

pub fn make_histogram_identities(array: &[f64]) {
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
    println!("\n\n# Histogram for Phred-scaled accuracies:");
    // print every entry in the vector, except the last one which is done separately
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

pub fn make_histogram_phaseblocks(array: &[i64]) {
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
    println!("\n\n# Histogram for phaseblock lengths:");
    // print every entry in the vector, except the last one which is done separately
    for (index, entry) in counts.iter().dropping_back(1).enumerate() {
        println!(
            "{: >14} {}",
            format!(
                "{}-{}",
                index as i64 * stepsize,
                (index + 1) * stepsize as usize
            ),
            "∎".repeat(entry / dotsize)
        );
    }
    println!(
        "{: >14} {}",
        format!("{}+", (counts.len() - 1) * stepsize as usize),
        "∎".repeat(counts.last().unwrap() / dotsize)
    );
}

pub fn make_histogram_exons(array: &[usize]) {
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
    println!("\n\n# Histogram for number of exons:");
    // print the second entry in the vector. The first entry is 0 exons, which is not used (empty)
    // 1 exon is renamed to unspliced
    println!(
        "{: >9} {}",
        format!("unspliced"),
        "∎".repeat(counts[1] / dotsize)
    );
    // print every entry in the vector, except the first two (first one empty, and second one already done) and last one (done later)
    for (index, entry) in counts.iter().skip(2).dropping_back(1).enumerate() {
        println!(
            "{: >9} {}",
            format!("{} exons", (index + 2)),
            "∎".repeat(entry / dotsize)
        );
    }
    println!(
        "{: >9} {}",
        format!("{}+ exons", (counts.len() - 1)),
        "∎".repeat(counts.last().unwrap() / dotsize)
    );
}