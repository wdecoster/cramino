pub fn splice_metrics(mut exon_counts: Vec<usize>) {
    exon_counts.sort_unstable();
    let num_reads = exon_counts.len();
    println!("Median number of exons\t{:.2}", median(&exon_counts));
    println!(
        "Mean number of exons\t{:.2}",
        exon_counts.iter().sum::<usize>() / num_reads
    )
}

pub fn median(array: &Vec<usize>) -> f32 {
    if (array.len() % 2) == 0 {
        let ind_left = array.len() / 2 - 1;
        let ind_right = array.len() / 2;
        (array[ind_left] + array[ind_right]) as f32 / 2.0
    } else {
        array[(array.len() / 2)] as f32
    }
}
