pub fn splice_metrics(exon_counts: &Vec<usize>) {
    let num_reads = exon_counts.len();
    println!(
        "Median number of exons\t{:.2}",
        crate::calculations::median_splice(exon_counts)
    );
    println!(
        "Mean number of exons\t{:.2}",
        (exon_counts.iter().sum::<usize>() as f32) / (num_reads as f32)
    );
    // number of reads with 1 exon (unspliced)
    let num_single_exon = exon_counts.iter().filter(|&&x| x == 1).count();
    println!(
        "Fraction unspliced reads\t{:.2}",
        (num_single_exon as f32) / (num_reads as f32)
    );
}
