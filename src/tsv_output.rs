use crate::metrics;

pub fn print_tsv_output(metrics: &metrics::Metrics) {
    // Prepare headers and values separately
    let mut headers = Vec::new();
    let mut values = Vec::new();
    
    // File info
    headers.push("file_name");
    values.push(metrics.file_info.name.clone());
    headers.push("file_path");
    values.push(metrics.file_info.path.clone());
    headers.push("creation_time");
    values.push(metrics.file_info.creation_time.clone());
    
    // Alignment stats
    headers.push("num_alignments");
    values.push(metrics.alignment_stats.num_alignments.to_string());
    headers.push("percent_from_total");
    values.push(format!("{:.2}", metrics.alignment_stats.percent_from_total));
    headers.push("num_reads");
    values.push(metrics.alignment_stats.num_reads.to_string());
    
    // Read stats
    headers.push("yield_gb");
    values.push(format!("{:.2}", metrics.read_stats.yield_gb));
    headers.push("mean_coverage");
    values.push(format!("{:.2}", metrics.read_stats.mean_coverage));
    headers.push("yield_gb_long");
    values.push(format!("{:.2}", metrics.read_stats.yield_gb_long));
    headers.push("n50");
    values.push(metrics.read_stats.n50.to_string());
    headers.push("n75");
    values.push(metrics.read_stats.n75.to_string());
    headers.push("median_length");
    values.push(format!("{:.2}", metrics.read_stats.median_length));
    headers.push("mean_length");
    values.push(format!("{:.2}", metrics.read_stats.mean_length));
    
    // Identity stats (if available)
    if let Some(identity_stats) = &metrics.identity_stats {
        headers.push("median_identity");
        values.push(format!("{:.2}", identity_stats.median_identity));
        headers.push("mean_identity");
        values.push(format!("{:.2}", identity_stats.mean_identity));
        headers.push("modal_identity");
        values.push(format!("{:.1}", identity_stats.modal_identity));
    }
    
    // Phase stats (if available)
    if let Some(phase_stats) = &metrics.phase_stats {
        headers.push("fraction_phased");
        values.push(format!("{:.2}", phase_stats.fraction_phased));
        headers.push("num_phaseblocks");
        values.push(phase_stats.num_phaseblocks.to_string());
        headers.push("total_bases_phased_gb");
        values.push(format!("{:.2}", phase_stats.total_bases_phased_gb));
        headers.push("median_phaseblock_length");
        values.push(format!("{:.2}", phase_stats.median_phaseblock_length));
        headers.push("n50_phaseblock_length");
        values.push(phase_stats.n50_phaseblock_length.to_string());
    }
    
    // Splice stats (if available)
    if let Some(splice_stats) = &metrics.splice_stats {
        headers.push("median_exons");
        values.push(splice_stats.median_exons.to_string());
        headers.push("mean_exons");
        values.push(format!("{:.2}", splice_stats.mean_exons));
        headers.push("fraction_unspliced");
        values.push(format!("{:.2}", splice_stats.fraction_unspliced));
    }
    
    // Print headers and values as TSV
    println!("{}", headers.join("\t"));
    println!("{}", values.join("\t"));
}