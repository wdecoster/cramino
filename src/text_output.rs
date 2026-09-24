use crate::metrics;

pub fn print_text_output(metrics: &metrics::Metrics) {
    // Print file info
    println!("File name\t{}", metrics.file_info.name);

    // Print alignment stats
    println!(
        "Number of alignments\t{}",
        metrics.alignment_stats.num_alignments
    );
    println!(
        "% from total reads\t{:.2}",
        metrics.alignment_stats.percent_from_total
    );
    println!("Number of reads\t{}", metrics.alignment_stats.num_reads);

    // Print read stats
    println!("Yield [Gb]\t{:.2}", metrics.read_stats.yield_gb);
    match metrics.read_stats.mean_coverage {
        Some(coverage) => println!("Mean coverage\t{:.2}", coverage),
        None => println!("Mean coverage\tNA"),
    }
    println!(
        "Yield [Gb] (>25kb)\t{:.2}",
        metrics.read_stats.yield_gb_long
    );
    println!("N50\t{}", metrics.read_stats.n50);
    println!("N75\t{}", metrics.read_stats.n75);
    println!("Median length\t{:.2}", metrics.read_stats.median_length);
    println!("Mean length\t{:.2}", metrics.read_stats.mean_length);
    println!("N50 aligned\t{}", metrics.read_stats.n50_aligned);
    println!("N75 aligned\t{}", metrics.read_stats.n75_aligned);
    println!(
        "Median length aligned\t{:.2}",
        metrics.read_stats.median_length_aligned
    );
    println!(
        "Mean length aligned\t{:.2}",
        metrics.read_stats.mean_length_aligned
    );
    println!();

    // Print identity stats if available
    if let Some(identity_stats) = &metrics.identity_stats {
        println!("Median identity\t{:.2}", identity_stats.median_identity);
        println!("Mean identity\t{:.2}", identity_stats.mean_identity);
        println!("Modal identity\t{:.1}", identity_stats.modal_identity);
        println!();
    }
    if let Some(identity_stats) = &metrics.estimated_identity_stats {
        println!(
            "Median estimated identity\t{:.2}",
            identity_stats.median_estimated_identity
        );
        println!(
            "Mean estimated identity\t{:.2}",
            identity_stats.mean_estimated_identity
        );
        println!(
            "Modal estimated identity\t{:.1}",
            identity_stats.modal_estimated_identity
        );
        println!();
    }

    // Print phase stats if available
    if let Some(phase_stats) = &metrics.phase_stats {
        println!("Fraction reads phased\t{:.2}", phase_stats.fraction_phased);
        println!("Number of phaseblocks\t{}", phase_stats.num_phaseblocks);
        println!(
            "Total bases phased [Gb]\t{:.2}",
            phase_stats.total_bases_phased_gb
        );
        println!(
            "Median phaseblock length\t{:.2}",
            phase_stats.median_phaseblock_length
        );
        println!(
            "N50 phaseblock length\t{}",
            phase_stats.n50_phaseblock_length
        );
        println!();
    }

    // Print karyotype stats if available
    if let Some(karyotype_stats) = &metrics.karyotype_stats {
        if !karyotype_stats.is_empty() {
            println!("\n\n# Normalized read count per chromosome\n");
            for chrom_data in karyotype_stats {
                println!(
                    "{}\t{:.2}",
                    chrom_data.chromosome, chrom_data.normalized_count
                );
            }
        } else {
            println!("\n\n# Warning - no contigs found in BAM file!\n");
        }
        println!();
    }

    // Print splice stats if available
    if let Some(splice_stats) = &metrics.splice_stats {
        println!("Median number of exons\t{}", splice_stats.median_exons);
        println!("Mean number of exons\t{:.2}", splice_stats.mean_exons);
        println!(
            "Fraction unspliced reads\t{:.2}",
            splice_stats.fraction_unspliced
        );
        println!();
    }
    // Print file info
    println!("Path\t{}", metrics.file_info.path);
    println!("Creation time\t{}", metrics.file_info.creation_time);
}
