use rust_htslib::bam::{self, Read};
use std::collections::HashMap;

pub fn make_karyotype(tids: &Vec<i32>, bamp: String) {
    let bam = bam::Reader::from_path(bamp).unwrap();
    let header = bam::Header::from_template(bam.header());
    let head_view = bam::HeaderView::from_header(&header);

    let mut tidcount: HashMap<i32, usize> = HashMap::new();
    for tid in tids {
        *tidcount.entry(*tid).or_default() += 1;
    }
    let mut chroms = vec![];
    let mut norm_count = vec![];
    for (tid, count) in tidcount.into_iter() {
        let chrom = std::str::from_utf8(head_view.tid2name(tid.try_into().unwrap())).unwrap();
        let chrom_length = head_view.target_len(tid.try_into().unwrap()).unwrap();
        chroms.push(chrom);
        norm_count.push((count as f32) / (chrom_length as f32));
    }
    let median_count = crate::calculations::median_phaseblocks(norm_count.clone());
    let mut zipped = chroms.iter().zip(norm_count).collect::<Vec<_>>();
    zipped.sort_by_key(|&(&val, _)| val);
    if !zipped.is_empty() {
        println!("\n\n# Normalized read count per chromosome\n");
        for (chrom, count) in zipped {
            println!("{}\t{:.2}", chrom, count / median_count)
        }
    } else {
        println!("\n\n# Warning - no contigs found in BAM file!\n");
    }
}
