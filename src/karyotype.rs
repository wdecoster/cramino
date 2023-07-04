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
        if !chrom.contains('_') {
            let chrom_length = head_view.target_len(tid.try_into().unwrap()).unwrap();
            chroms.push(chrom);
            norm_count.push((count as f32) / (chrom_length as f32));
        }
    }
    let median_count = median(norm_count.clone());
    let mut zipped = chroms.iter().zip(norm_count).collect::<Vec<_>>();
    zipped.sort_by_key(|&(&val, _)| val);
    println!("\n\n# Normalized read count per chromosome\n");
    for (chrom, count) in zipped {
        println!("{}\t{:.2}", chrom, count / median_count)
    }
}

pub fn median(mut array: Vec<f32>) -> f32 {
    array.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());
    if (array.len() % 2) == 0 {
        let ind_left = array.len() / 2 - 1;
        let ind_right = array.len() / 2;
        (array[ind_left] + array[ind_right]) / 2.0
    } else {
        array[array.len() / 2]
    }
}
