use itertools::izip;
use log::error;

pub fn phase_metrics(
    tids: &[i32],
    starts: Vec<i64>,
    ends: Vec<i64>,
    phasesets: &Vec<Option<u32>>,
) -> Vec<i64> {
    let mut phased_reads = izip!(tids, starts, ends, phasesets)
        .filter(|(_, _, _, p)| p.is_some())
        .collect::<Vec<_>>();
    phased_reads.sort_unstable();

    let num_phased_reads = phased_reads.len();
    if num_phased_reads == 0 {
        error!("Not a single phased read found!");
        return vec![];
    }

    let mut phased_reads_iter = phased_reads.into_iter();
    let (mut chrom1, mut start1, mut block_end, mut phaseset1) = phased_reads_iter.next().unwrap();
    let mut phaseblocks = vec![];
    for (chrom, start, end, phaseset) in phased_reads_iter {
        if chrom == chrom1 && phaseset == phaseset1 {
            block_end = end;
            continue;
        } else {
            phaseblocks.push(block_end - start1);
            chrom1 = chrom;
            start1 = start;
            block_end = end;
            phaseset1 = phaseset;
        }
    }
    phaseblocks.push(block_end - start1);
    phaseblocks
}

pub fn median(array: &[i64]) -> f64 {
    if array.len().is_multiple_of(2) {
        let ind_left = array.len() / 2 - 1;
        let ind_right = array.len() / 2;
        (array[ind_left] + array[ind_right]) as f64 / 2.0
    } else {
        array[array.len() / 2] as f64
    }
}

pub fn get_n50(lengths: &[i64], nb_bases_total: i64) -> i64 {
    let mut acc = 0;
    for val in lengths.iter() {
        acc += *val;
        if acc > nb_bases_total / 2 {
            return *val;
        }
    }

    lengths[lengths.len() - 1]
}
