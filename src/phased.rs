use itertools::izip;
use log::error;

/// Length of each phaseblock, from the leftmost start to the rightmost end of the alignments
/// (primary and supplementary) with the same phaseset on the same chromosome. Reads of
/// different phasesets can overlap, so the reads are grouped per phaseset rather than split
/// at each change of phaseset along the genome.
pub fn phase_metrics(
    tids: &[i32],
    starts: Vec<i64>,
    ends: Vec<i64>,
    phasesets: &[Option<u32>],
) -> Vec<i64> {
    let mut phased_reads = izip!(tids, phasesets, starts, ends)
        .filter_map(|(tid, phaseset, start, end)| phaseset.map(|p| (*tid, p, start, end)))
        .collect::<Vec<_>>();
    if phased_reads.is_empty() {
        error!("Not a single phased read found!");
        return vec![];
    }
    phased_reads.sort_unstable();

    let mut phaseblocks = vec![];
    let mut reads = phased_reads.into_iter();
    let (mut block_tid, mut block_phaseset, mut block_start, mut block_end) = reads.next().unwrap();
    for (tid, phaseset, start, end) in reads {
        if tid == block_tid && phaseset == block_phaseset {
            // reads are sorted by start, but a later read can end before an earlier one
            block_end = block_end.max(end);
        } else {
            phaseblocks.push(block_end - block_start);
            (block_tid, block_phaseset, block_start, block_end) = (tid, phaseset, start, end);
        }
    }
    phaseblocks.push(block_end - block_start);
    phaseblocks
}

/// Median phaseblock length, which requires a sorted array
pub fn median(array: &[i64]) -> f64 {
    if array.len().is_multiple_of(2) {
        let ind_left = array.len() / 2 - 1;
        let ind_right = array.len() / 2;
        (array[ind_left] + array[ind_right]) as f64 / 2.0
    } else {
        array[array.len() / 2] as f64
    }
}

/// N50 of the phaseblock lengths, which requires the lengths sorted in descending order
pub fn get_n50(lengths: &[i64], nb_bases_total: i64) -> i64 {
    let mut acc = 0;
    for val in lengths.iter() {
        acc += *val;
        if 2 * acc >= nb_bases_total {
            return *val;
        }
    }

    lengths[lengths.len() - 1]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_median_of_sorted_phaseblocks() {
        // phaseblocks have to be sorted by the caller, as they are collected in genomic order
        let mut phaseblocks = vec![100, 5000, 300];
        phaseblocks.sort_unstable_by(|a, b| b.cmp(a));
        assert_eq!(median(&phaseblocks), 300.0);
    }

    #[test]
    fn test_n50_of_sorted_phaseblocks() {
        // sorted descending, the accumulated length passes half of 5400 at 5000
        let mut phaseblocks = vec![100, 5000, 300];
        let total = phaseblocks.iter().sum::<i64>();
        phaseblocks.sort_unstable_by(|a, b| b.cmp(a));
        assert_eq!(get_n50(&phaseblocks, total), 5000);
    }

    #[test]
    fn test_n50_when_cumulative_sum_hits_exactly_half() {
        assert_eq!(get_n50(&[5, 3, 2], 10), 5);
    }

    #[test]
    fn test_phaseblock_extends_to_rightmost_end() {
        // the second read starts later but ends earlier than the first
        let blocks = phase_metrics(&[0, 0], vec![0, 100], vec![10000, 200], &[Some(1), Some(1)]);
        assert_eq!(blocks, vec![10000]);
    }

    #[test]
    fn test_overlapping_phasesets_are_not_split() {
        // reads of phaseset 1 and 2 alternate along the genome
        let blocks = phase_metrics(
            &[0, 0, 0, 0],
            vec![0, 1500, 2000, 2500],
            vec![1000, 2500, 3000, 3500],
            &[Some(1), Some(2), Some(1), Some(2)],
        );
        assert_eq!(blocks, vec![3000, 2000]);
    }

    #[test]
    fn test_same_phaseset_on_other_chromosome_is_other_block() {
        let blocks = phase_metrics(
            &[0, 1, 1],
            vec![0, 0, 50],
            vec![100, 100, 400],
            &[Some(1), Some(1), None],
        );
        assert_eq!(blocks, vec![100, 100]);
    }
}
