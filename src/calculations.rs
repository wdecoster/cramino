use std::collections::HashMap;

/// N50 (or other percentile) of lengths sorted in descending order
pub fn get_n<T: Copy + Into<u128>>(lengths: &[T], nb_bases_total: u128, percentile: f64) -> u128 {
    // Handle empty array case
    if lengths.is_empty() {
        return 0; // Return 0 for N50/N75 when no reads match the criteria
    }

    let mut acc = 0;
    for val in lengths.iter().map(|v| (*v).into()) {
        acc += val;
        if acc as f64 >= nb_bases_total as f64 * percentile {
            return val;
        }
    }

    lengths[lengths.len() - 1].into()
}

/// Median of a sorted array, 0 if the array is empty
pub fn median<T: Into<f64> + Copy>(array: &[T]) -> f64 {
    if array.is_empty() {
        return 0.0;
    }
    if array.len().is_multiple_of(2) {
        let ind_left = array.len() / 2 - 1;
        let ind_right = array.len() / 2;
        (array[ind_left].into() + array[ind_right].into()) / 2.0
    } else {
        array[array.len() / 2].into()
    }
}

/// Median number of exons per read.
/// The exon counts are collected in the order in which the reads appear in the file,
/// so they have to be sorted here before taking the middle value.
pub fn median_splice(array: &[usize]) -> usize {
    let mut array = array.to_vec();
    array.sort_unstable();
    if array.len().is_multiple_of(2) {
        let ind_left = array.len() / 2 - 1;
        let ind_right = array.len() / 2;
        (array[ind_left] + array[ind_right]) / 2
    } else {
        array[array.len() / 2]
    }
}

pub fn modal_accuracy(array: &[f64]) -> f64 {
    // this doesn't work for f64s, so first I multiply by 10 and then divide by 10 at the end to get the original value again
    // it gets converted to an int, so some resolution is lost, but the floating point differences don't really matter anyway
    let inflate = 10.0;
    let frequencies =
        array
            .iter()
            .map(|x| (x * inflate) as i32)
            .fold(HashMap::new(), |mut freqs, value| {
                *freqs.entry(value).or_insert(0) += 1;
                freqs
            });
    // ties are broken by the highest value, so that the result doesn't depend on the
    // (random) iteration order of the HashMap
    let mode = frequencies
        .into_iter()
        .max_by_key(|&(value, count)| (count, value))
        .map(|(value, _)| value);
    mode.expect("Failed getting the modal accuracy!") as f64 / inflate
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_median() {
        // the arrays are sorted by the caller
        assert_eq!(median(&[4.7, 3.2, 1.5]), 3.2);
        assert_eq!(median(&[7.8, 5.6, 3.4, 1.2]), 4.5);
        assert_eq!(median(&[1u32]), 1.0);
        assert_eq!(median::<u32>(&[]), 0.0);
    }

    #[test]
    fn test_median_splice_unsorted() {
        // exon counts arrive in the order of the reads in the file, so the median
        // has to be independent of that order
        assert_eq!(median_splice(&[1, 9, 2]), 2);
        assert_eq!(median_splice(&[9, 2, 1]), 2);
        assert_eq!(median_splice(&[2, 1, 9]), 2);
    }

    #[test]
    fn test_median_splice_even() {
        // the two middle values of the sorted array are 2 and 4
        assert_eq!(median_splice(&[7, 2, 1, 4]), 3);
    }

    #[test]
    fn test_n50_when_cumulative_sum_hits_exactly_half() {
        // the first read holds exactly half of the bases, and half of the bases
        // are in reads of at least this length, so the N50 is 5 (not 3)
        assert_eq!(get_n(&[5u32, 3, 2], 10, 0.50), 5);
    }

    #[test]
    fn test_modal_accuracy_tie() {
        // 90.0 and 99.9 are both found twice
        for _ in 0..20 {
            assert_eq!(modal_accuracy(&[90.0, 99.9, 90.0, 99.9, 95.0]), 99.9);
        }
    }

    #[test]
    fn test_modal_accuracy() {
        let array = [1.1, 2.2, 2.2, 3.3, 4.4];
        let expected = 2.2;
        let result = modal_accuracy(&array);
        assert_eq!(result, expected, "The modal accuracy calculation failed!");
    }
}
