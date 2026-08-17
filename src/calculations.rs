use std::collections::HashMap;

pub fn get_n(lengths: &[u128], nb_bases_total: u128, percentile: f64) -> u128 {
    // Handle empty array case
    if lengths.is_empty() {
        return 0; // Return 0 for N50/N75 when no reads match the criteria
    }

    let mut acc = 0;
    for val in lengths.iter() {
        acc += *val;
        if acc as f64 > nb_bases_total as f64 * percentile {
            return *val;
        }
    }

    lengths[lengths.len() - 1]
}

pub fn median<T: Into<f64> + Copy>(array: &[T]) -> f64 {
    if array.len().is_multiple_of(2) {
        let ind_left = array.len() / 2 - 1;
        let ind_right = array.len() / 2;
        (array[ind_left].into() + array[ind_right].into()) / 2.0
    } else {
        array[array.len() / 2].into()
    }
}

pub fn median_length(array: &[u128]) -> f64 {
    if array.len().is_multiple_of(2) {
        let ind_left = array.len() / 2 - 1;
        let ind_right = array.len() / 2;

        (array[ind_left] + array[ind_right]) as f64 / 2.0
    } else {
        array[array.len() / 2] as f64
    }
}

pub fn median_phaseblocks(mut array: Vec<f32>) -> f32 {
    array.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());
    if array.len().is_multiple_of(2) {
        let ind_left = array.len().checked_div(2).unwrap().saturating_sub(1);
        let ind_right = array.len().checked_div(2).unwrap_or(0);
        if (ind_left == 0) & (ind_right == 0) {
            return 0.0;
        }
        (array[ind_left] + array[ind_right]) / 2.0
    } else {
        array[array.len() / 2]
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
    let mode = frequencies
        .into_iter()
        .max_by_key(|&(_, count)| count)
        .map(|(value, _)| value);
    mode.expect("Failed getting the modal accuracy!") as f64 / inflate
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_median_odd() {
        let v1 = vec![3.2, 1.5, 4.7];
        assert_eq!(median_phaseblocks(v1), 3.2);
    }

    #[test]
    fn test_median_even() {
        let v2 = vec![1.2, 3.4, 5.6, 7.8];
        assert_eq!(median_phaseblocks(v2), 4.5);
    }

    #[test]
    fn test_median_single_element() {
        let v3 = vec![1.0];
        assert_eq!(median_phaseblocks(v3), 1.0);
    }

    #[test]
    fn test_median_no_element() {
        let v3 = vec![];
        assert_eq!(median_phaseblocks(v3), 0.0);
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
    fn test_modal_accuracy() {
        let array = [1.1, 2.2, 2.2, 3.3, 4.4];
        let expected = 2.2;
        let result = modal_accuracy(&array);
        assert_eq!(result, expected, "The modal accuracy calculation failed!");
    }
}
