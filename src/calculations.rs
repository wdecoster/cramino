pub fn get_n(lengths: &Vec<u64>, nb_bases_total: u64, percentile: f64) -> u64 {
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
    if (array.len() % 2) == 0 {
        let ind_left = array.len() / 2 - 1;
        let ind_right = array.len() / 2;
        (array[ind_left].into() + array[ind_right].into()) / 2.0
    } else {
        array[array.len() / 2].into()
    }
}

pub fn median_length(array: &[u64]) -> f64 {
    if (array.len() % 2) == 0 {
        let ind_left = array.len() / 2 - 1;
        let ind_right = array.len() / 2;

        (array[ind_left] + array[ind_right]) as f64 / 2.0
    } else {
        array[array.len() / 2] as f64
    }
}

/// Returns the median of an array of normalized read counts.
///
/// The array is assumed to be a slice of normalized read counts for each
/// chromosome after having been aligned using minimap2.
///
/// # Examples
///
/// ```rust, ignore
/// # use crate::calculations::median_phaseblocks;
/// // Array with odd number of elements
/// let v1 = vec![3.2, 1.5, 4.7];
/// assert_eq!(median_phaseblocks(v1), 3.2);
///
/// // Array with even number of elements
/// let v2 = vec![1.2, 3.4, 5.6, 7.8];
/// assert_eq!(median_phaseblocks(v2), 4.5);
///
/// // Array with a single element
/// let v3 = vec![1.0];
/// assert_eq!(median_phaseblocks(v3), 1.0);
/// ```
pub fn median_phaseblocks(mut array: Vec<f32>) -> f32 {
    array.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());
    if (array.len() % 2) == 0 {
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

pub fn median_splice(array: &Vec<usize>) -> usize {
    if (array.len() % 2) == 0 {
        let ind_left = array.len() / 2 - 1;
        let ind_right = array.len() / 2;
        (array[ind_left] + array[ind_right]) / 2
    } else {
        array[array.len() / 2]
    }
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
}
