pub fn get_n50(lengths: &Vec<u64>, nb_bases_total: u64) -> u64 {
    let mut acc = 0;
    for val in lengths.iter() {
        acc += *val;
        if acc > nb_bases_total / 2 {
            return *val;
        }
    }

    lengths[lengths.len() - 1]
}

pub fn get_n75(lengths: &Vec<u64>, nb_bases_total: u64) -> u64 {
    let mut acc = 0;
    for val in lengths.iter() {
        acc += *val;
        if acc as f64 > nb_bases_total as f64 * 0.75 {
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

pub fn median_phaseblocks(mut array: Vec<f32>) -> f32 {
    array.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());
    if (array.len() % 2) == 0 {
        let ind_left = array.len() / 2 - 1;
        let ind_right = array.len() / 2;
        (array[ind_left] + array[ind_right]) / 2.0
    } else {
        array[array.len() / 2]
    }
}
