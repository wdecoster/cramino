pub fn get_n50(lengths: &Vec<u32>, nb_bases_total: u32) -> u32 {
    let mut acc = 0;
    for val in lengths.iter() {
        acc += *val;
        if acc > nb_bases_total / 2 {
            return *val;
        }
    }

    lengths[lengths.len() - 1]
}

/// Take the median of the lengths of the STRs, relative to the reference genome
/// If the vector is empty then return NAN
pub fn median(array: &Vec<u32>) -> f32 {
    if (array.len() % 2) == 0 {
        let ind_left = array.len() / 2 - 1;
        let ind_right = array.len() / 2;
        (array[ind_left] + array[ind_right]) as f32 / 2.0
    } else {
        array[(array.len() / 2)] as f32
    }
}
