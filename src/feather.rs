use arrow::datatypes::{DataType, Field, Schema};
use std::fs::File;
use std::sync::Arc;

use arrow::{
    self,
    array::{ArrayRef, Float64Array, UInt64Array},
    ipc::writer::FileWriter,
    record_batch::RecordBatch,
};

/// Saves the length and identity of each alignment, in the order of the input.
/// With --ubam, the identities are estimated from the base qualities.
/// Reads without an identity (NaN) are written as null
pub fn save_as_arrow(filename: String, lengths: Vec<u64>, identities: &[f64]) {
    let identities: Vec<Option<f64>> = identities
        .iter()
        .map(|identity| (!identity.is_nan()).then_some(*identity))
        .collect();
    let schema = Arc::new(Schema::new(vec![
        Field::new("identities", DataType::Float64, true),
        Field::new("lengths", DataType::UInt64, false),
    ]));
    let columns: Vec<ArrayRef> = vec![
        Arc::new(Float64Array::from(identities)),
        Arc::new(UInt64Array::from(lengths)),
    ];
    let batch = RecordBatch::try_new(schema.clone(), columns).expect("create arrow batch error");
    let buffer = File::create(filename).expect("create file error");

    let mut writer = FileWriter::try_new(buffer, &schema).expect("create arrow file writer error");

    writer.write(&batch).expect("write arrow batch error");
    writer.finish().expect("finish write arrow error");
}
