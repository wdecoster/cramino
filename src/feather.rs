use arrow::datatypes::{DataType, Field, Schema};
use std::fs::File;
use std::sync::Arc;

use arrow::{
    self,
    array::{Float64Array, UInt64Array},
    ipc::writer::FileWriter,
    record_batch::RecordBatch,
};

pub fn save_as_arrow(filename: String, lengths: Vec<u64>, identities: Vec<f64>) {
    let identities_array = Arc::new(Float64Array::from(identities)) as _;
    let lengths_array = Arc::new(UInt64Array::from(lengths)) as _;
    let batch =
        RecordBatch::try_from_iter([("identities", identities_array), ("lengths", lengths_array)])
            .unwrap();

    let schema = Schema::new(vec![
        Field::new("identities", DataType::Float64, false),
        Field::new("lengths", DataType::UInt64, false),
    ]);
    let buffer = File::create(filename).expect("create file error");

    let mut writer = FileWriter::try_new(buffer, &schema).expect("create arrow file writer error");

    writer.write(&batch).expect("write arrow batch error");
    writer.finish().expect("finish write arrow error");
}

pub fn save_as_arrow_ubam(filename: String, lengths: Vec<u64>) {
    let lengths_array = Arc::new(UInt64Array::from(lengths)) as _;
    let batch = RecordBatch::try_from_iter([("lengths", lengths_array)]).unwrap();

    let schema = Schema::new(vec![Field::new("lengths", DataType::UInt64, false)]);
    let buffer = File::create(filename).expect("create file error");

    let mut writer = FileWriter::try_new(buffer, &schema).expect("create arrow file writer error");

    writer.write(&batch).expect("write arrow batch error");
    writer.finish().expect("finish write arrow error");
}
