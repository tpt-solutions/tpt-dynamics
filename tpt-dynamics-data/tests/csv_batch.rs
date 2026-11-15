//! `CsvReader::next_batch` must give exactly the records the iterator gives, for
//! any batch size, including quoted fields, blank lines and error positions.

use std::io::Cursor;

use tpt_dynamics_data::csv::{CsvOptions, CsvReader};

fn sample(rows: usize) -> String {
    let mut s = String::from("id,name,value\n");
    for i in 0..rows {
        if i % 17 == 5 {
            s.push('\n'); // blank line, skipped
        }
        if i % 7 == 3 {
            s.push_str(&format!("{i},\"quoted, {i}\",{}.5\n", i * 2));
        } else {
            s.push_str(&format!("{i},plain,{}\n", i * 2));
        }
    }
    s
}

fn collect_iter(text: &str) -> Vec<(u64, Vec<String>)> {
    CsvReader::new(Cursor::new(text.as_bytes().to_vec()), CsvOptions::default())
        .unwrap()
        .map(|r| {
            let r = r.unwrap();
            (r.line(), r.fields().to_vec())
        })
        .collect()
}

fn collect_batches(text: &str, batch: usize) -> Vec<(u64, Vec<String>)> {
    let mut reader =
        CsvReader::new(Cursor::new(text.as_bytes().to_vec()), CsvOptions::default()).unwrap();
    let mut out = Vec::new();
    loop {
        let records = reader.next_batch(batch).unwrap();
        if records.is_empty() {
            break;
        }
        assert!(records.len() <= batch);
        out.extend(records.into_iter().map(|r| (r.line(), r.fields().to_vec())));
    }
    out
}

#[test]
fn batches_match_the_iterator_for_every_batch_size() {
    let text = sample(500);
    let expected = collect_iter(&text);
    assert_eq!(expected.len(), 500);
    for batch in [1, 2, 3, 7, 64, 499, 500, 501, 1 << 20] {
        assert_eq!(collect_batches(&text, batch), expected, "batch size {batch}");
    }
}

#[test]
fn field_count_error_reports_the_same_line() {
    let text = "a,b\n1,2\n3,4\n5\n6,7\n";
    let mut reader =
        CsvReader::new(Cursor::new(text.as_bytes().to_vec()), CsvOptions::default()).unwrap();
    let err = reader.next_batch(100).unwrap_err().to_string();
    let mut iter =
        CsvReader::new(Cursor::new(text.as_bytes().to_vec()), CsvOptions::default()).unwrap();
    let iter_err = loop {
        match iter.next() {
            Some(Err(e)) => break e.to_string(),
            Some(Ok(_)) => continue,
            None => panic!("expected an error"),
        }
    };
    assert_eq!(err, iter_err);
}

#[test]
fn empty_input_gives_empty_batches() {
    let mut reader = CsvReader::new(
        Cursor::new(Vec::new()),
        CsvOptions { has_header: false, ..CsvOptions::default() },
    )
    .unwrap();
    assert!(reader.next_batch(10).unwrap().is_empty());
    assert!(reader.next_batch(10).unwrap().is_empty());
}
