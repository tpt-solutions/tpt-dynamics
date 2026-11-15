# Regenerates the Parquet fixtures used by tests/parquet_import.rs.
# Usage: python generate.py <output-dir>   (requires pyarrow)
import pyarrow as pa, pyarrow.parquet as pq, os, sys
out = sys.argv[1]
n = 1000
base = {
    "i64": pa.array([i * 7 - 300 for i in range(n)], pa.int64()),
    "i32": pa.array([i - 500 for i in range(n)], pa.int32()),
    "f64": pa.array([i * 0.25 - 10.0 for i in range(n)], pa.float64()),
    "f32": pa.array([i * 0.5 for i in range(n)], pa.float32()),
    "flag": pa.array([i % 3 == 0 for i in range(n)], pa.bool_()),
    "cat": pa.array([i % 5 for i in range(n)], pa.int64()),  # low cardinality -> dictionary
}
t = pa.table(base)
for codec in ["NONE", "SNAPPY", "GZIP"]:
    for dict_on in [True, False]:
        name = f"flat_{codec.lower()}_{'dict' if dict_on else 'plain'}.parquet"
        pq.write_table(t, os.path.join(out, name), compression=codec, use_dictionary=dict_on, data_page_version="1.0", row_group_size=300)
# optional column with nulls
nullable = pa.table({
    "x": pa.array([None if i % 4 == 0 else float(i) for i in range(20)], pa.float64()),
    "k": pa.array([None if i % 6 == 0 else i for i in range(20)], pa.int64()),
})
pq.write_table(nullable, os.path.join(out, "nulls_plain.parquet"), use_dictionary=False, compression="NONE")
# optional column with no actual nulls
nonull = pa.table({"x": pa.array([float(i) for i in range(20)], pa.float64())})
pq.write_table(nonull, os.path.join(out, "optional_no_nulls.parquet"), compression="NONE")
# unsupported: strings
pq.write_table(pa.table({"s": pa.array(["a", "bb"])}), os.path.join(out, "unsupported_string.parquet"))
print("ok")
# pyarrow defaults (no overrides): what most real-world writers produce.
t2 = pa.table({"i64": pa.array([i * 7 - 300 for i in range(n)], pa.int64()),
               "i32": pa.array([i - 500 for i in range(n)], pa.int32()),
               "f64": pa.array([i * 0.25 - 10.0 for i in range(n)], pa.float64()),
               "f32": pa.array([i * 0.5 for i in range(n)], pa.float32()),
               "flag": pa.array([i % 3 == 0 for i in range(n)], pa.bool_()),
               "cat": pa.array([i % 5 for i in range(n)], pa.int64())})
pq.write_table(t2, os.path.join(out, "default_writer.parquet"))
