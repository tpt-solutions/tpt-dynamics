Event analysis project (created by `tptdyn init hep-analysis`)

  python make_events.py     writes events.csv (synthetic; standard library only)
  analyze.txt               command sheet with analyze and histogram examples

The sheet writes histogram files into hists/. Input can also be a Parquet file
or a tptdyn columnar/row-store file; see `tptdyn analyze --help`.
