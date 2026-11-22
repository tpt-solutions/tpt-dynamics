"""Same analysis as `tptdyn analyze data/events-10gb.csv --filter "pt > 5"
--hist pt,100,0,100 --hist eta,100,-5,5`, using ROOT RDataFrame.

Usage: python root_analyze.py <csv> <threads>   (threads=1 runs single-threaded)
Prints wall time, row counts and histogram entry counts.
"""
import sys
import time

import ROOT

csv, threads = sys.argv[1], int(sys.argv[2])
if threads > 1:
    ROOT.EnableImplicitMT(threads)

t0 = time.perf_counter()
df = ROOT.RDF.FromCSV(csv, True)
n_rows = df.Count()
kept = df.Filter("pt > 5", "pt>5")
n_kept = kept.Count()
h_pt = kept.Histo1D(ROOT.RDF.TH1DModel("pt", "pt", 100, 0, 100), "pt")
h_eta = kept.Histo1D(ROOT.RDF.TH1DModel("eta", "eta", 100, -5, 5), "eta")
rows = n_rows.GetValue()
kept_rows = n_kept.GetValue()
pt_entries = h_pt.GetEntries()
eta_entries = h_eta.GetEntries()
elapsed = time.perf_counter() - t0

print(f"threads={threads} wall_s={elapsed:.1f}")
print(f"rows={rows} kept={kept_rows} pt_entries={pt_entries:.0f} eta_entries={eta_entries:.0f}")
