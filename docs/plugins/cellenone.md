# cellenONE

Parses the exports a cellenONE writes during a single-cell preparation and maps
each isolated cell to the well it was dispensed into.

```python
import proteolyzer as pz

mapper = pz.cellenone.CoordinatesMapping(
    root_dir="cellenone_files", label_type="mTRAQ", plex=2
)

metadata = mapper.map_data()   # one row per isolated cell, with well and label
stats = mapper.map_stats()     # environment readings from the instrument logs

mapper.save("prep_results/")
```

A root directory is one preparation. The conditions inside it are all part of
the same preparation, and different experiments have different numbers of them.

## What it reads

Files are classified by their directory *and* their name, into mutually
exclusive buckets: cell geometry, dispense events, labelling runs and pickups.
Operators name labelling logs inconsistently, so classifying on the filename
alone silently misfiles them — which is why the classification is recorded in
the provenance.

Cells are keyed on their position and imaging channel. The `Teg` column names
the channel a row was measured in, so a cell imaged in transmission and a
fluorescence channel appears more than once; the channels are pivoted into
per-channel columns rather than left as extra rows.

Label dispenses are collapsed into `Droplets` and `Dispenses` counts. Where the
instrument log does not say how many drops were dispensed, the count stays
missing rather than becoming zero.

## What `save` writes

`metadata.parquet`, `instrument_stats.parquet`, and a `provenance.jsonl` entry
naming the package version, the configuration, which logs were classified as
what, and how many wells clashed.

The metadata frame on its own says none of that. Keeping it matters because the
classification is the record of what was actually picked up.

## Well clashes

Two cells assigned to the same well, or a well assigned two labels of the same
plex, are flagged rather than silently resolved. Coordinate matching reports
the distance it matched at, so a cell matched loosely is visible as such.
