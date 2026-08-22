# cellenONE

The instrument writes a directory per run: a log for every step of the
preparation, an isolated-cells table beside the imaging step, and a great many
images. What a consumer wants out of it is one row per cell — where it landed,
which label it got, how big it was, and what the chamber was doing at the time.

```python
from proteolyzer.cellenone import CoordinatesMapping

run = CoordinatesMapping(files, label_type="mTRAQ", plex=3)

run.map_data()        # a row per cell: position, well, label, geometry
run.map_stats()       # the chamber readings, stamped with the step
run.parse_droplets()  # every droplet dispensed, by step
```

`files` is what the run directory holds, as paths or as uploaded file objects
with a `name`.

## Which file is which step

Nothing in the directory reliably says. Folders are named after the step on a
good day and numbered `01`–`06` on a normal one, and a directory that has been
worked in holds test runs, second attempts, and steps named after their reagent
rather than their purpose. So the step is worked out from the file itself, in
this order:

1. **The name**, where it says something.
2. **The `Run ID:` line inside the log**, which is the folder name over again —
   so a renamed folder still says what it held. The method summary next to it is
   not usable for this: it lists every reagent the run touched, and a quench log
   names the label it is quenching.
3. **What the file is made of**: the cells table by its columns, pickup by having
   one well per cell, labelling by the reagent row it draws from, DMSO and quench
   by going down before and after the label.

A step that cannot be named is left as a plain dispense record rather than
pushed into one that is already accounted for. `run.stages` says what each file
was taken for, and passing that dictionary back in overrides it.

## Two things worth knowing

**A step that ran twice ran again because the first attempt was abandoned.**
Where both attempts recorded the same cell, the later one is what became of it;
where the redo covered part of the slide, the cells it left alone keep the first
attempt's record. `run.contributed` says how many cells each file is left
accounting for.

**A geoprops row is a cell *per imaging channel*.** A slide imaged in
transmission and one fluorescence channel has two rows for every cell on it, so
counting rows doubles the cell count and halves every geometry average. The
transmission row is the cell, and each other channel joins onto it under its own
name: `Diameter.Green` beside `Diameter`. Nothing is discarded — in a sorting
experiment the fluorescence is what the experiment was for.
