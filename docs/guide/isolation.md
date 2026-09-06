# DIA isolation windows

A DIA method fragments everything inside a window of m/z, so the MS2 spectrum a
precursor is identified from depends on where the window edges fell relative to
its isotopes. Where part of the envelope lands outside the window, the fragments
come from a different slice of the signal than the MS1 quantitation was measured
over — a property of the method and the precursor rather than of the sample, and
worth ruling out before a quantitative difference is read as biology.

Two functions ask that question, and they are the same walk: one answers with a
verdict, the other with the window and the room it left.

## Was the envelope isolated whole?

```python
from proteolyzer.core import envelope_split

report.frame["Envelope"] = envelope_split(report.frame, windows)
```

`'Intact'`, `'Split'`, or missing where the design covers no window the precursor
could have come from — aligned to the report's own index, so it can be assigned
straight back onto it.

A precursor counts as intact when *some one* window holds its whole envelope,
rather than by looking up "its" window first: schemes overlap, and a precursor
isolated whole by any of them was isolated whole.

## Which window, and how much room?

```python
from proteolyzer.core import envelope_room

room = envelope_room(report.frame, windows)
room["Window"]   # a position in `windows`, -1 where none isolated it
room["Room"]     # m/z from M+2 to that window's upper edge
```

`Room` is negative where part of the envelope was fragmented elsewhere and
missing where no window isolated the precursor. The sign of it *is*
`envelope_split`'s verdict — that function is a thin wrapper over this one — so
there is one copy of the rule rather than two.

!!! warning "`Window` is a position, not a label"

    It is `windows.iloc[w]`, never `windows.loc[w]`. The two agree only while the
    design still carries its original index, and a caller that filtered its
    windows has one that does not — where `.loc` is a `KeyError`, or the wrong
    window with nothing said. `-1` is also why it has to be a position: a label
    could be `-1`.

This is what to reach for when the question is not "was this split" but "which
window should be widened": the window index is what a plot of the design needs,
and re-deriving it from `ISOTOPE_STEP` and `ENVELOPE_ISOTOPES` is what
`streamlit-DO-MS` was doing before v0.7.0.

## What the caller supplies

Both take the frame rather than a `Report`, and read four things:

| from | columns |
| --- | --- |
| the report | `Precursor.Mz`, `Precursor.Charge`, and `IM` where the design has mobility |
| the design | `Start Mass [m/z]` / `End Mass [m/z]`, and `Start IM [1/K0]` / `End IM [1/K0]` where it has them |

`Precursor.Mz` is not in every project's column subset — that subset being a
pipeline's rather than a dashboard's — so a frame without the columns is answered
with nothing rather than a guess. The window design is a property of the
instrument method, not of the search, so it arrives in a file of its own: DIA-NN's
`--export-windows` writes one beside every raw file it read.

Windows are taken as m/z by ion mobility rectangles, which is the same
approximation a window overlay is drawn with. Mobility is only read when the
design carries it — a precursor's isotopes share its charge and so its mobility.
