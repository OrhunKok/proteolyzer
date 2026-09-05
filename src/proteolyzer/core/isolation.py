"""What a DIA isolation-window design did to a precursor's isotopic envelope.

A DIA method fragments everything inside a window of m/z, so the MS2 spectrum a
precursor is identified from depends on where the window edges fell relative to
its isotopes. Where part of the envelope lands outside the window, the fragments
come from a different slice of the signal than the MS1 quantitation was measured
over -- which is a property of the design and the precursor, not of the sample,
and worth ruling out before a quantitative difference is attributed to one.

The columns this reads are not all ones the core loads by default:
``Precursor.Mz`` is not in the DIA-NN report subset, because that subset is a
pipeline's rather than a dashboard's. So the caller supplies the frame, and a
frame without the columns is answered with nothing rather than a guess.

Moved here from streamlit-DO-MS, which had the only copy of it and could ask the
question only from inside a dashboard. See DECISIONS.md.
"""

import numpy as np
import pandas as pd

#: The gap between one isotope of a peptide and the next, as a neutral mass: 13C
#: less 12C. Divided by the charge it becomes a step in m/z.
ISOTOPE_STEP: float = 1.00336

#: How much of the envelope has to be isolated together for the precursor to be
#: fragmented as one thing. M0 to M+2 carries all but a little of the signal.
ENVELOPE_ISOTOPES: int = 2

#: What the frame is asked for. Everything else it may carry is ignored.
PRECURSOR_COLUMNS: frozenset[str] = frozenset({"Precursor.Mz", "Precursor.Charge"})

#: The window design's own columns, which come from the instrument method rather
#: than from a search engine, so they are named as the method writes them.
WINDOW_MASS_COLUMNS: tuple[str, str] = ("Start Mass [m/z]", "End Mass [m/z]")
WINDOW_MOBILITY_COLUMNS: tuple[str, str] = ("Start IM [1/K0]", "End IM [1/K0]")


def envelope_room(
    report: pd.DataFrame,
    windows: pd.DataFrame,
    isotopes: int = ENVELOPE_ISOTOPES,
) -> pd.DataFrame:
    """Which window isolated each precursor, and how much room it left its envelope.

    'Window' is a *position* in ``windows``, so it is ``windows.iloc[w]`` and not
    ``windows.loc[w]`` -- the two part company the moment a caller hands over a
    design it filtered, and ``.loc`` is then either a ``KeyError`` or, where the
    labels happen to overlap the positions, the wrong window with nothing said.
    It is -1 where no window isolated the precursor, which is a sentinel a
    position can carry and a label could not.

    'Room' is the m/z from M+2 to that window's upper edge, negative where part
    of the envelope was fragmented elsewhere, NaN where no window isolated it.
    Both are aligned to the report's index.

    Windows overlap, so more than one may isolate a precursor by m/z and, where
    the design carries mobility, by ion mobility too. Walking them in ascending
    upper edge and letting the last one to claim a precursor win picks the same
    window `envelope_split` used to call 'whole' on -- the one with the most
    room, because ``max(ends[isolated]) >= top`` is the same statement as
    ``any(isolated & (top <= ends))`` -- without building a precursors-by-windows
    matrix to get there.
    """

    def numeric(column: str) -> pd.Series:
        # Charge and mobility arrive categorical from the loader, and to_numeric
        # will not take a categorical directly.
        return pd.to_numeric(report[column].astype("object"), errors="coerce")

    empty = pd.DataFrame(
        {
            "Window": pd.Series(-1, index=report.index, dtype="int64"),
            "Room": pd.Series(np.nan, index=report.index, dtype="float64"),
        }
    )

    if not PRECURSOR_COLUMNS.issubset(report.columns) or windows.empty:
        return empty

    # Only worth keying on mobility when the design carries it -- otherwise it is
    # measured per identification and barely dedupes at all.
    mobile = set(WINDOW_MOBILITY_COLUMNS).issubset(windows.columns) and (
        "IM" in report.columns
    )

    # One row per precursor rather than per identification: the answer depends on
    # the m/z, the charge and the mobility, none of which vary by run or channel,
    # and a report can carry a million rows against a few thousand precursors.
    columns = {name: numeric(name) for name in ("Precursor.Mz", "Precursor.Charge")}
    if mobile:
        columns["IM"] = numeric("IM")

    frame = pd.DataFrame(columns, index=report.index)
    unique = frame.drop_duplicates()

    mz = unique["Precursor.Mz"].to_numpy(dtype=float)
    charge = unique["Precursor.Charge"].to_numpy(dtype=float)
    top = mz + isotopes * ISOTOPE_STEP / charge

    starts = windows[WINDOW_MASS_COLUMNS[0]].to_numpy(dtype=float)
    ends = windows[WINDOW_MASS_COLUMNS[1]].to_numpy(dtype=float)

    mobility = low = high = None
    if mobile and (unique["IM"] > 0).any():
        # A PASEF window is a patch of the m/z by mobility plane, so both have to
        # match. Isotopes share the precursor's charge and so its mobility.
        mobility = unique["IM"].to_numpy(dtype=float)
        low = windows[WINDOW_MOBILITY_COLUMNS[0]].to_numpy(dtype=float)
        high = windows[WINDOW_MOBILITY_COLUMNS[1]].to_numpy(dtype=float)

    window_idx = np.full(len(unique), -1, dtype=np.int64)
    room = np.full(len(unique), np.nan, dtype=np.float64)

    for i in np.argsort(ends):
        isolated = (mz >= starts[i]) & (mz <= ends[i])
        if mobility is not None and low is not None and high is not None:
            isolated &= (mobility >= low[i]) & (mobility <= high[i])
        window_idx[isolated] = i
        room[isolated] = ends[i] - top[isolated]

    unique = unique.assign(Window=window_idx, Room=room)

    # Back onto every identification of each precursor. A left merge keeps the
    # report's row order, so the values line up with the index they came from.
    spread = frame.merge(unique, on=list(frame.columns), how="left")

    return pd.DataFrame(
        {"Window": spread["Window"].to_numpy(), "Room": spread["Room"].to_numpy()},
        index=report.index,
    )


def envelope_split(
    report: pd.DataFrame,
    windows: pd.DataFrame,
    isotopes: int = ENVELOPE_ISOTOPES,
) -> pd.Series:
    """Whether each precursor's isotopic envelope was isolated whole or in parts.

    A precursor is isolated by the window its monoisotopic m/z falls in, and its
    isotopes sit ``ISOTOPE_STEP / charge`` above that, one per 13C. Where an
    isotope lands past the window's upper edge it is fragmented in a different
    window, or in none, so the MS2 spectrum is missing part of the envelope the
    MS1 quantitation was measured from.

    A precursor counts as intact when some one window holds its whole envelope,
    rather than by looking up 'its' window first: schemes overlap, and a
    precursor isolated whole by any of them was isolated whole. Windows are taken
    as m/z by ion mobility rectangles, which is the same approximation a window
    overlay is drawn with.

    Returns ``'Intact'``, ``'Split'``, or ``None`` where the design covers no
    window the precursor could have come from -- aligned to the report's own
    index, so the answer can be assigned straight back onto it.

    A thin wrapper over `envelope_room`: the sign of its `Room` is this verdict,
    positive meaning the window's upper edge reached M+2 or past it.
    """

    room = envelope_room(report, windows, isotopes)

    # Built as an object array and wrapped, rather than assigned into a Series,
    # so the dtype handed back is the one every release up to v0.6.0 handed
    # back. Assigning into an object Series leaves it object, where wrapping
    # lets pandas infer -- which is `str` on pandas 3. Same values either way,
    # and a consumer that has already pinned a dtype should not have to find
    # out from a plot that it moved.
    isolated = (room["Window"] >= 0).to_numpy()
    envelope = np.full(len(report), None, dtype=object)
    envelope[isolated] = np.where(
        room["Room"].to_numpy()[isolated] >= 0, "Intact", "Split"
    )

    return pd.Series(envelope, index=report.index)
