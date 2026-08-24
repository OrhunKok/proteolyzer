"""What a DIA window design did to a precursor's envelope.

The cases came with the function from streamlit-DO-MS, where they ran as a script
of assertions rather than as tests. The window design here is deliberately small
-- two windows meeting at 700, each over its own slice of ion mobility -- because
every question worth asking is about an edge.
"""

import pandas as pd

from proteolyzer.core import envelope_split


def windows() -> pd.DataFrame:
    """Two windows meeting at 700, each covering its own slice of mobility."""
    return pd.DataFrame(
        {
            "Start Mass [m/z]": [400.0, 700.0],
            "End Mass [m/z]": [700.0, 1000.0],
            "Start IM [1/K0]": [0.7, 0.9],
            "End IM [1/K0]": [0.9, 1.1],
        }
    )


def precursors() -> pd.DataFrame:
    """One precursor per case, in the order the cases below read them."""
    return pd.DataFrame(
        {
            "Precursor.Mz": [500.0, 699.9, 699.9, 999.5, 1200.0, 500.0],
            # Categorical because that is how the loader hands charge over.
            "Precursor.Charge": pd.Categorical([2, 1, 4, 2, 2, 2]),
            "IM": [0.8, 0.8, 0.8, 1.0, 1.0, 5.0],
        }
    )


def test_a_precursor_well_inside_its_window_is_whole():
    assert envelope_split(precursors(), windows())[0] == "Intact"


def test_an_isotope_over_the_edge_splits_the_envelope():
    """M0 is in the window and M+1 is not, so MS2 saw part of what MS1 measured."""
    assert envelope_split(precursors(), windows())[1] == "Split"


def test_charge_decides_how_far_the_envelope_reaches():
    """Same m/z as the case above at charge 4: the isotopes sit closer, and it
    still crosses the edge. The spacing is ISOTOPE_STEP / charge, not a constant."""
    assert envelope_split(precursors(), windows())[2] == "Split"


def test_the_last_windows_upper_edge_splits_envelopes_too():
    """Nothing is above the top window, so an isotope past it is fragmented in no
    window at all -- which is the same missing signal, not a special case."""
    assert envelope_split(precursors(), windows())[3] == "Split"


def test_a_precursor_no_window_covers_is_left_out_rather_than_guessed_at():
    status = envelope_split(precursors(), windows())
    assert pd.isna(status[4]) and pd.isna(status[5])


def test_mobility_counts_as_much_as_mz_where_the_design_gives_it():
    """Right m/z, wrong mobility: some other window's job. Drop the mobility and
    the same precursor is answered on m/z alone."""
    without = precursors().drop(columns="IM")
    assert pd.isna(envelope_split(precursors(), windows())[5])
    assert envelope_split(without, windows())[5] == "Intact"


def test_no_window_design_claims_nothing():
    empty = windows().iloc[:0]
    assert all(pd.isna(value) for value in envelope_split(precursors(), empty))


def test_a_frame_without_the_columns_is_answered_with_nothing():
    """``Precursor.Mz`` is not in the core's DIA-NN subset, that subset being a
    pipeline's rather than a dashboard's, so a caller may well not have it."""
    lacking = precursors().drop(columns="Precursor.Mz")
    answer = envelope_split(lacking, windows())

    assert len(answer) == len(lacking)
    assert all(pd.isna(value) for value in answer)


def test_the_same_precursor_identified_again_gets_the_same_answer():
    """Deduplicating to ask the question once has to give every identification of
    that precursor the same answer, on the index it arrived with."""
    repeated = precursors().iloc[[1, 1, 0]].set_axis([7, 8, 9])
    answer = envelope_split(repeated, windows())

    assert list(answer) == ["Split", "Split", "Intact"]
    assert list(answer.index) == [7, 8, 9]
