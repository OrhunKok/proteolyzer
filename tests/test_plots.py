import shutil

import matplotlib
import numpy as np
import pandas as pd
import pytest

matplotlib.use("Agg")

from proteolyzer.plots import PlotBase, VolcanoPlot  # noqa: E402

# The default "science" theme renders text with LaTeX, so anything that lays
# text out needs a TeX installation. Those tests run on the stock theme.
needs_latex = pytest.mark.skipif(
    shutil.which("latex") is None, reason="requires a LaTeX installation"
)


@pytest.fixture(autouse=True)
def _fresh_figure():
    """Plots draw onto plt.gca(), so each test needs a clean canvas."""
    yield
    matplotlib.pyplot.close("all")


@pytest.fixture
def volcano_data() -> pd.DataFrame:
    rng = np.random.default_rng(0)
    return pd.DataFrame(
        {
            "log2FC": np.concatenate([rng.normal(2, 0.1, 5), rng.normal(-2, 0.1, 5)]),
            "pvalue": np.concatenate([np.full(5, 1e-6), np.full(5, 0.5)]),
            "Gene": [f"GENE{i}" for i in range(10)],
        }
    )


def test_plot_base_without_a_figure_is_inert(caplog):
    """Regression: attribute access on an unplotted PlotBase recursed forever."""
    base = PlotBase()
    assert base.ax is None
    base.save("unused.png")
    base.show()
    assert "No figure has been generated yet." in caplog.text
    assert "No plot has been generated yet." in caplog.text


def test_unknown_attributes_still_raise():
    base, missing = PlotBase(), "nonsense"
    with pytest.raises(AttributeError, match="no attribute 'nonsense'"):
        getattr(base, missing)


def test_unknown_theme_is_reported_but_does_not_fail(caplog):
    base = PlotBase(theme="not-a-real-theme")
    with base.plot_theme():
        pass
    assert "Could not apply theme" in caplog.text


def test_plot_theme_propagates_errors_from_its_body():
    with pytest.raises(ZeroDivisionError), PlotBase().plot_theme():
        raise ZeroDivisionError("from inside the theme context")


def test_filter_kwargs_keeps_only_accepted_arguments():
    def target(a, b=None):
        return a, b

    assert PlotBase().filter_kwargs(target, {"a": 1, "c": 2}) == {"a": 1}


def test_volcano_plot_classifies_regulation(volcano_data):
    plot = VolcanoPlot(volcano_data, x="log2FC", y="pvalue", label="Gene")

    counts = plot.data["Regulation"].value_counts()
    assert counts["up"] == 5
    assert counts["notsig"] == 5
    assert plot.ax is not None
    # -log10 transform applied, threshold drawn at -log10(0.05).
    assert plot.data["pvalue"].max() == pytest.approx(6.0)
    assert plot.signif == pytest.approx(-np.log10(0.05))


def test_volcano_plot_detects_pre_transformed_pvalues(volcano_data):
    transformed = volcano_data.assign(pvalue=-np.log10(volcano_data["pvalue"]))
    plot = VolcanoPlot(transformed, x="log2FC", y="pvalue")
    assert plot.data["pvalue"].max() == pytest.approx(6.0)


def test_volcano_plot_delegates_to_the_axes(volcano_data):
    plot = VolcanoPlot(volcano_data, x="log2FC", y="pvalue")
    plot.set_xlabel("fold change")
    assert plot.ax.get_xlabel() == "fold change"


def test_volcano_plot_can_label_points(volcano_data):
    plot = VolcanoPlot(
        volcano_data, x="log2FC", y="pvalue", label="Gene", theme="default"
    )
    plot.label_points()
    assert any(text.get_text().startswith("GENE") for text in plot.ax.texts)


def test_volcano_plot_labels_only_the_most_extreme_points(volcano_data, caplog):
    plot = VolcanoPlot(
        volcano_data, x="log2FC", y="pvalue", label="Gene", theme="default"
    )
    plot.label_points(signif_filter=False, max_label=4)
    assert "Too many points to label" in caplog.text
    assert len([t for t in plot.ax.texts if t.get_text().startswith("GENE")]) == 4


def test_volcano_plot_saves_and_closes(volcano_data, tmp_path):
    plot = VolcanoPlot(volcano_data, x="log2FC", y="pvalue", theme="default")
    out = tmp_path / "volcano.png"
    plot.save(out)
    assert out.exists()


@needs_latex
def test_volcano_plot_renders_with_the_default_science_theme(volcano_data, tmp_path):
    plot = VolcanoPlot(volcano_data, x="log2FC", y="pvalue")
    plot.add_data_point_count()
    out = tmp_path / "volcano_science.png"
    plot.save(out)
    assert out.exists()


# ------------------------------------------------------- volcano refinements


def test_a_p_value_of_zero_is_drawn_rather_than_lost(caplog):
    """Regression: -log10(0) is inf, counted in the box but not on the axes."""
    data = pd.DataFrame({"lfc": [2.0, -2.0], "p": [0.0, 1e-8]})

    plot = VolcanoPlot(data, x="lfc", y="p", theme="default")

    assert np.isfinite(plot.data["p"]).all()
    # Held at the most extreme finite value, so it is visibly the top point.
    assert plot.data["p"].tolist() == [8.0, 8.0]
    assert "p-value(s) of exactly 0" in caplog.text
    assert plot.data["Regulation"].tolist() == ["up", "down"]


def test_all_p_values_zero_still_plots():
    plot = VolcanoPlot(
        pd.DataFrame({"lfc": [1.0, -1.0], "p": [0.0, 0.0]}),
        x="lfc",
        y="p",
        theme="default",
    )
    assert np.isfinite(plot.data["p"]).all()


def test_an_effect_threshold_excludes_small_changes(volcano_data):
    """A tiny fold change with a great p-value is not a hit."""
    data = pd.concat(
        [
            volcano_data,
            pd.DataFrame({"log2FC": [0.05], "pvalue": [1e-9], "Gene": ["X"]}),
        ],
        ignore_index=True,
    )

    without = VolcanoPlot(data, x="log2FC", y="pvalue", theme="default")
    with_threshold = VolcanoPlot(
        data, x="log2FC", y="pvalue", effect_threshold=1.0, theme="default"
    )

    assert without.data.set_index("Gene").loc["X", "Regulation"] == "up"
    assert with_threshold.data.set_index("Gene").loc["X", "Regulation"] == "notsig"


def test_the_effect_threshold_is_drawn(volcano_data):
    plot = VolcanoPlot(
        volcano_data, x="log2FC", y="pvalue", effect_threshold=1.5, theme="default"
    )
    vertical = [
        line.get_xdata()
        for line in plot.ax.get_lines()
        if len(line.get_xdata()) == 2 and line.get_xdata()[0] == line.get_xdata()[1]
    ]
    assert sorted(xdata[0] for xdata in vertical) == [-1.5, 0.0, 1.5]


@pytest.mark.parametrize("transformed", [True, False])
def test_the_transform_can_be_stated_explicitly(transformed):
    """The heuristic cannot tell these apart; the caller can."""
    data = pd.DataFrame({"lfc": [1.0, -1.0], "y": [0.5, 0.9]})

    plot = VolcanoPlot(data, x="lfc", y="y", transformed=transformed, theme="default")

    if transformed:
        assert plot.data["y"].tolist() == [0.5, 0.9]
    else:
        assert plot.data["y"].tolist() == pytest.approx([0.30103, 0.045757], rel=1e-3)


def test_the_count_box_is_placed_off_the_x_axis_only(volcano_data):
    """Regression: the label offset was computed from the y coordinate."""
    plot = VolcanoPlot(volcano_data, x="log2FC", y="pvalue", theme="default")
    plot._add_delta_count_box(box_position=(0.1, 0.9), box_width=0.3)

    labels = [t for t in plot.ax.texts if "$>$0" in t.get_text()]
    assert labels
    # 0.1 + 0.3 * 0.1, not 0.1 + (0.9 - 0.1) * 0.1.
    assert labels[-1].get_position()[0] == pytest.approx(0.13)


def test_counting_points_on_empty_axes_is_reported(caplog):
    import matplotlib.pyplot as plt

    from proteolyzer.plots.relational import RelPlot

    plot = RelPlot(theme="default")
    plot.ax = plt.figure().add_subplot()

    plot.add_data_point_count()

    assert "needs a scatter on the axes" in caplog.text
    assert len(plot.ax.texts) == 0


def test_a_symmetric_x_axis_can_be_asked_for(volcano_data):
    plot = VolcanoPlot(
        volcano_data, x="log2FC", y="pvalue", symmetric_x=True, theme="default"
    )
    low, high = plot.ax.get_xlim()
    assert low == pytest.approx(-high)


def test_control_options_are_not_forwarded_to_seaborn(volcano_data):
    """Regression: options read out of **kwargs reached scatterplot and raised."""
    plot = VolcanoPlot(
        volcano_data,
        x="log2FC",
        y="pvalue",
        symmetric_x=True,
        delta_text_size=9,
        theme="default",
    )
    assert plot.delta_text_size == 9


def test_seaborn_keywords_still_reach_the_plot(volcano_data):
    plot = VolcanoPlot(
        volcano_data, x="log2FC", y="pvalue", marker="^", theme="default"
    )
    assert plot.ax.collections
