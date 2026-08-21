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
