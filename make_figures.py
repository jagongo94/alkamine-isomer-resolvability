#!/usr/bin/env python3
"""
Figures for the alkamine isomer-resolvability analysis.

Panel a   MS/MS cosine against structural Tanimoto, one point per distinct
          spectrum, coloured by the structural similarity of its best partner.
Panels b+ mirror plots of representative pairs, listed in examples.tsv.

Run run_analysis.py first, then:

    python make_figures.py --mgf-dir data --results results --out figures
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

import alkamine as ak
from run_analysis import LIBRARIES

# structural similarity bins: label, lower bound, upper bound, colour
BINS = [
    ("extremely similar (>0.9)", 0.9, 1.01, "#1b3f8b"),
    ("very highly similar (0.8-0.9)", 0.8, 0.9, "#5b9bd5"),
    ("highly similar (0.6-0.8)", 0.6, 0.8, "#e8c33c"),
    ("similar (0.4-0.6)", 0.4, 0.6, "#e8903c"),
    ("moderately similar (0.2-0.4)", 0.2, 0.4, "#d1462f"),
    ("weakly similar (<0.2)", 0.0, 0.2, "#6b3024"),
]
PINK, GREY, TOP, BOTTOM, UNMATCHED = "#e06ba8", "#9a9a94", "#2b5fa8", "#e8732a", "#bdbcb6"

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 7,
    "axes.edgecolor": "#9a998f", "axes.linewidth": 0.6,
    "xtick.color": "#5f5e58", "ytick.color": "#5f5e58",
    "pdf.fonttype": 42, "svg.fonttype": "none",
})


def short_name(name: str) -> str:
    """Strip the library's annotation suffix from a compound name."""
    return re.sub(r"\s*\(known isomers.*$", "", str(name)).split(" [M")[0].strip()


def save(figure, path: Path) -> None:
    for extension in ("svg", "pdf", "png"):
        figure.savefig(path.with_suffix(f".{extension}"), bbox_inches="tight",
                       facecolor="white", dpi=400 if extension == "png" else None)
    plt.close(figure)
    print(f"  {path.name}.[svg|pdf|png]")


def panel_a(table: pd.DataFrame, out: Path) -> None:
    """Cosine against Tanimoto, one point per distinct spectrum."""
    figure, axes = plt.subplots(figsize=(8.2, 5.0), facecolor="white")
    below = table[table.cosine < ak.COSINE_CUTOFF]
    above = table[table.cosine >= ak.COSINE_CUTOFF]

    axes.scatter(below.cosine, below.tanimoto, s=9, c=GREY, lw=0, rasterized=True)
    flat = above[above.same_flat_structure]
    axes.scatter(flat.cosine, flat.tanimoto, s=13, c=PINK, lw=0, rasterized=True)
    for _, low, high, colour in BINS:
        chunk = above[(~above.same_flat_structure)
                      & (above.tanimoto >= low) & (above.tanimoto < high)]
        axes.scatter(chunk.cosine, chunk.tanimoto, s=13, c=colour, lw=0, rasterized=True)

    for level in (0.2, 0.4, 0.6, 0.8, 0.9):
        axes.axhline(level, ls="--", color="#c9c8c2", lw=0.5, zorder=0)
    axes.axvline(ak.COSINE_CUTOFF, ls="--", color="#5f5e58", lw=0.9)

    axes.set_xlim(0, 1.02)
    axes.set_ylim(0, 1.04)
    axes.set_xlabel("MS/MS similarity (cosine, fragment ions only)", fontsize=10)
    axes.set_ylabel(f"Structural similarity (Tanimoto, Morgan radius {ak.MORGAN_RADIUS})",
                    fontsize=10)
    axes.tick_params(labelsize=9)
    axes.spines[["top", "right"]].set_visible(False)
    axes.set_axisbelow(True)
    axes.set_title("MS/MS similarity against structural similarity",
                   loc="left", fontsize=11.5, color="#1f1f1d", pad=16)
    axes.text(0, 1.025,
              "each point is one distinct spectrum, placed by its best "
              f"same-formula partner (n = {len(table):,})",
              transform=axes.transAxes, fontsize=8.5, color="#5f5e58")

    handles = [Line2D([], [], marker="o", ls="", ms=7, color=PINK,
                      label=f"same flat structure (n = {len(flat):,})")]
    for label, low, high, colour in BINS:
        count = len(above[(~above.same_flat_structure)
                          & (above.tanimoto >= low) & (above.tanimoto < high)])
        handles.append(Line2D([], [], marker="o", ls="", ms=7, color=colour,
                              label=f"{label}  n = {count:,}"))
    handles.append(Line2D([], [], marker="o", ls="", ms=7, color=GREY,
                          label=f"cosine < {ak.COSINE_CUTOFF:.1f} (n = {len(below):,})"))
    axes.legend(handles=handles, frameon=False, fontsize=9, loc="upper left",
                bbox_to_anchor=(1.015, 1.0), labelspacing=0.75, borderaxespad=0.0)
    save(figure, out / "Panel_a_cosine_vs_tanimoto")


def matched_mask(mz_a: np.ndarray, mz_b: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Flag the fragments of each spectrum that pair with one in the other."""
    hit_a = np.zeros(len(mz_a), bool)
    hit_b = np.zeros(len(mz_b), bool)
    for i, mass in enumerate(mz_a):
        free = [j for j in np.where(np.abs(mz_b - mass) <= ak.FRAGMENT_TOL)[0]
                if not hit_b[j]]
        if free:
            hit_a[i] = True
            hit_b[free[0]] = True
    return hit_a, hit_b


def mirror(axes, row: pd.Series, spectra: dict[int, ak.Spectrum]) -> None:
    """Draw one head-to-tail comparison of a spectrum and its best partner."""
    def fragments(index: int):
        spectrum = spectra[index]
        keep = np.abs(spectrum.mz - spectrum.precursor_mz) > ak.PRECURSOR_WINDOW
        mz, intensity = spectrum.mz[keep], spectrum.intensity[keep]
        return mz, 100 * intensity / intensity.max()

    mz_top, int_top = fragments(int(row.spectrum))
    mz_bottom, int_bottom = fragments(int(row.partner))
    hit_top, hit_bottom = matched_mask(mz_top, mz_bottom)

    axes.vlines(mz_top[~hit_top], 0, int_top[~hit_top], color=UNMATCHED, lw=0.9)
    axes.vlines(mz_top[hit_top], 0, int_top[hit_top], color=TOP, lw=1.0)
    axes.vlines(mz_bottom[~hit_bottom], 0, -int_bottom[~hit_bottom], color=UNMATCHED, lw=0.9)
    axes.vlines(mz_bottom[hit_bottom], 0, -int_bottom[hit_bottom], color=BOTTOM, lw=1.0)
    axes.axhline(0, color="#7f7f7f", lw=0.7, zorder=0)

    everything = np.r_[mz_top, mz_bottom]
    span = max(everything.max() - everything.min(), 1)
    axes.set_xlim(everything.min() - span * 0.04, everything.max() + span * 0.06)
    axes.set_ylim(-175, 175)
    axes.set_yticks([-100, 0, 100])
    axes.set_yticklabels(["100", "0", "100"], fontsize=6)
    axes.tick_params(axis="x", labelsize=6)
    for spine in axes.spines.values():
        spine.set_color("#9a998f")
        spine.set_linewidth(0.6)
    axes.set_xlabel("m/z", fontsize=6.5, labelpad=1)
    axes.text(0.012, 0.975, short_name(row["name"]), transform=axes.transAxes,
              fontsize=6.2, color=TOP, va="top")
    axes.text(0.012, 0.025, short_name(row.partner_name), transform=axes.transAxes,
              fontsize=6.2, color=BOTTOM, va="bottom")
    axes.set_title(f"cosine {row.cosine:.2f} · Tanimoto {row.tanimoto:.2f} · "
                   f"{int(row.matched_peaks)} matched",
                   loc="left", fontsize=6.4, color="#1f1f1d", pad=3)


def choose(table: pd.DataFrame, name_a: str, name_b: str,
           target: float | None, low_cosine: bool) -> pd.Series:
    """Find the example pair, preferring the most fragment-rich comparison."""
    def contains(column, text):
        return column.str.contains(text, case=False, regex=False)

    short_a = table["name"].map(short_name)
    short_b = table.partner_name.map(short_name)
    hit = table[(contains(short_a, name_a) & contains(short_b, name_b))
                | (contains(short_a, name_b) & contains(short_b, name_a))]
    hit = hit[hit.cosine < ak.COSINE_CUTOFF] if low_cosine else hit[hit.cosine >= ak.COSINE_CUTOFF]
    if hit.empty:
        raise SystemExit(f"no pair found for {name_a} / {name_b}")
    if target is not None and not np.isnan(target):
        return hit.iloc[(hit.cosine - target).abs().argsort()].iloc[0]
    return hit.sort_values("matched_peaks", ascending=False).iloc[0]


def mirror_panels(table: pd.DataFrame, spectra: dict[int, ak.Spectrum],
                  examples: pd.DataFrame, out: Path) -> None:
    for panel, group in examples.groupby("panel", sort=False):
        rows = (len(group) + 1) // 2
        figure = plt.figure(figsize=(7.0, 1.85 * rows), facecolor="white")
        grid = figure.add_gridspec(rows, 2, hspace=0.72, wspace=0.20)
        for position, (_, example) in enumerate(group.iterrows()):
            axes = figure.add_subplot(grid[position // 2, position % 2])
            mirror(axes, choose(table, example.compound_a, example.compound_b,
                                example.target_cosine, bool(example.low_cosine)), spectra)
            if position % 2 == 0:
                axes.set_ylabel("Intensity (%)", fontsize=6.5)
        figure.suptitle(group.bin_label.iloc[0], x=0.005, ha="left",
                        fontsize=8.4, color="#1f1f1d", fontweight="bold",
                        y=1.0 + 0.055 / rows)
        save(figure, out / f"Panel_{panel}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mgf-dir", default="data", type=Path)
    parser.add_argument("--results", default="results", type=Path)
    parser.add_argument("--examples", default="examples.tsv", type=Path)
    parser.add_argument("--out", default="figures", type=Path)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    table = pd.read_csv(args.results / "best_partner.csv")
    spectra_list = ak.read_libraries(
        {name: args.mgf_dir / filename for name, filename in LIBRARIES.items()})
    spectra = {s.index: s for s in spectra_list}

    print("writing figures ...")
    panel_a(table, args.out)
    mirror_panels(table, spectra, pd.read_csv(args.examples, sep="\t"), args.out)


if __name__ == "__main__":
    main()
