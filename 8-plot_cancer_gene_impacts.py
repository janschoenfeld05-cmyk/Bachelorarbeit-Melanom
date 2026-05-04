#!/usr/bin/env python3
"""
Erstellt Grafiken fuer die Analyse von Protoonkogenen und Tumorsuppressorgenen.

Die Plots lesen die Ergebnisse aus `results/cancer_gene_impact_analysis` ein
und erzeugen eine kleine Sammlung gut lesbarer Abbildungen:
1. einen Ueberblick ueber die Mutationsinterpretation pro Genrolle
2. die wichtigsten Protoonkogene
3. die wichtigsten Tumorsuppressorgene
4. einen kombinierten Ueberblick mit allen drei Ansichten
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import numpy as np
import pandas as pd

BASE_DIR = Path(__file__).resolve().parent.parent

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-codex")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp/xdg-cache-codex")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

matplotlib.rcParams["font.family"] = "DejaVu Sans"
matplotlib.rcParams["axes.unicode_minus"] = False


ROLE_ORDER = ["Proto-Onkogen", "Tumorsuppressorgen", "Kontextabhaengig"]
ROLE_PROTO = "Proto-Onkogen"
ROLE_TSG = "Tumorsuppressorgen"
ROLE_DUAL = "Kontextabhaengig"
ROLE_DISPLAY_LABELS = {
    ROLE_PROTO: "Protoonkogene",
    ROLE_TSG: "Tumorsuppressorgene",
    ROLE_DUAL: "kontextabhängige",
}
ROLE_COLORS = {
    "Proto-Onkogen": "#D1495B",
    "Tumorsuppressorgen": "#2D6A9F",
    "Kontextabhaengig": "#7A8B99",
}
CONSISTENCY_COLORS = {
    "role_consistent": "#2A9D8F",
    "unclear": "#E9C46A",
    "role_inconsistent": "#E76F51",
    "context_dependent": "#6D597A",
    "neutral_or_noncoding": "#C9CED6",
}
CONSISTENCY_LABELS = {
    "role_consistent": "rollenkonsistent",
    "unclear": "Unklar proteinverändernd",
    "role_inconsistent": "rolleninkonsistent",
    "context_dependent": "Kontextabhängig",
    "neutral_or_noncoding": "Neutral/nicht-kodierend",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Erstellt Grafiken fuer die Krebsgen-Impact-Analyse."
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=BASE_DIR / "results" / "cancer_gene_impact_analysis",
        help="Ordner mit den CSV-Ergebnissen der Krebsgen-Analyse.",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=BASE_DIR / "results" / "cancer_gene_impact_analysis" / "figures",
        help="Ausgabeordner fuer die Grafiken.",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=15,
        help="Wie viele Gene pro Genklasse in den Top-Plots gezeigt werden sollen.",
    )
    return parser.parse_args()


def load_tables(input_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Laedt die benoetigten Ergebnistabellen."""
    summary = pd.read_csv(input_dir / "cancer_gene_summary.csv")
    detailed = pd.read_csv(input_dir / "cancer_gene_mutations_detailed.csv")
    report = pd.read_csv(input_dir / "signature_share_summary.csv") if False else pd.DataFrame()
    _ = report
    return summary, detailed, pd.DataFrame()


def prepare_numeric_columns(summary: pd.DataFrame, detailed: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Konvertiert relevante Spalten von Text in numerische Werte."""
    summary = summary.copy()
    detailed = detailed.copy()

    summary_numeric_cols = [
        "mutation_count",
        "patient_count",
        "sample_count",
        "role_consistent_count",
        "unclear_count",
        "role_inconsistent_count",
        "context_dependent_count",
        "neutral_or_noncoding_count",
        "high_impact_count",
        "moderate_impact_count",
        "low_impact_count",
        "modifier_impact_count",
    ]
    for column in summary_numeric_cols:
        summary[column] = pd.to_numeric(summary[column], errors="coerce").fillna(0)

    return summary, detailed


def save_figure(fig: plt.Figure, path: Path) -> None:
    """Speichert jede Abbildung als PNG und SVG."""
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(path.with_suffix(".svg"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_top_genes(ax: plt.Axes, data: pd.DataFrame, role: str, top_n: int) -> None:
    """
    Zeigt Top-Gene als ueberlagerte Balken:
    - hell: alle Mutationen im Gen
    - farbig: rollen-konsistente Mutationen
    """
    subset = (
        data[data["gene_role"] == role]
        .sort_values(["role_consistent_count", "mutation_count", "patient_count"], ascending=False)
        .head(top_n)
        .sort_values("role_consistent_count", ascending=True)
    )

    ax.barh(
        subset["gene"],
        subset["mutation_count"],
        color="#D9DEE5",
        edgecolor="white",
        height=0.8,
        label="Alle Mutationen",
    )
    ax.barh(
        subset["gene"],
        subset["role_consistent_count"],
        color=ROLE_COLORS[role],
        edgecolor="white",
        height=0.8,
        label="rollenkonsistent",
    )

    for _, row in subset.iterrows():
        ax.text(
            row["mutation_count"] + 0.8,
            row["gene"],
            f"P={int(row['patient_count'])}",
            va="center",
            ha="left",
            fontsize=8,
            color="#4A4A4A",
        )

    if role == ROLE_PROTO:
        title = f"Top {len(subset)} Onkogene"
    elif role == ROLE_TSG:
        title = f"Top {len(subset)} Tumorsuppressorgene"
    else:
        title = f"Top {len(subset)} {role}"

    ax.set_title(title, fontsize=13, weight="bold")
    ax.set_xlabel("Mutationen")
    ax.set_ylabel("")
    ax.grid(axis="x", linestyle=":", linewidth=0.8, alpha=0.5)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def plot_role_consistency_overview(ax: plt.Axes, detailed: pd.DataFrame) -> None:
    """Stellt die Interpretationsklassen pro Genrolle als gestapelte Balken dar."""
    counts = (
        detailed.groupby(["gene_role", "role_consistency"])
        .size()
        .reset_index(name="count")
        .pivot(index="gene_role", columns="role_consistency", values="count")
        .fillna(0)
    )
    counts = counts.reindex(ROLE_ORDER).fillna(0)
    display_labels = [ROLE_DISPLAY_LABELS.get(role, role) for role in counts.index]

    left = np.zeros(len(counts), dtype=float)
    for consistency in CONSISTENCY_LABELS:
        if consistency not in counts.columns:
            values = np.zeros(len(counts), dtype=float)
        else:
            values = counts[consistency].to_numpy(dtype=float)
        ax.barh(
            display_labels,
            values,
            left=left,
            color=CONSISTENCY_COLORS[consistency],
            edgecolor="white",
            height=0.8,
            label=CONSISTENCY_LABELS[consistency],
        )
        left = left + values

    ax.set_title("Mutationsinterpretation nach Genrolle", fontsize=13, weight="bold")
    ax.set_xlabel("Anzahl Mutationen")
    ax.set_ylabel("")
    ax.grid(axis="x", linestyle=":", linewidth=0.8, alpha=0.5)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def make_individual_plots(summary: pd.DataFrame, detailed: pd.DataFrame, outdir: Path, top_n: int) -> None:
    """Erstellt die einzelnen Hauptplots."""
    fig, ax = plt.subplots(figsize=(10, 6.8), constrained_layout=True)
    plot_role_consistency_overview(ax, detailed)
    ax.legend(loc="lower right", frameon=False, fontsize=9)
    save_figure(fig, outdir / "role_consistency_overview")

    fig, ax = plt.subplots(figsize=(10, 7.5), constrained_layout=True)
    plot_top_genes(ax, summary, ROLE_PROTO, top_n)
    ax.legend(loc="lower right", frameon=False, fontsize=9)
    save_figure(fig, outdir / "top_proto_oncogenes")

    fig, ax = plt.subplots(figsize=(10, 7.5), constrained_layout=True)
    plot_top_genes(ax, summary, ROLE_TSG, top_n)
    ax.legend(loc="lower right", frameon=False, fontsize=9)
    save_figure(fig, outdir / "top_tumor_suppressor_genes")


def make_overview_plot(summary: pd.DataFrame, detailed: pd.DataFrame, outdir: Path, top_n: int) -> None:
    """Kombiniert die wichtigsten Ansichten in einer Uebersichtsabbildung."""
    fig, axes = plt.subplots(1, 3, figsize=(22, 8), constrained_layout=True)

    plot_role_consistency_overview(axes[0], detailed)
    plot_top_genes(axes[1], summary, ROLE_PROTO, top_n)
    plot_top_genes(axes[2], summary, ROLE_TSG, top_n)

    consistency_handles = [
        Patch(facecolor=CONSISTENCY_COLORS[key], label=CONSISTENCY_LABELS[key])
        for key in CONSISTENCY_LABELS
    ]
    gene_handles = [
        Patch(facecolor="#D9DEE5", label="Alle Mutationen"),
        Patch(facecolor=ROLE_COLORS[ROLE_PROTO], label="rollenkonsistent"),
    ]

    axes[0].legend(handles=consistency_handles, loc="lower right", frameon=False, fontsize=8)
    axes[1].legend(handles=gene_handles, loc="lower right", frameon=False, fontsize=8)
    axes[2].legend(handles=gene_handles, loc="lower right", frameon=False, fontsize=8)

    fig.suptitle(
        "Auswirkungen somatischer Mutationen auf Protoonkogene und Tumorsuppressorgene",
        fontsize=16,
        weight="bold",
    )
    save_figure(fig, outdir / "cancer_gene_impact_overview")


def main() -> None:
    args = parse_args()
    summary, detailed, _ = load_tables(args.input_dir)
    summary, detailed = prepare_numeric_columns(summary, detailed)

    args.outdir.mkdir(parents=True, exist_ok=True)
    make_individual_plots(summary, detailed, args.outdir, args.top_n)
    make_overview_plot(summary, detailed, args.outdir, args.top_n)

    print(f"Fertig. Grafiken liegen unter: {args.outdir}")


if __name__ == "__main__":
    main()
