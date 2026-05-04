#!/usr/bin/env python3
"""
Erstellt eine breit formatierte SBS96-Heatmap aus einer Wide-CSV.

Erwartetes Eingabeformat:
- erste Spalte: `patient_id`
- restliche 96 Spalten: SBS96-Kanaele mit relativen Haeufigkeiten

Die Darstellung ist auf Lesbarkeit fuer groessere Kohorten optimiert:
- breites Querformat
- keine einzelnen Kanalbeschriftungen
- stattdessen nur die sechs Substitutionsgruppen auf der x-Achse
- vertikale Trennlinien zwischen den Gruppen
"""

from __future__ import annotations

import argparse
import os
import re
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-codex")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp/xdg-cache-codex")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

matplotlib.rcParams["font.family"] = "DejaVu Sans"
matplotlib.rcParams["axes.unicode_minus"] = False


SUBSTITUTION_PATTERN = re.compile(r"\[([ACGT]>[ACGT])\]")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Erstellt eine breite SBS96-Heatmap mit Gruppenbeschriftungen."
    )
    parser.add_argument(
        "--input-csv",
        type=Path,
        default=BASE_DIR / "results" / "krebspatienten_sbs96" / "sbs96_per_patient.csv",
        help="Pfad zur SBS96-Matrix im Wide-Format.",
    )
    parser.add_argument(
        "--output-prefix",
        type=Path,
        default=BASE_DIR / "results" / "krebspatienten_sbs96" / "sbs96_heatmap",
        help="Ausgabepfad ohne Dateiendung fuer PNG und SVG.",
    )
    parser.add_argument(
        "--cmap",
        choices=["viridis", "YlOrRd"],
        default="viridis",
        help="Farbskala fuer die Heatmap.",
    )
    return parser.parse_args()


def load_sbs96_matrix(path: Path) -> pd.DataFrame:
    """
    Laedt die SBS96-Tabelle und setzt die erste Spalte als Patientenindex.

    Die Werte werden numerisch erzwungen, damit `imshow` spaeter sicher mit
    einer reinen float-Matrix arbeiten kann.
    """
    dataframe = pd.read_csv(path)
    if dataframe.empty:
        raise ValueError(f"Die Eingabedatei ist leer: {path}")
    if dataframe.shape[1] < 2:
        raise ValueError("Die Datei muss eine ID-Spalte und SBS96-Werte enthalten.")

    index_column = dataframe.columns[0]
    matrix = dataframe.set_index(index_column)
    matrix = matrix.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    return matrix


def get_substitution_blocks(columns: list[str]) -> list[tuple[str, int, int]]:
    """
    Bestimmt zusammenhaengende Bloecke gleicher Substitutionsgruppen.

    Rueckgabeformat:
    - Label der Substitution, z. B. `C>T`
    - Startindex des Blocks (inklusive)
    - Endindex des Blocks (inklusive)
    """
    blocks: list[tuple[str, int, int]] = []
    current_label: str | None = None
    block_start = 0

    for index, column_name in enumerate(columns):
        match = SUBSTITUTION_PATTERN.search(column_name)
        label = match.group(1) if match else column_name
        if current_label is None:
            current_label = label
            block_start = index
            continue
        if label != current_label:
            blocks.append((current_label, block_start, index - 1))
            current_label = label
            block_start = index

    if current_label is not None:
        blocks.append((current_label, block_start, len(columns) - 1))

    return blocks


def plot_heatmap(matrix: pd.DataFrame, output_prefix: Path, cmap: str) -> None:
    blocks = get_substitution_blocks(list(matrix.columns))
    data = matrix.to_numpy(dtype=float)

    fig, ax = plt.subplots(figsize=(24, 10), constrained_layout=True)
    image = ax.imshow(
        data,
        aspect="auto",
        interpolation="nearest",
        cmap=cmap,
    )

    # Patientenachsenlabel bleibt erhalten, aber die einzelnen Namen werden
    # ausgeblendet, weil sie bei grossen Kohorten visuell keinen Mehrwert bieten.
    ax.set_ylabel("Patienten", fontsize=18, fontweight="bold")
    ax.set_yticks([])

    # Statt 96 Einzelkanaelen werden nur die sechs Substitutionsgruppen mittig
    # unter ihren jeweiligen 16 Kanaelen beschriftet.
    tick_positions = [(start + end) / 2 for _, start, end in blocks]
    tick_labels = [label for label, _, _ in blocks]
    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels, fontsize=15, fontweight="bold")
    ax.set_xlabel("Substitutionsgruppen", fontsize=18, fontweight="bold")
    ax.set_title("SBS96-Mutationsprofile aller Patienten", fontsize=24, weight="bold", pad=16)

    # Trennlinien helfen, die sechs Substitutionsgruppen klar voneinander abzugrenzen.
    for _, _, end in blocks[:-1]:
        ax.axvline(end + 0.5, color="white", linewidth=1.2, alpha=0.9)

    colorbar = fig.colorbar(image, ax=ax, pad=0.02)
    colorbar.set_label("Relative Häufigkeit", fontsize=16, fontweight="bold")
    colorbar.ax.tick_params(labelsize=12)

    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_prefix.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(output_prefix.with_suffix(".svg"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    matrix = load_sbs96_matrix(args.input_csv)
    plot_heatmap(matrix=matrix, output_prefix=args.output_prefix, cmap=args.cmap)
    print(f"Fertig. Heatmap gespeichert unter: {args.output_prefix}")


if __name__ == "__main__":
    main()
