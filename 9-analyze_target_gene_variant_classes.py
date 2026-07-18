#!/usr/bin/env python3
"""
Analysiert Variant_Classifications in ausgewaehlten Melanom-Genen aus TCGA-MAF-Dateien.

Das Skript ist fuer einen rekursiven Ordner mit `.maf` / `.maf.gz` Dateien
ausgelegt und beantwortet gezielt die Frage, welche Mutationsarten in einer
kleinen, biologisch relevanten Zielgenliste auftreten.

Analyseschritte:
1. rekursives Einlesen aller MAF-Dateien
2. Filtern auf eine Zielgenliste ueber `Hugo_Symbol`
3. Extraktion ausgewaehlter MAF-Felder pro Mutation
4. Berechnung der absoluten und relativen Haeufigkeiten je Variant_Classification
5. Berechnung der Anzahl betroffener Patienten pro Gen
6. heuristische Einordnung in funktionell relevante vs. nicht funktionell relevante Mutationen
7. Export tabellarischer Zusammenfassungen und zweier Abbildungen

Ausgaben:
- `filtered_target_gene_mutations.csv`
- `variant_classification_per_gene.csv`
- `functional_impact_summary.csv`
- `top_mutations_per_gene.csv`
- `variant_classification_stacked.png` / `.svg`
- `functional_relevance_fraction.png` / `.svg`
- `report.txt`
"""

from __future__ import annotations

import argparse
import csv
import gzip
import os
from collections import Counter
from pathlib import Path
from typing import Iterable, Iterator, TextIO

from script_common import configure_runtime_env, detect_base_dir

BASE_DIR = detect_base_dir(__file__)
configure_runtime_env("bachelorarbeit")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
import numpy as np
import pandas as pd

matplotlib.rcParams["font.family"] = "DejaVu Sans"
matplotlib.rcParams["axes.unicode_minus"] = False


DEFAULT_TARGET_GENES: list[str] = [
    "BRAF",
    "NRAS",
    "TP53",
    "CDKN2A",
    "NF1",
    "PTEN",
    "RB1",
    "RAC1",
    "MAP2K1",
    "MAP2K2",
]

# Binaere Heuristik wie im Prompt gewuenscht: proteinveraendernde bzw.
# klar splice-/truncation-nahe Klassen gelten als funktionell relevant.
FUNCTIONALLY_RELEVANT_CLASSES: set[str] = {
    "Missense_Mutation",
    "Nonsense_Mutation",
    "Frame_Shift_Del",
    "Frame_Shift_Ins",
    "Splice_Site",
    "In_Frame_Del",
    "In_Frame_Ins",
    "Nonstop_Mutation",
    "Translation_Start_Site",
    "Start_Codon_SNP",
    "Start_Codon_Del",
    "Start_Codon_Ins",
    "De_novo_Start_InFrame",
    "De_novo_Start_OutOfFrame",
    "Stop_Codon_Del",
    "Stop_Codon_Ins",
}

NON_FUNCTIONAL_LABEL = "not_functional_relevant"
FUNCTIONAL_LABEL = "functional_relevant"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Analysiert Variant_Classifications und funktionelle Relevanz "
            "in ausgewaehlten Zielgenen aus rekursiv abgelegten MAF-Dateien."
        )
    )
    parser.add_argument(
        "--maf-root",
        type=Path,
        default=BASE_DIR / "data" / "TCGA-SKCM_open_maf" / "files",
        help="Wurzelordner mit rekursiv abgelegten `.maf` oder `.maf.gz` Dateien.",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=BASE_DIR / "results" / "target_gene_variant_analysis",
        help="Ausgabeordner fuer Tabellen, Abbildungen und Report.",
    )
    parser.add_argument(
        "--genes",
        type=str,
        default=",".join(DEFAULT_TARGET_GENES),
        help="Komma-separierte Liste der Zielgene.",
    )
    parser.add_argument(
        "--top-n-mutations",
        type=int,
        default=10,
        help="Wie viele Top-Proteinveraenderungen pro Gen exportiert werden sollen.",
    )
    return parser.parse_args()


def parse_gene_list(raw_genes: str) -> list[str]:
    genes = [(gene or "").strip().upper() for gene in raw_genes.split(",")]
    deduplicated: list[str] = []
    seen: set[str] = set()
    for gene in genes:
        if gene and gene not in seen:
            deduplicated.append(gene)
            seen.add(gene)
    if not deduplicated:
        raise ValueError("Die Zielgenliste ist leer.")
    return deduplicated


def iter_maf_files(maf_root: Path) -> list[Path]:
    if maf_root.is_file():
        return [maf_root]
    if not maf_root.exists():
        raise FileNotFoundError(f"MAF-Pfad nicht gefunden: {maf_root}")

    maf_files = sorted(maf_root.rglob("*.maf")) + sorted(maf_root.rglob("*.maf.gz"))
    if not maf_files:
        raise FileNotFoundError(f"Keine `.maf` oder `.maf.gz` Dateien unter {maf_root} gefunden.")
    return maf_files


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open("rt")


def iter_maf_rows(maf_path: Path) -> Iterator[dict[str, str]]:
    with open_text(maf_path) as handle:
        data_lines: Iterable[str] = (line for line in handle if not line.startswith("#"))
        reader = csv.DictReader(data_lines, delimiter="\t")
        yield from reader


def infer_patient_id(case_id: str, sample_barcode: str) -> str:
    """
    Nutzt bevorzugt `case_id`; bei TCGA-Barcodes faellt sonst auf den
    Patientenanteil des Sample-Barcodes zurueck.
    """
    case_id = case_id.strip()
    sample_barcode = sample_barcode.strip()
    if case_id:
        return case_id
    if sample_barcode.startswith("TCGA-") and len(sample_barcode) >= 12:
        return sample_barcode[:12]
    if sample_barcode:
        return sample_barcode
    return "UNKNOWN_PATIENT"


def classify_functional_relevance(variant_classification: str) -> str:
    """
    Binaere Heuristik gemaess Prompt:
    proteinveraendernde / truncierende / splice-site-nahe Klassen gelten als relevant,
    alle uebrigen Klassen werden als nicht funktionell relevant behandelt.
    """
    if variant_classification in FUNCTIONALLY_RELEVANT_CLASSES:
        return FUNCTIONAL_LABEL
    return NON_FUNCTIONAL_LABEL


def normalize_field(value: str | None, fallback: str = "NA") -> str:
    normalized = (value or "").strip()
    return normalized if normalized else fallback


def extract_target_gene_mutations(
    maf_root: Path,
    target_genes: list[str],
) -> tuple[pd.DataFrame, list[Path], int]:
    maf_files = iter_maf_files(maf_root)
    target_gene_set = set(target_genes)
    records: list[dict[str, str]] = []
    total_rows = 0

    for maf_path in maf_files:
        for row in iter_maf_rows(maf_path):
            total_rows += 1
            gene = normalize_field(row.get("Hugo_Symbol"), fallback="").upper()
            if gene not in target_gene_set:
                continue

            sample_barcode = normalize_field(row.get("Tumor_Sample_Barcode"))
            case_id = normalize_field(row.get("case_id"), fallback="")
            patient_id = infer_patient_id(case_id, sample_barcode)
            variant_classification = normalize_field(row.get("Variant_Classification"))
            records.append(
                {
                    "gene": gene,
                    "variant_classification": variant_classification,
                    "variant_type": normalize_field(row.get("Variant_Type")),
                    "sample_id": sample_barcode,
                    "patient_id": patient_id,
                    "case_id": case_id or patient_id,
                    "hgvsp_short": normalize_field(row.get("HGVSp_Short")),
                    "t_alt_count": normalize_field(row.get("t_alt_count"), fallback=""),
                    "t_ref_count": normalize_field(row.get("t_ref_count"), fallback=""),
                    "functional_relevance": classify_functional_relevance(variant_classification),
                    "maf_file": maf_path.name,
                }
            )

    columns = [
        "gene",
        "variant_classification",
        "variant_type",
        "sample_id",
        "patient_id",
        "case_id",
        "hgvsp_short",
        "t_alt_count",
        "t_ref_count",
        "functional_relevance",
        "maf_file",
    ]
    mutations = pd.DataFrame.from_records(records, columns=columns)
    return mutations, maf_files, total_rows


def build_variant_classification_table(
    mutations: pd.DataFrame,
    target_genes: list[str],
) -> pd.DataFrame:
    columns = [
        "gene",
        "variant_classification",
        "mutation_count",
        "mutation_fraction_within_gene",
        "patient_count",
        "sample_count",
    ]
    if mutations.empty:
        return pd.DataFrame(columns=columns)

    counts = (
        mutations.groupby(["gene", "variant_classification"], as_index=False)
        .agg(
            mutation_count=("gene", "size"),
            patient_count=("patient_id", "nunique"),
            sample_count=("sample_id", "nunique"),
        )
    )
    totals = mutations.groupby("gene").size().rename("total_mutations")
    counts["mutation_fraction_within_gene"] = counts.apply(
        lambda row: row["mutation_count"] / float(totals.loc[row["gene"]]),
        axis=1,
    )

    gene_order = {gene: index for index, gene in enumerate(target_genes)}
    counts["gene_order"] = counts["gene"].map(gene_order)
    class_order = (
        counts.groupby("variant_classification")["mutation_count"]
        .sum()
        .sort_values(ascending=False)
        .index
    )
    class_rank = {name: index for index, name in enumerate(class_order)}
    counts["class_order"] = counts["variant_classification"].map(class_rank)
    counts = counts.sort_values(["gene_order", "class_order", "variant_classification"]).drop(
        columns=["gene_order", "class_order"]
    )
    return counts.reset_index(drop=True)[columns]


def build_functional_impact_summary(
    mutations: pd.DataFrame,
    target_genes: list[str],
) -> pd.DataFrame:
    rows: list[dict[str, float | int | str]] = []
    for gene in target_genes:
        gene_mutations = mutations.loc[mutations["gene"] == gene]
        total_mutations = int(len(gene_mutations))
        affected_patients = int(gene_mutations["patient_id"].nunique()) if total_mutations else 0
        affected_samples = int(gene_mutations["sample_id"].nunique()) if total_mutations else 0
        functional_mutations = int((gene_mutations["functional_relevance"] == FUNCTIONAL_LABEL).sum())
        non_functional_mutations = total_mutations - functional_mutations
        functional_fraction = (
            float(functional_mutations / total_mutations) if total_mutations else 0.0
        )
        rows.append(
            {
                "gene": gene,
                "total_mutations": total_mutations,
                "affected_patients": affected_patients,
                "affected_samples": affected_samples,
                "functional_relevant_mutations": functional_mutations,
                "not_functional_relevant_mutations": non_functional_mutations,
                "functional_relevant_fraction": functional_fraction,
                "not_functional_relevant_fraction": 1.0 - functional_fraction if total_mutations else 0.0,
            }
        )
    return pd.DataFrame(rows)


def top_variant_classification(series: pd.Series) -> str:
    counts = series.value_counts()
    return str(counts.index[0]) if not counts.empty else "NA"


def build_top_mutations_table(
    mutations: pd.DataFrame,
    target_genes: list[str],
    top_n: int,
) -> pd.DataFrame:
    columns = [
        "gene",
        "rank_within_gene",
        "hgvsp_short",
        "variant_classification",
        "mutation_count",
        "patient_count",
        "sample_count",
    ]
    if mutations.empty:
        return pd.DataFrame(columns=columns)

    grouped = (
        mutations.assign(
            hgvsp_short=mutations["hgvsp_short"].replace({"": "protein_change_not_annotated"})
        )
        .groupby(["gene", "hgvsp_short"], as_index=False)
        .agg(
            mutation_count=("gene", "size"),
            patient_count=("patient_id", "nunique"),
            sample_count=("sample_id", "nunique"),
            variant_classification=("variant_classification", top_variant_classification),
        )
    )

    gene_order = {gene: index for index, gene in enumerate(target_genes)}
    grouped["gene_order"] = grouped["gene"].map(gene_order)
    grouped = grouped.sort_values(
        ["gene_order", "mutation_count", "patient_count", "sample_count", "hgvsp_short"],
        ascending=[True, False, False, False, True],
    )
    grouped["rank_within_gene"] = grouped.groupby("gene").cumcount() + 1
    grouped = grouped.loc[grouped["rank_within_gene"] <= top_n]
    grouped = grouped.drop(columns=["gene_order"])
    return grouped.reset_index(drop=True)[columns]


def write_dataframe(dataframe: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    dataframe.to_csv(path, index=False)


def plot_variant_classification_stacked(
    variant_table: pd.DataFrame,
    target_genes: list[str],
    out_prefix: Path,
) -> None:
    if variant_table.empty:
        return

    pivot = (
        variant_table.pivot(
            index="gene",
            columns="variant_classification",
            values="mutation_fraction_within_gene",
        )
        .fillna(0.0)
        .reindex(target_genes)
        .fillna(0.0)
    )
    pivot = pivot.loc[pivot.sum(axis=1) > 0]
    if pivot.empty:
        return

    class_order = (
        variant_table.groupby("variant_classification")["mutation_count"]
        .sum()
        .sort_values(ascending=False)
        .index
        .tolist()
    )
    colors = plt.cm.tab20(np.linspace(0, 1, max(len(class_order), 1)))

    fig, ax = plt.subplots(figsize=(14, 7), constrained_layout=True)
    bottom = np.zeros(len(pivot), dtype=float)
    x = np.arange(len(pivot))

    for color, variant_classification in zip(colors, class_order):
        values = pivot[variant_classification].to_numpy()
        ax.bar(
            x,
            values,
            bottom=bottom,
            label=variant_classification,
            color=color,
            edgecolor="white",
            linewidth=0.5,
        )
        bottom += values

    ax.set_xticks(x)
    ax.set_xticklabels(pivot.index, rotation=35, ha="right")
    ax.set_ylim(0, 1)
    ax.yaxis.set_major_formatter(PercentFormatter(xmax=1, decimals=0))
    ax.set_ylabel("Anteil der Variant_Classification innerhalb des Gens")
    ax.set_xlabel("Gen")
    ax.set_title("Variant_Classifications pro Zielgen", fontsize=15, weight="bold")
    ax.grid(axis="y", linestyle=":", linewidth=0.7, alpha=0.5)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, 1.18),
        ncol=3,
        frameon=False,
        fontsize=9,
    )

    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_prefix.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(out_prefix.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)


def plot_functional_relevance_fraction(
    functional_summary: pd.DataFrame,
    out_prefix: Path,
) -> None:
    summary = functional_summary.loc[functional_summary["total_mutations"] > 0].copy()
    if summary.empty:
        return

    fig, ax = plt.subplots(figsize=(11, 6), constrained_layout=True)
    x = np.arange(len(summary))
    values = summary["functional_relevant_fraction"].to_numpy(dtype=float)
    bars = ax.bar(x, values, color="#d95f02", edgecolor="white", linewidth=0.8)

    ax.set_xticks(x)
    ax.set_xticklabels(summary["gene"], rotation=35, ha="right")
    ax.set_ylim(0, 1.05)
    ax.yaxis.set_major_formatter(PercentFormatter(xmax=1, decimals=0))
    ax.set_ylabel("Anteil funktionell relevanter Mutationen")
    ax.set_xlabel("Gen")
    ax.set_title("Funktionell relevante Mutationen pro Zielgen", fontsize=15, weight="bold")
    ax.grid(axis="y", linestyle=":", linewidth=0.7, alpha=0.5)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    for bar, row in zip(bars, summary.itertuples()):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 0.02,
            f"{row.functional_relevant_mutations}/{row.total_mutations}",
            ha="center",
            va="bottom",
            fontsize=9,
        )

    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_prefix.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(out_prefix.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)


def format_top_items(counter: Counter[str], top_n: int = 3) -> str:
    if not counter:
        return "keine"
    return "; ".join(f"{name}:{count}" for name, count in counter.most_common(top_n))


def write_report(
    path: Path,
    maf_files: list[Path],
    total_rows: int,
    target_genes: list[str],
    mutations: pd.DataFrame,
    functional_summary: pd.DataFrame,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)

    lines = [
        "Zielgenanalyse fuer UV-relevante Melanomgene",
        "===========================================",
        "",
        f"Verarbeitete MAF-Dateien: {len(maf_files)}",
        f"Alle eingelesenen MAF-Zeilen: {total_rows}",
        f"Mutationen in Zielgenen: {len(mutations)}",
        f"Zielgene: {', '.join(target_genes)}",
        "",
        "Definition funktionell relevant:",
        "- relevant: Missense, Nonsense, Frameshift, Splice_Site sowie weitere klar proteinveraendernde Klassen",
        "- nicht relevant: alle uebrigen Variant_Classifications im Sinn dieser binaren Heuristik",
        "",
    ]

    if mutations.empty:
        lines.append("Es wurden keine Mutationen in den ausgewaehlten Zielgenen gefunden.")
        path.write_text("\n".join(lines))
        return

    lines.append("Genweise Zusammenfassung:")
    for row in functional_summary.itertuples():
        gene_mutations = mutations.loc[mutations["gene"] == row.gene]
        class_counter = Counter(gene_mutations["variant_classification"])
        top_hgvsp_counter = Counter(gene_mutations["hgvsp_short"].replace({"NA": "protein_change_not_annotated"}))
        lines.extend(
            [
                f"- {row.gene}: {row.total_mutations} Mutationen, "
                f"{row.affected_patients} Patienten, "
                f"{row.functional_relevant_fraction:.2%} funktionell relevant",
                f"  Top Variant_Classifications: {format_top_items(class_counter)}",
                f"  Top Proteinwechsel: {format_top_items(top_hgvsp_counter)}",
            ]
        )

    path.write_text("\n".join(lines))


def main() -> None:
    args = parse_args()
    target_genes = parse_gene_list(args.genes)
    args.outdir.mkdir(parents=True, exist_ok=True)

    mutations, maf_files, total_rows = extract_target_gene_mutations(args.maf_root, target_genes)

    if not mutations.empty:
        mutations = mutations.sort_values(
            ["gene", "patient_id", "sample_id", "variant_classification", "hgvsp_short"],
            ascending=[True, True, True, True, True],
        ).reset_index(drop=True)

    variant_table = build_variant_classification_table(mutations, target_genes)
    functional_summary = build_functional_impact_summary(mutations, target_genes)
    top_mutations = build_top_mutations_table(mutations, target_genes, args.top_n_mutations)

    write_dataframe(mutations, args.outdir / "filtered_target_gene_mutations.csv")
    write_dataframe(variant_table, args.outdir / "variant_classification_per_gene.csv")
    write_dataframe(functional_summary, args.outdir / "functional_impact_summary.csv")
    write_dataframe(top_mutations, args.outdir / "top_mutations_per_gene.csv")

    plot_variant_classification_stacked(
        variant_table,
        target_genes=target_genes,
        out_prefix=args.outdir / "variant_classification_stacked",
    )
    plot_functional_relevance_fraction(
        functional_summary,
        out_prefix=args.outdir / "functional_relevance_fraction",
    )
    write_report(
        args.outdir / "report.txt",
        maf_files=maf_files,
        total_rows=total_rows,
        target_genes=target_genes,
        mutations=mutations,
        functional_summary=functional_summary,
    )

    print(f"Fertig. Ergebnisse unter: {args.outdir}")
    print(f"Verarbeitete MAF-Dateien: {len(maf_files)}")
    print(f"Mutationen in Zielgenen: {len(mutations)}")


if __name__ == "__main__":
    main()
