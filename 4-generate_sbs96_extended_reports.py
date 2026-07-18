#!/usr/bin/env python3
"""
Erzeugt zusätzliche SBS96-Reports und Grafiken im Stil der alten Auswertung.

Eingaben:
- `merged_snvs_with_context.csv` aus `build_krebspatienten_sbs96_dataset.py`
- optional die bestehende `sbs96_per_patient.csv`

Ausgaben unter `results/neu`:
- `sbs96_aggregated.csv`, `sbs96_aggregated.svg`, `sbs96_summary.txt`
- `sbs96_per_sample.csv`, `sbs96_sample_summary.csv`, `all_samples_overview.svg`, `all_samples_report.txt`
- `all_<N>/...` mit patientenbasierten Zusammenfassungen, Clustern und Grafiken
- `krebspatienten_sbs96/sbs96_heatmap.png` und `sbs96_heatmap.svg`
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

from script_common import configure_runtime_env, detect_base_dir

BASE_DIR = detect_base_dir(__file__)
configure_runtime_env("bachelorarbeit")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")
os.environ.setdefault("LOKY_MAX_CPU_COUNT", "1")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

matplotlib.rcParams["font.family"] = "DejaVu Sans"
matplotlib.rcParams["axes.unicode_minus"] = False


SBS_SUBSTITUTIONS = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
CLUSTER_COLORS = ["#f4a261", "#2a9d8f", "#6d597a", "#e76f51", "#457b9d", "#8d99ae", "#e9c46a", "#90be6d"]
SUBSTITUTION_COLORS = {
    "C>A": "#4C78A8",
    "C>G": "#F58518",
    "C>T": "#E45756",
    "T>A": "#72B7B2",
    "T>C": "#54A24B",
    "T>G": "#EECA3B",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Erzeugt zusätzliche SBS96-Visualisierungen und Reports."
    )
    parser.add_argument(
        "--merged-snvs",
        type=Path,
        default=BASE_DIR / "results" / "krebspatienten_sbs96" / "merged_snvs_with_context.csv",
        help="CSV mit allen SNVs inklusive case_id, patient_id und sbs96_channel.",
    )
    parser.add_argument(
        "--sbs96-relative",
        type=Path,
        default=BASE_DIR / "results" / "krebspatienten_sbs96" / "sbs96_per_patient.csv",
        help="Relative SBS96-Profile pro Sample/Patient für die Heatmap.",
    )
    parser.add_argument(
        "--results-root",
        type=Path,
        default=BASE_DIR / "results",
        help="Wurzelordner für die Zusatzreports.",
    )
    parser.add_argument(
        "--top-n-overview",
        type=int,
        default=25,
        help="Wie viele Top-Samples im Textreport hervorgehoben werden sollen.",
    )
    parser.add_argument(
        "--max-table-rows",
        type=int,
        default=150,
        help="Maximale Zeilenzahl in den großen Tabellen-Grafiken für bessere Lesbarkeit.",
    )
    return parser.parse_args()


def load_data(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    if df.empty:
        raise ValueError(f"Keine Daten in {path} gefunden.")
    return df


def summarize_entity(entity_counts: pd.DataFrame, entity_name: str) -> pd.DataFrame:
    """
    Erstellt Zusammenfassungen mit Gesamt-SBS, Substitutionsfraktionen und Top-Kanälen.
    """
    pivot = entity_counts.pivot(index=entity_name, columns="sbs96_channel", values="count").fillna(0).astype(int)

    substitution_map = pd.Series(
        {column: column.split("[", 1)[1].split("]")[0] for column in pivot.columns},
        name="substitution",
    )
    substitution_counts = {
        substitution: pivot.loc[:, substitution_map[substitution_map == substitution].index].sum(axis=1)
        for substitution in SBS_SUBSTITUTIONS
    }
    total_sbs = pivot.sum(axis=1)

    summary = pd.DataFrame(
        {
            entity_name: pivot.index,
            "total_sbs": total_sbs.values,
        }
    )
    for substitution in SBS_SUBSTITUTIONS:
        counts = substitution_counts[substitution]
        summary[f"{substitution}_count"] = counts.values
        summary[f"{substitution}_fraction"] = np.where(total_sbs.values > 0, counts.values / total_sbs.values, 0.0)

    top_channels = pivot.apply(lambda row: row.sort_values(ascending=False).head(3), axis=1)
    for rank in range(3):
        summary[f"top{rank+1}_channel"] = top_channels.apply(lambda s, r=rank: s.index[r] if len(s) > r else "", axis=1).values
        summary[f"top{rank+1}_count"] = top_channels.apply(
            lambda s, r=rank: int(s.iloc[r]) if len(s) > r and pd.notna(s.iloc[r]) else 0,
            axis=1,
        ).values

    return summary.sort_values("total_sbs", ascending=False).reset_index(drop=True), pivot


def aggregated_sbs96(entity_counts: pd.DataFrame) -> pd.DataFrame:
    agg = entity_counts.groupby("sbs96_channel", as_index=False)["count"].sum().sort_values("count", ascending=False)
    total = agg["count"].sum()
    agg["fraction"] = np.where(total > 0, agg["count"] / total, 0.0)
    return agg


def write_dataframe(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


def pairwise_euclidean_distance(matrix: np.ndarray) -> np.ndarray:
    """Berechnet eine vollständige euklidische Distanzmatrix."""
    squared_norms = np.sum(matrix * matrix, axis=1, keepdims=True)
    distances_sq = squared_norms + squared_norms.T - 2.0 * matrix @ matrix.T
    distances_sq = np.maximum(distances_sq, 0.0)
    return np.sqrt(distances_sq)


def run_kmeans_numpy(matrix: np.ndarray, k: int, random_state: int = 42, max_iter: int = 100) -> tuple[np.ndarray, np.ndarray, float]:
    """Einfache K-Means-Implementierung ohne sklearn."""
    rng = np.random.default_rng(random_state)
    centroid_indices = rng.choice(len(matrix), size=k, replace=False)
    centroids = matrix[centroid_indices].copy()
    labels = np.zeros(len(matrix), dtype=int)

    for _ in range(max_iter):
        distances = np.sum((matrix[:, None, :] - centroids[None, :, :]) ** 2, axis=2)
        new_labels = np.argmin(distances, axis=1)
        if np.array_equal(new_labels, labels):
            break
        labels = new_labels
        for cluster_id in range(k):
            mask = labels == cluster_id
            if np.any(mask):
                centroids[cluster_id] = matrix[mask].mean(axis=0)
            else:
                centroids[cluster_id] = matrix[rng.integers(0, len(matrix))]

    inertia = float(np.sum((matrix - centroids[labels]) ** 2))
    return labels, centroids, inertia


def compute_pca_projection(matrix: np.ndarray, n_components: int = 2) -> tuple[np.ndarray, np.ndarray]:
    """
    Berechnet eine einfache PCA-Projektion via SVD.

    Die SBS96-Fraktionen werden dabei lediglich zentriert, aber nicht weiter
    skaliert, damit die bereits normierten Profile in ihrer relativen Struktur
    erhalten bleiben.
    """
    if matrix.ndim != 2 or matrix.shape[0] == 0:
        return np.zeros((0, n_components)), np.zeros(n_components)

    centered = matrix - matrix.mean(axis=0, keepdims=True)
    _, singular_values, right_vectors = np.linalg.svd(centered, full_matrices=False)
    projection = centered @ right_vectors[:n_components].T

    if len(matrix) <= 1:
        return projection, np.zeros(n_components)

    variances = (singular_values ** 2) / (len(matrix) - 1)
    total_variance = float(variances.sum())
    if total_variance <= 0:
        explained = np.zeros(n_components)
    else:
        explained = variances[:n_components] / total_variance
    return projection, explained


def ordered_sbs96_columns(columns: list[str]) -> list[str]:
    """Ordnet SBS96-Kanäle nach den sechs Substitutionsgruppen."""
    ordered_columns: list[str] = []
    for substitution in SBS_SUBSTITUTIONS:
        sub_cols = sorted(column for column in columns if f"[{substitution}]" in column)
        ordered_columns.extend(sub_cols)
    return ordered_columns


def silhouette_score_numpy(matrix: np.ndarray, labels: np.ndarray) -> float:
    """Berechnet den Silhouette-Score aus einer Distanzmatrix."""
    unique_labels = np.unique(labels)
    if len(unique_labels) < 2:
        return 0.0

    distances = pairwise_euclidean_distance(matrix)
    scores: list[float] = []
    for index in range(len(matrix)):
        same_cluster = labels == labels[index]
        same_cluster[index] = False

        if np.any(same_cluster):
            a = float(distances[index, same_cluster].mean())
        else:
            a = 0.0

        b_candidates = []
        for other_label in unique_labels:
            if other_label == labels[index]:
                continue
            mask = labels == other_label
            if np.any(mask):
                b_candidates.append(float(distances[index, mask].mean()))
        b = min(b_candidates) if b_candidates else 0.0

        denom = max(a, b)
        scores.append((b - a) / denom if denom > 0 else 0.0)
    return float(np.mean(scores))


def plot_cluster_pca(
    fractions: pd.DataFrame,
    labels: np.ndarray,
    total_sbs: np.ndarray,
    path: Path,
) -> None:
    """Projiziert die vollständigen SBS96-Profile per PCA in zwei Dimensionen."""
    matrix = fractions.to_numpy(dtype=float, copy=False)
    coords, explained = compute_pca_projection(matrix, n_components=2)
    fig, ax = plt.subplots(figsize=(9.5, 7.4), constrained_layout=True)
    sizes = np.clip(np.sqrt(total_sbs) * 1.2, 16, 110)

    for cluster_id in sorted(np.unique(labels)):
        mask = labels == cluster_id
        color = CLUSTER_COLORS[(cluster_id - 1) % len(CLUSTER_COLORS)]
        ax.scatter(
            coords[mask, 0],
            coords[mask, 1],
            s=sizes[mask],
            color=color,
            alpha=0.78,
            edgecolors="#333333",
            linewidths=0.35,
            label=f"Cluster {cluster_id}",
        )
        centroid = coords[mask].mean(axis=0)
        ax.scatter(
            centroid[0],
            centroid[1],
            s=220,
            color=color,
            marker="X",
            edgecolors="white",
            linewidths=1.0,
            zorder=5,
        )

    ax.set_title("PCA-Projektion der SBS96-Cluster", fontsize=16, weight="bold")
    ax.set_xlabel(f"PC1 ({explained[0] * 100:.1f}% erklärte Varianz)")
    ax.set_ylabel(f"PC2 ({explained[1] * 100:.1f}% erklärte Varianz)")
    ax.grid(True, linestyle=":", alpha=0.45)
    ax.legend(frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(path.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_cluster_mean_profiles(fractions: pd.DataFrame, labels: np.ndarray, path: Path) -> None:
    """Zeigt das mittlere SBS96-Profil je Cluster als gruppierte Balkendarstellung."""
    ordered_columns = ordered_sbs96_columns(list(fractions.columns))
    group_centers: list[float] = []
    group_labels: list[str] = []
    group_boundaries: list[float] = []
    next_start = 0
    for substitution in SBS_SUBSTITUTIONS:
        sub_cols = [column for column in ordered_columns if f"[{substitution}]" in column]
        if not sub_cols:
            continue
        start = next_start
        end = start + len(sub_cols) - 1
        group_centers.append((start + end) / 2)
        group_labels.append(substitution)
        next_start = end + 1
        group_boundaries.append(end + 0.5)

    fractions = fractions[ordered_columns]
    columns = ordered_columns
    x_positions = np.arange(len(columns))
    substitution_labels = [column.split("[", 1)[1].split("]")[0] for column in columns]
    bar_colors = [SUBSTITUTION_COLORS[substitution] for substitution in substitution_labels]

    unique_labels = sorted(np.unique(labels))
    fig, axes = plt.subplots(
        len(unique_labels),
        1,
        figsize=(18, max(4.6, 2.9 * len(unique_labels) + 1.2)),
        sharex=True,
        constrained_layout=False,
    )
    if len(unique_labels) == 1:
        axes = [axes]

    ymax = 0.0
    cluster_means: dict[int, pd.Series] = {}
    for cluster_id in unique_labels:
        mask = labels == cluster_id
        cluster_mean = fractions.loc[mask].mean(axis=0)
        cluster_means[cluster_id] = cluster_mean
        ymax = max(ymax, float(cluster_mean.max()))

    for axis, cluster_id in zip(axes, unique_labels):
        cluster_mean = cluster_means[cluster_id]
        axis.bar(x_positions, cluster_mean.to_numpy(), color=bar_colors, width=0.92, edgecolor="none")
        for boundary in group_boundaries[:-1]:
            axis.axvline(boundary, color="#FFFFFF", linewidth=1.1, alpha=0.95)
        axis.set_ylim(0, ymax * 1.12 if ymax > 0 else 1.0)
        axis.set_ylabel("Mittlerer Anteil")
        axis.set_title(
            f"Cluster {cluster_id} (n={int((labels == cluster_id).sum())})",
            fontsize=13,
            weight="bold",
            loc="left",
        )
        axis.grid(axis="y", linestyle=":", linewidth=0.7, alpha=0.45)
        axis.set_axisbelow(True)
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)

    axes[-1].set_xticks(group_centers)
    axes[-1].set_xticklabels(group_labels, fontsize=11, fontweight="bold")
    axes[-1].set_xlabel("Substitutionsgruppen")

    handles = [
        plt.Rectangle((0, 0), 1, 1, color=SUBSTITUTION_COLORS[substitution], label=substitution)
        for substitution in SBS_SUBSTITUTIONS
    ]
    fig.subplots_adjust(top=0.88, hspace=0.15)
    fig.suptitle("Mittlere SBS96-Profile pro Cluster", fontsize=16, weight="bold", y=0.985)
    fig.legend(
        handles=handles,
        loc="upper center",
        ncol=len(SBS_SUBSTITUTIONS),
        frameon=False,
        bbox_to_anchor=(0.5, 0.955),
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(path.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def write_cluster_mean_profiles(fractions: pd.DataFrame, labels: np.ndarray, path: Path) -> None:
    """Schreibt die mittleren SBS96-Profile pro Cluster als lange Tabelle."""
    ordered_columns = ordered_sbs96_columns(list(fractions.columns))
    fractions = fractions[ordered_columns]

    rows: list[dict[str, object]] = []
    for cluster_id in sorted(np.unique(labels)):
        mask = labels == cluster_id
        cluster_mean = fractions.loc[mask].mean(axis=0)
        cluster_size = int(mask.sum())
        for channel_order, channel in enumerate(ordered_columns, start=1):
            substitution = channel.split("[", 1)[1].split("]")[0]
            rows.append(
                {
                    "cluster": int(cluster_id),
                    "cluster_size": cluster_size,
                    "channel_order": channel_order,
                    "sbs96_channel": channel,
                    "substitution": substitution,
                    "mean_fraction": float(cluster_mean[channel]),
                }
            )
    write_dataframe(pd.DataFrame(rows), path)


def plot_silhouette_model_selection(model_records: list[tuple[int, float, float]], best_k: int, path: Path) -> None:
    """Visualisiert den Silhouette-Score für die getesteten Clusterzahlen."""
    ks = [record[0] for record in model_records]
    scores = [record[1] for record in model_records]

    fig, ax = plt.subplots(figsize=(8.2, 5.4), constrained_layout=True)
    ax.plot(ks, scores, color="#2A9D8F", marker="o", linewidth=2.2, markersize=7)
    for k, score in zip(ks, scores):
        ax.text(k, score + 0.012, f"{score:.3f}", ha="center", va="bottom", fontsize=9)
    ax.axvline(best_k, color="#E76F51", linestyle="--", linewidth=1.3, alpha=0.85)
    ax.scatter([best_k], [scores[ks.index(best_k)]], color="#E76F51", s=80, zorder=5, label=f"gewähltes k = {best_k}")
    ax.set_title("Silhouetten-Score für k = 2 bis 8", fontsize=15, weight="bold")
    ax.set_xlabel("Anzahl Cluster (k)")
    ax.set_ylabel("Silhouetten-Score")
    ax.set_xticks(ks)
    ax.grid(axis="y", linestyle=":", linewidth=0.8, alpha=0.45)
    ax.legend(frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(path.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_aggregated_spectrum(agg: pd.DataFrame, path: Path) -> None:
    subs = agg.assign(substitution=agg["sbs96_channel"].str.extract(r"\[(\w>\w)\]")[0]).groupby("substitution", as_index=False)["count"].sum()
    subs["fraction"] = subs["count"] / subs["count"].sum()
    subs = subs.set_index("substitution").reindex(SBS_SUBSTITUTIONS).reset_index()

    fig, ax = plt.subplots(figsize=(9, 5), constrained_layout=True)
    colors = plt.cm.viridis(np.linspace(0.15, 0.9, len(subs)))
    bars = ax.bar(subs["substitution"], subs["fraction"] * 100, color=colors, edgecolor="white")
    ax.set_title("Aggregiertes Substitutionsspektrum", fontsize=15, weight="bold")
    ax.set_xlabel("Substitution")
    ax.set_ylabel("Anteil (%)")
    ax.set_axisbelow(True)
    ax.grid(axis="y", linestyle=":", linewidth=0.8, alpha=0.5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for bar, value in zip(bars, subs["fraction"] * 100):
        ax.text(bar.get_x() + bar.get_width() / 2, value + 0.3, f"{value:.2f}%", ha="center", va="bottom", fontsize=9)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(path.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_heatmap(sbs96_relative_path: Path, outdir: Path) -> None:
    df = pd.read_csv(sbs96_relative_path, index_col=0)
    df = df.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    # Nach C>T-Dominanz sortieren, damit die Heatmap lesbarer wird.
    c_to_t_cols = [col for col in df.columns if "[C>T]" in col]
    if c_to_t_cols:
        df = df.assign(_sort=df[c_to_t_cols].sum(axis=1)).sort_values("_sort", ascending=False).drop(columns="_sort")

    fig, ax = plt.subplots(figsize=(24, 10), constrained_layout=True)
    sns.heatmap(df, cmap="viridis", cbar_kws={"label": "Relative Häufigkeit"}, xticklabels=False, yticklabels=False, ax=ax)
    ax.set_title("SBS96-Heatmap aller Profile", fontsize=24, weight="bold", pad=16)
    ax.set_xlabel("Substitutionsgruppen", fontsize=18, fontweight="bold")
    ax.set_ylabel("Patienten", fontsize=18, fontweight="bold")
    ax.set_xticks([7.5, 23.5, 39.5, 55.5, 71.5, 87.5])
    ax.set_xticklabels(SBS_SUBSTITUTIONS, fontsize=15, fontweight="bold")
    for boundary in [15.5, 31.5, 47.5, 63.5, 79.5]:
        ax.axvline(boundary, color="white", linewidth=1.2, alpha=0.9)
    colorbar = ax.collections[0].colorbar
    if colorbar is not None:
        colorbar.set_label("Relative Häufigkeit", fontsize=16, fontweight="bold")
        colorbar.ax.tick_params(labelsize=12)
    outdir.mkdir(parents=True, exist_ok=True)
    fig.savefig(outdir / "sbs96_heatmap.svg", bbox_inches="tight")
    fig.savefig(outdir / "sbs96_heatmap.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def make_table_plot(
    summary: pd.DataFrame,
    id_col: str,
    title: str,
    value_cols: list[str],
    label_map: dict[str, str],
    path: Path,
    extra_col: str | None = None,
    max_rows: int | None = None,
) -> None:
    """
    Zeichnet eine kompakte Heatmap-Tabelle ähnlich den alten SVG-Übersichten.
    """
    table = summary[[id_col] + ([extra_col] if extra_col else []) + value_cols].copy()
    if max_rows is not None and len(table) > max_rows:
        table = table.head(max_rows).copy()
    rows = len(table)
    fig_height = max(7, rows * 0.22 + 1.8)
    fig, ax = plt.subplots(figsize=(16, fig_height))

    numeric_part = table.drop(columns=[id_col] + ([extra_col] if extra_col else []))
    sns.heatmap(
        numeric_part,
        cmap="coolwarm_r",
        cbar=True,
        annot=False,
        yticklabels=False,
        xticklabels=[label_map.get(col, col) for col in numeric_part.columns],
        ax=ax,
    )
    ax.set_title(title, fontsize=16, weight="bold", pad=12)
    ax.set_xlabel("")
    ax.set_ylabel("")

    # Textspalte links außerhalb der Heatmap ergänzen.
    for row_index, row in table.iterrows():
        ax.text(-0.4, row_index + 0.5, str(row[id_col]), ha="right", va="center", fontsize=8, transform=ax.transData)
        if extra_col:
            ax.text(-0.15, row_index + 0.5, str(row[extra_col]), ha="center", va="center", fontsize=8, transform=ax.transData)

    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)


def choose_k(profile_matrix: pd.DataFrame) -> tuple[int, list[tuple[int, float, float]]]:
    records: list[tuple[int, float, float]] = []
    best_k = 2
    best_score = -1.0
    matrix = profile_matrix.to_numpy(dtype=float, copy=False)
    max_k = min(8, len(matrix) - 1)
    for k in range(2, max_k + 1):
        labels, _, inertia = run_kmeans_numpy(matrix, k=k, random_state=42)
        score = silhouette_score_numpy(matrix, labels)
        records.append((k, score, inertia))
        if score > best_score:
            best_score = score
            best_k = k
    return best_k, records


def run_clustering(
    patient_counts: pd.DataFrame,
    patient_summary: pd.DataFrame,
    outdir: Path,
    max_table_rows: int | None = None,
) -> None:
    fractions = patient_counts.div(patient_counts.sum(axis=1), axis=0).fillna(0.0)
    matrix = fractions.to_numpy(dtype=float, copy=False)
    best_k, model_records = choose_k(fractions)
    raw_labels, _, _ = run_kmeans_numpy(matrix, k=best_k, random_state=42)
    labels = raw_labels + 1
    silhouette = silhouette_score_numpy(matrix, raw_labels)
    summary_by_case = patient_summary.set_index("case_id").loc[fractions.index].copy()
    cluster_series = pd.Series(labels, index=fractions.index, name="cluster")

    clusters = summary_by_case.join(cluster_series, how="inner").reset_index()
    clusters = clusters.rename(columns={"case_id": "patient"})
    clusters = clusters[["patient", "cluster", "total_sbs"]].sort_values(["cluster", "total_sbs"], ascending=[True, False])
    write_dataframe(clusters, outdir / "patient_clusters.csv")

    cluster_lines = [
        "K-Means Clustering der Patienten auf Basis der SBS96-Fraktionen",
        "===========================================================",
        "",
        f"Patienten mit Profil: {len(fractions)}",
        f"Gewaehlte Clusterzahl k: {best_k}",
        f"Silhouette-Score: {silhouette:.4f}",
        "",
        "Modellvergleich:",
    ]
    for k, score, inertia in model_records:
        cluster_lines.append(f"k={k}\tSilhouette={score:.4f}\tInertia={inertia:.6f}")

    for cluster_id in sorted(clusters["cluster"].unique()):
        patient_ids = clusters.loc[clusters["cluster"] == cluster_id, "patient"]
        cluster_counts = patient_counts.loc[patient_ids]
        cluster_fracs = cluster_counts.div(cluster_counts.sum(axis=1), axis=0).fillna(0.0)
        cluster_totals = clusters.loc[clusters["cluster"] == cluster_id, "total_sbs"]
        top_channels = cluster_fracs.mean(axis=0).sort_values(ascending=False).head(8)
        top_patients = clusters.loc[clusters["cluster"] == cluster_id].head(10)
        cluster_lines.extend(
            [
                "",
                f"Cluster {cluster_id}",
                f"Groesse: {len(patient_ids)} Patienten",
                f"Mittelwert SBS/Patient: {cluster_totals.mean():.2f}" if len(cluster_totals) else "",
                f"Median SBS/Patient: {cluster_totals.median():.2f}" if len(cluster_totals) else "",
                "Top Kanaele: " + ", ".join(f"{channel} ({value:.2%})" for channel, value in top_channels.items()),
                "Top Patienten nach SBS: " + ", ".join(f"{row.patient} ({int(row.total_sbs)})" for row in top_patients.itertuples()),
            ]
        )
    write_text(outdir / "cluster_summary.txt", "\n".join(line for line in cluster_lines if line is not None))

    # Scatter analog zur alten Darstellung: C>T vs T>C, Punktgröße ~ Gesamt-SBS
    sample_fracs = summary_by_case
    fig, ax = plt.subplots(figsize=(9, 7), constrained_layout=True)
    c_t = sample_fracs["C>T_fraction"].to_numpy()
    t_c = sample_fracs["T>C_fraction"].to_numpy()
    sizes = np.clip(np.sqrt(sample_fracs["total_sbs"].to_numpy()) * 1.2, 12, 90)
    cluster_values = cluster_series.loc[sample_fracs.index].to_numpy()
    for cluster_id in sorted(np.unique(cluster_values)):
        mask = cluster_values == cluster_id
        ax.scatter(
            c_t[mask],
            t_c[mask],
            s=sizes[mask],
            color=CLUSTER_COLORS[(cluster_id - 1) % len(CLUSTER_COLORS)],
            alpha=0.75,
            edgecolors="#333333",
            linewidths=0.4,
            label=f"Cluster {cluster_id}",
        )
    ax.set_title("Clusterprojektion (C>T vs T>C)", fontsize=16, weight="bold")
    ax.set_xlabel("C>T-Anteil")
    ax.set_ylabel("T>C-Anteil")
    ax.grid(True, linestyle=":", alpha=0.5)
    ax.legend(frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.savefig(outdir / "patient_clusters_scatter.svg", bbox_inches="tight")
    fig.savefig(outdir / "patient_clusters_scatter.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    plot_cluster_pca(
        fractions=fractions,
        labels=labels,
        total_sbs=sample_fracs["total_sbs"].to_numpy(dtype=float),
        path=outdir / "patient_clusters_pca",
    )
    plot_cluster_mean_profiles(fractions=fractions, labels=labels, path=outdir / "patient_cluster_mean_profiles")
    write_cluster_mean_profiles(fractions=fractions, labels=labels, path=outdir / "cluster_mean_sbs96_profiles.csv")
    plot_silhouette_model_selection(model_records=model_records, best_k=best_k, path=outdir / "patient_clusters_silhouette")

    cluster_overview = summary_by_case.join(cluster_series, how="inner").reset_index()
    cluster_overview = cluster_overview.sort_values(["cluster", "total_sbs"], ascending=[True, False])
    make_table_plot(
        summary=cluster_overview,
        id_col="case_id",
        title="Patientencluster nach SBS96-Profil",
        value_cols=["total_sbs"] + [f"{sub}_fraction" for sub in SBS_SUBSTITUTIONS],
        label_map={
            "total_sbs": "SBS gesamt",
            "C>A_fraction": "C>A",
            "C>G_fraction": "C>G",
            "C>T_fraction": "C>T",
            "T>A_fraction": "T>A",
            "T>C_fraction": "T>C",
            "T>G_fraction": "T>G",
        },
        path=outdir / "patient_clusters_overview",
        extra_col="cluster",
        max_rows=max_table_rows,
    )


def main() -> None:
    args = parse_args()
    print("Lade SNV-Datensatz...", flush=True)
    merged = load_data(args.merged_snvs)
    merged["count"] = 1

    root = args.results_root
    root.mkdir(parents=True, exist_ok=True)

    # Gesamt- und Sample-Ebene
    print("Berechne Sample-Zusammenfassungen...", flush=True)
    sample_counts = merged.groupby(["patient_id", "sbs96_channel"], as_index=False)["count"].sum()
    sample_summary, sample_pivot = summarize_entity(sample_counts, "patient_id")
    sample_summary = sample_summary.rename(columns={"patient_id": "sample"})
    sample_pivot = sample_pivot.sort_index()
    agg = aggregated_sbs96(sample_counts)

    write_dataframe(agg, root / "sbs96_aggregated.csv")
    plot_aggregated_spectrum(agg, root / "sbs96_aggregated")

    sample_wide = sample_pivot.reset_index().rename(columns={"patient_id": "sample"})
    sample_wide.insert(1, "total_sbs", sample_pivot.sum(axis=1).astype(int).values)
    write_dataframe(sample_wide.rename(columns={"sample": "id"}), root / "sbs96_per_sample.csv")
    write_dataframe(sample_summary, root / "sbs96_sample_summary.csv")

    subtype_lines = []
    total_agg = agg["count"].sum()
    for substitution in SBS_SUBSTITUTIONS:
        sub_count = agg.loc[agg["sbs96_channel"].str.contains(f"\\[{substitution}\\]"), "count"].sum()
        subtype_lines.append(f"{substitution}\t{sub_count}\t{sub_count / total_agg:.2%}")
    write_text(root / "subtype_summary.txt", "\n".join(subtype_lines))

    top_samples = sample_summary.head(args.top_n_overview)
    sample_report_lines = [
        "SBS96 Analyse aller lokal verfuegbaren Proben",
        "===========================================",
        "",
        f"Proben: {len(sample_summary)}",
        f"Gesamte SBS: {int(sample_summary['total_sbs'].sum())}",
        f"Mittelwert SBS/Probe: {sample_summary['total_sbs'].mean():.2f}",
        f"Median SBS/Probe: {sample_summary['total_sbs'].median():.2f}",
        "",
        "Aggregierte Hauptsubtypen:",
        *subtype_lines,
        "",
        f"Top {len(top_samples)} Proben nach SBS:",
    ]
    for row in top_samples.itertuples():
        sample_report_lines.append(
            f"{row.sample}\t{int(row.total_sbs)}\t"
            f"{[(row.top1_channel, row.top1_count), (row.top2_channel, row.top2_count), (row.top3_channel, row.top3_count)]}"
        )
    write_text(root / "all_samples_report.txt", "\n".join(sample_report_lines))
    write_text(root / "sbs96_summary.txt", "\n".join(sample_report_lines[:12]))

    print("Erzeuge Sample-Übersichten...", flush=True)
    make_table_plot(
        summary=sample_summary,
        id_col="sample",
        title="Alle Proben: SBS96-Übersicht",
        value_cols=["total_sbs"] + [f"{sub}_fraction" for sub in SBS_SUBSTITUTIONS],
        label_map={
            "total_sbs": "SBS gesamt",
            "C>A_fraction": "C>A",
            "C>G_fraction": "C>G",
            "C>T_fraction": "C>T",
            "T>A_fraction": "T>A",
            "T>C_fraction": "T>C",
            "T>G_fraction": "T>G",
        },
        path=root / "all_samples_overview",
        max_rows=args.max_table_rows,
    )

    # Patienten-/Case-Ebene
    print("Berechne Patienten-Zusammenfassungen...", flush=True)
    patient_counts = merged.groupby(["case_id", "sbs96_channel"], as_index=False)["count"].sum()
    patient_summary, patient_pivot = summarize_entity(patient_counts, "case_id")
    patient_summary = patient_summary.rename(columns={"case_id": "patient"})
    patient_wide = patient_pivot.reset_index().rename(columns={"case_id": "patient"})
    patient_wide.insert(1, "total_sbs", patient_pivot.sum(axis=1).astype(int).values)

    maf_count = merged["maf_file"].nunique()
    patient_outdir = root / f"all_{maf_count}"
    patient_outdir.mkdir(parents=True, exist_ok=True)

    write_dataframe(agg, patient_outdir / "sbs96_aggregated.csv")
    write_dataframe(sample_wide.rename(columns={"sample": "id"}), patient_outdir / "sbs96_per_sample.csv")
    write_dataframe(patient_wide, patient_outdir / "sbs96_per_patient.csv")
    write_dataframe(patient_summary, patient_outdir / "patient_summary.csv")
    write_text(patient_outdir / "subtype_summary.txt", "\n".join(subtype_lines))

    report_lines = [
        "SBS96 Analyse aller lokal verfuegbaren Patienten",
        "=============================================",
        "",
        f"MAF-Dateien: {maf_count}",
        f"Proben: {len(sample_summary)}",
        f"Patienten: {len(patient_summary)}",
        f"Gesamte SBS: {int(patient_summary['total_sbs'].sum())}",
        f"Mittelwert SBS/Patient: {patient_summary['total_sbs'].mean():.2f}",
        f"Median SBS/Patient: {patient_summary['total_sbs'].median():.2f}",
        f"Min SBS/Patient: {int(patient_summary['total_sbs'].min())}",
        f"Max SBS/Patient: {int(patient_summary['total_sbs'].max())}",
        "",
        "Aggregierte Hauptsubtypen:",
        *subtype_lines,
        "",
        f"Top {min(20, len(patient_summary))} Patienten nach SBS:",
    ]
    for row in patient_summary.head(20).itertuples():
        report_lines.append(
            f"{row.patient}\t{int(row.total_sbs)}\t"
            f"{[(row.top1_channel, row.top1_count), (row.top2_channel, row.top2_count), (row.top3_channel, row.top3_count)]}"
        )
    write_text(patient_outdir / "report.txt", "\n".join(report_lines))

    print("Erzeuge Patienten-Übersichten und Cluster...", flush=True)
    make_table_plot(
        summary=patient_summary,
        id_col="patient",
        title="Patientenübersicht SBS96",
        value_cols=["total_sbs"] + [f"{sub}_fraction" for sub in SBS_SUBSTITUTIONS],
        label_map={
            "total_sbs": "SBS gesamt",
            "C>A_fraction": "C>A",
            "C>G_fraction": "C>G",
            "C>T_fraction": "C>T",
            "T>A_fraction": "T>A",
            "T>C_fraction": "T>C",
            "T>G_fraction": "T>G",
        },
        path=patient_outdir / "patient_overview",
        max_rows=args.max_table_rows,
    )

    run_clustering(
        patient_pivot,
        patient_summary.rename(columns={"patient": "case_id"}),
        patient_outdir,
        max_table_rows=args.max_table_rows,
    )

    print("Erzeuge Heatmap...", flush=True)
    plot_heatmap(args.sbs96_relative, root / "krebspatienten_sbs96")

    print(f"Fertig. Zusatzreports unter: {root}")
    print(f"Patientenordner: {patient_outdir}")


if __name__ == "__main__":
    main()
