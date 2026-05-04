#!/usr/bin/env python3
"""
Fuehrt MAF-Dateien aus einem Patientenordner zusammen und berechnet SBS96-Profile.

Das Skript erledigt folgende Schritte:
1. rekursives Einlesen aller `.maf` / `.maf.gz` Dateien aus einem Eingabeordner
2. Filtern auf echte SNVs (hier im GDC-MAF als `Variant_Type == SNP`)
3. Bestimmen des Nukleotid-Kontexts aus dem Referenzgenom hg38
4. Ueberfuehren aller gueltigen SNVs in einen gemeinsamen Datensatz
5. Berechnen eines 96er Mutationsprofils (SBS96) pro Patient
6. Schaetzen der Anteile ausgewaehlter COSMIC-SBS-Signaturen pro Patient

Die Ausgabe besteht aus:
- `merged_snvs_with_context.csv`: alle gueltigen SNVs aller Patienten in einer Tabelle
- `sbs96_per_patient.csv`: 96er Profil pro Patient im Wide-Format
- `signature_shares_per_patient.csv`: geschaetzte Anteile fuer SBS1, SBS2, SBS5, SBS7, SBS13 und SBS38
- `signature_share_summary.csv`: Mittelwert und Median der Signaturanteile
- `signature_share_cohort_mean.svg`: eigenes Diagramm fuer das Kohortenmittel
- `signature_shares_stacked.svg`: Abbildung der Signaturanteile pro Patient
- `summary.txt`: kompakte Zusammenfassung des Laufs

Beispiel:
    python scripts/build_krebspatienten_sbs96_dataset.py
"""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import re
import urllib.request
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, TextIO

BASE_DIR = Path(__file__).resolve().parent.parent

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-codex")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp/xdg-cache-codex")

import matplotlib
import numpy as np
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
from scipy.optimize import nnls

matplotlib.rcParams["font.family"] = "DejaVu Sans"
matplotlib.rcParams["axes.unicode_minus"] = False


SBS_SUBSTITUTIONS: list[str] = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
SBS96_CHANNELS: list[str] = [
    f"{left}[{sub}]{right}"
    for sub in SBS_SUBSTITUTIONS
    for left in "ACGT"
    for right in "ACGT"
]
SBS96_SET: set[str] = set(SBS96_CHANNELS)
COMPLEMENT = str.maketrans("ACGT", "TGCA")
COSMIC_CANDIDATE_URLS: list[str] = [
    "https://cancer.sanger.ac.uk/signatures/documents/2124/COSMIC_v3.4_SBS_GRCh38.txt",
    "https://cancer.sanger.ac.uk/signatures/documents/2047/COSMIC_v3.3.1_SBS_GRCh38.txt",
    "https://cancer.sanger.ac.uk/signatures/documents/453/COSMIC_v3.2_SBS_GRCh38.txt",
]
REPORTED_SIGNATURES: list[str] = ["SBS1", "SBS2", "SBS5", "SBS13", "SBS38"]
AGGREGATED_SIGNATURES: dict[str, list[str]] = {
    # In aktuellen COSMIC-Versionen ist SBS7 in mehrere UV-Subsignaturen aufgeteilt.
    "SBS7": ["SBS7a", "SBS7b", "SBS7c", "SBS7d"],
}
PLOT_SIGNATURE_ORDER: list[str] = ["SBS1", "SBS2", "SBS5", "SBS7", "SBS13", "SBS38", "other_signatures"]
PLOT_COLORS: dict[str, str] = {
    "SBS1": "#4C78A8",
    "SBS2": "#F58518",
    "SBS5": "#54A24B",
    "SBS7": "#E45756",
    "SBS13": "#72B7B2",
    "SBS38": "#EECA3B",
    "other_signatures": "#B8B8B8",
}
SUMMARY_SIGNATURES: list[str] = ["SBS1", "SBS2", "SBS5", "SBS7", "SBS13", "SBS38"]


@dataclass(frozen=True)
class MergedSnvRecord:
    """Repraesentiert eine einzelne, gueltige SNV nach allen Filterschritten."""

    patient_id: str
    case_id: str
    patient_folder: str
    maf_file: str
    hugo_symbol: str
    variant_classification: str
    chromosome: str
    start_position: int
    end_position: int
    ref_allele_maf: str
    alt_allele_maf: str
    ref_allele_hg38: str
    hg38_trinucleotide_context: str
    maf_context_11mer: str
    sbs_substitution: str
    sbs96_channel: str
    ref_matches_hg38: bool


def reverse_complement(sequence: str) -> str:
    """Bildet das Reverse-Complement einer DNA-Sequenz."""
    return sequence.translate(COMPLEMENT)[::-1]


def is_valid_base(base: str) -> bool:
    """Nur A, C, G, T sind fuer die SBS96-Klassifikation gueltig."""
    return len(base) == 1 and base in "ACGT"


def normalize_sbs96(ref: str, alt: str, trinucleotide_context: str) -> tuple[str, str] | tuple[None, None]:
    """
    Ueberfuehrt eine SNV in die kanonische pyrimidin-zentrierte SBS96-Schreibweise.

    Beispiel:
    - G>A in GCA wird gespiegelt zu C>T in TGC
    - Rueckgabe: ("T[C>T]G", "C>T")
    """
    ref = ref.upper()
    alt = alt.upper()
    trinucleotide_context = trinucleotide_context.upper()

    if len(trinucleotide_context) != 3 or any(base not in "ACGT" for base in trinucleotide_context):
        return None, None
    if not (is_valid_base(ref) and is_valid_base(alt)) or ref == alt:
        return None, None
    if trinucleotide_context[1] != ref:
        return None, None

    # SBS96 wird immer aus Sicht eines Pyrimidins (C oder T) beschrieben.
    if ref in "AG":
        trinucleotide_context = reverse_complement(trinucleotide_context)
        ref = trinucleotide_context[1]
        alt = alt.translate(COMPLEMENT)

    if ref not in "CT":
        return None, None

    substitution = f"{ref}>{alt}"
    channel = f"{trinucleotide_context[0]}[{substitution}]{trinucleotide_context[2]}"
    if channel not in SBS96_SET:
        return None, None

    return channel, substitution


class IndexedFasta:
    """
    Minimaler FASTA-Reader mit Random-Access ueber eine `.fai`-Datei.

    Dadurch kann fuer jede Mutation gezielt die Referenzbase bzw. der
    Trinukleotid-Kontext abgefragt werden, ohne die gesamte FASTA in den
    Arbeitsspeicher zu laden.
    """

    def __init__(self, fasta_path: Path, fai_path: Path | None = None) -> None:
        self.fasta_path = fasta_path
        self.fai_path = fai_path or Path(f"{fasta_path}.fai")
        self.index = self._read_fai(self.fai_path)
        self.handle = self.fasta_path.open("rb")
        self.alias_map = self._build_alias_map()

    def close(self) -> None:
        self.handle.close()

    @staticmethod
    def _read_fai(path: Path) -> dict[str, tuple[int, int, int, int]]:
        if not path.exists():
            raise FileNotFoundError(f"FASTA-Index (.fai) nicht gefunden: {path}")

        index: dict[str, tuple[int, int, int, int]] = {}
        with path.open() as handle:
            for line in handle:
                if not line.strip():
                    continue
                name, length, offset, line_bases, line_width = line.rstrip("\n").split("\t")[:5]
                index[name] = (int(length), int(offset), int(line_bases), int(line_width))
        return index

    def _build_alias_map(self) -> dict[str, str]:
        """
        Erlaubt sowohl RefSeq-Accession-Namen als auch `chr1`, `chr2`, ..., `chrX`.
        """
        aliases: dict[str, str] = {}
        refseq_map = {
            "1": "NC_000001.11",
            "2": "NC_000002.12",
            "3": "NC_000003.12",
            "4": "NC_000004.12",
            "5": "NC_000005.10",
            "6": "NC_000006.12",
            "7": "NC_000007.14",
            "8": "NC_000008.11",
            "9": "NC_000009.12",
            "10": "NC_000010.11",
            "11": "NC_000011.10",
            "12": "NC_000012.12",
            "13": "NC_000013.11",
            "14": "NC_000014.9",
            "15": "NC_000015.10",
            "16": "NC_000016.10",
            "17": "NC_000017.11",
            "18": "NC_000018.10",
            "19": "NC_000019.10",
            "20": "NC_000020.11",
            "21": "NC_000021.9",
            "22": "NC_000022.11",
            "X": "NC_000023.11",
            "Y": "NC_000024.10",
            "M": "NC_012920.1",
            "MT": "NC_012920.1",
        }

        for contig_name in self.index:
            aliases[contig_name] = contig_name
            aliases[contig_name.upper()] = contig_name

        for chrom, accession in refseq_map.items():
            if accession in self.index:
                aliases[chrom] = accession
                aliases[f"chr{chrom}"] = accession
                aliases[f"CHR{chrom}"] = accession

        return aliases

    def resolve_chrom(self, chrom: str) -> str | None:
        if not chrom:
            return None
        return self.alias_map.get(chrom) or self.alias_map.get(chrom.upper())

    def fetch_base(self, chrom: str, pos1: int) -> str | None:
        """
        Liest exakt eine Base anhand der `.fai`-Offsets aus der FASTA.

        `pos1` ist 1-basiert, wie in MAF-Dateien ueblich.
        """
        resolved = self.resolve_chrom(chrom)
        if resolved is None:
            return None

        length, offset, line_bases, line_width = self.index[resolved]
        if pos1 < 1 or pos1 > length:
            return None

        pos0 = pos1 - 1
        byte_offset = offset + (pos0 // line_bases) * line_width + (pos0 % line_bases)
        self.handle.seek(byte_offset)
        base = self.handle.read(1).decode("ascii", errors="ignore").upper()
        return base if base in "ACGTN" else None

    def fetch_trinucleotide(self, chrom: str, center_pos1: int) -> str | None:
        left = self.fetch_base(chrom, center_pos1 - 1)
        center = self.fetch_base(chrom, center_pos1)
        right = self.fetch_base(chrom, center_pos1 + 1)
        if not left or not center or not right:
            return None
        trinucleotide = f"{left}{center}{right}"
        return trinucleotide if re.fullmatch(r"[ACGT]{3}", trinucleotide) else None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Liest MAF-Dateien ein, filtert auf SNVs, bestimmt den hg38-Kontext "
            "und erstellt ein SBS96-Profil pro Patient."
        )
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=BASE_DIR / "data" / "TCGA-SKCM_open_maf" / "files",
        help="Ordner mit Patienten-Unterordnern und .maf/.maf.gz-Dateien.",
    )
    parser.add_argument(
        "--hg38-fasta",
        type=Path,
        default=BASE_DIR / "GCF_000001405.40_GRCh38.p14_genomic.fna",
        help="Pfad zur hg38-FASTA-Datei.",
    )
    parser.add_argument(
        "--hg38-fai",
        type=Path,
        default=BASE_DIR / "GCF_000001405.40_GRCh38.p14_genomic.fna.fai",
        help="Pfad zur zugehoerigen hg38-.fai-Datei.",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=BASE_DIR / "results" / "krebspatienten_sbs96",
        help="Ausgabeordner fuer zusammengefuehrte Daten und SBS96-Profile.",
    )
    parser.add_argument(
        "--cosmic-file",
        type=Path,
        default=BASE_DIR / "reference" / "COSMIC_v3.4_SBS_GRCh38.txt",
        help="Lokaler Pfad zur COSMIC-SBS96-Referenzmatrix (GRCh38).",
    )
    parser.add_argument(
        "--download-cosmic-if-missing",
        action="store_true",
        help="Lade die COSMIC-SBS96-Referenzmatrix herunter, falls sie lokal fehlt.",
    )
    return parser.parse_args()


def iter_maf_files(input_dir: Path) -> list[Path]:
    """Sammelt alle MAF-Dateien rekursiv ein."""
    if input_dir.is_file():
        return [input_dir]
    if not input_dir.exists():
        raise FileNotFoundError(f"Eingabeordner nicht gefunden: {input_dir}")

    maf_files = sorted(input_dir.rglob("*.maf")) + sorted(input_dir.rglob("*.maf.gz"))
    if not maf_files:
        raise FileNotFoundError(f"Keine .maf/.maf.gz Dateien unter {input_dir} gefunden.")
    return maf_files


def detect_delimiter(header_line: str) -> str:
    """Unterstuetzt tab- und komma-separierte Signaturmatrizen."""
    return "\t" if header_line.count("\t") >= header_line.count(",") else ","


def maybe_download_cosmic(target_path: Path) -> bool:
    """Versucht mehrere bekannte offizielle COSMIC-URLs, bis eine erreichbar ist."""
    target_path.parent.mkdir(parents=True, exist_ok=True)
    for url in COSMIC_CANDIDATE_URLS:
        try:
            with urllib.request.urlopen(url, timeout=60) as response:
                target_path.write_bytes(response.read())
            return True
        except Exception:
            continue
    return False


def load_cosmic_signatures(cosmic_path: Path) -> dict[str, dict[str, float]]:
    """
    Laedt eine COSMIC-SBS96-Matrix.

    Erwartet werden 96 Kanaele in Zeilen und Signaturen in Spalten.
    Die Werte in der Matrix sind bereits relative Haeufigkeiten pro Signatur.
    """
    if not cosmic_path.exists():
        raise FileNotFoundError(f"COSMIC-Datei nicht gefunden: {cosmic_path}")

    with cosmic_path.open() as handle:
        first_line = handle.readline()
    delimiter = detect_delimiter(first_line)

    signatures: dict[str, dict[str, float]] = {}
    with cosmic_path.open() as handle:
        reader = csv.reader(handle, delimiter=delimiter)
        header = next(reader)
        signature_names = [cell.strip() for cell in header[1:] if cell.strip()]
        for signature_name in signature_names:
            signatures[signature_name] = {channel: 0.0 for channel in SBS96_CHANNELS}

        for row in reader:
            if not row:
                continue
            channel = row[0].strip()
            if channel not in SBS96_SET:
                continue
            for index, signature_name in enumerate(signature_names, start=1):
                try:
                    value = float(row[index]) if index < len(row) and row[index] else 0.0
                except ValueError:
                    value = 0.0
                signatures[signature_name][channel] = value

    return signatures


def open_text(path: Path) -> TextIO:
    """Oeffnet normale Textdateien und gzip-komprimierte MAF-Dateien transparent."""
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open("rt")


def choose_alt_allele(row: dict[str, str]) -> str:
    """
    Waehlt das mutierte Tumor-Allel.

    In MAF-Dateien ist `Tumor_Seq_Allele2` meist das relevante Alternativallel.
    Falls dieses aber leer ist oder der Referenz entspricht, nutzen wir
    `Tumor_Seq_Allele1` als Rueckfall.
    """
    ref = (row.get("Reference_Allele") or "").strip().upper()
    allele2 = (row.get("Tumor_Seq_Allele2") or "").strip().upper()
    allele1 = (row.get("Tumor_Seq_Allele1") or "").strip().upper()

    if is_valid_base(allele2) and allele2 != ref:
        return allele2
    if is_valid_base(allele1) and allele1 != ref:
        return allele1
    return allele2


def iter_maf_rows(maf_path: Path) -> Iterator[dict[str, str]]:
    """Liefert nur Datenzeilen, Kommentarzeilen (`#`) werden ignoriert."""
    with open_text(maf_path) as handle:
        data_lines: Iterable[str] = (line for line in handle if not line.startswith("#"))
        reader = csv.DictReader(data_lines, delimiter="\t")
        yield from reader


def record_from_row(row: dict[str, str], maf_path: Path, fasta: IndexedFasta) -> MergedSnvRecord | None:
    """
    Baut aus einer MAF-Zeile einen SNV-Datensatz inklusive hg38-Kontext.

    Ungueltige, nicht-kanonisierbare oder nicht-SNV-Ereignisse liefern `None`.
    """
    if (row.get("Variant_Type") or "").strip() != "SNP":
        return None

    ref_maf = (row.get("Reference_Allele") or "").strip().upper()
    alt = choose_alt_allele(row)
    chrom = (row.get("Chromosome") or "").strip()
    start_raw = (row.get("Start_Position") or "").strip()
    end_raw = (row.get("End_Position") or "").strip()

    if not (is_valid_base(ref_maf) and is_valid_base(alt) and alt != ref_maf):
        return None

    try:
        start_position = int(float(start_raw))
        end_position = int(float(end_raw)) if end_raw else start_position
    except ValueError:
        return None

    ref_hg38 = fasta.fetch_base(chrom, start_position)
    trinucleotide = fasta.fetch_trinucleotide(chrom, start_position)
    if ref_hg38 is None or ref_hg38 == "N" or trinucleotide is None:
        return None

    channel, substitution = normalize_sbs96(ref_hg38, alt, trinucleotide)
    if channel is None or substitution is None:
        return None

    patient_id = (
        (row.get("Tumor_Sample_Barcode") or "").strip()
        or (row.get("case_id") or "").strip()
        or maf_path.parent.name
    )
    case_id = (row.get("case_id") or "").strip() or maf_path.parent.name

    return MergedSnvRecord(
        patient_id=patient_id,
        case_id=case_id,
        patient_folder=maf_path.parent.name,
        maf_file=maf_path.name,
        hugo_symbol=(row.get("Hugo_Symbol") or "").strip(),
        variant_classification=(row.get("Variant_Classification") or "").strip(),
        chromosome=chrom,
        start_position=start_position,
        end_position=end_position,
        ref_allele_maf=ref_maf,
        alt_allele_maf=alt,
        ref_allele_hg38=ref_hg38,
        hg38_trinucleotide_context=trinucleotide,
        maf_context_11mer=(row.get("CONTEXT") or "").strip(),
        sbs_substitution=substitution,
        sbs96_channel=channel,
        ref_matches_hg38=(ref_maf == ref_hg38),
    )


def write_merged_dataset(path: Path, records: list[MergedSnvRecord]) -> None:
    """Schreibt alle gueltigen SNVs in einen gemeinsamen CSV-Datensatz."""
    fieldnames = [
        "patient_id",
        "case_id",
        "patient_folder",
        "maf_file",
        "hugo_symbol",
        "variant_classification",
        "chromosome",
        "start_position",
        "end_position",
        "ref_allele_maf",
        "alt_allele_maf",
        "ref_allele_hg38",
        "hg38_trinucleotide_context",
        "maf_context_11mer",
        "sbs_substitution",
        "sbs96_channel",
        "ref_matches_hg38",
    ]

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for record in records:
            writer.writerow(record.__dict__)


def write_sbs96_per_patient(path: Path, per_patient_profiles: dict[str, Counter[str]]) -> None:
    """
    Schreibt eine SBS96-Matrix im Wide-Format mit relativen Haeufigkeiten.

    Jede Zeile entspricht einem Patienten, jede der 96 Spalten einem SBS96-Kanal.
    Die Werte sind auf Summe 1 normiert, damit Patienten direkt vergleichbar sind.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["patient_id"] + SBS96_CHANNELS)
        for patient_id in sorted(per_patient_profiles):
            counts = per_patient_profiles[patient_id]
            total_snv = sum(counts.values())
            if total_snv == 0:
                frequencies = [0.0 for _ in SBS96_CHANNELS]
            else:
                frequencies = [counts.get(channel, 0) / total_snv for channel in SBS96_CHANNELS]
            writer.writerow([patient_id] + frequencies)


def cosine_similarity(vector_a: np.ndarray, vector_b: np.ndarray) -> float:
    """Misst, wie gut rekonstruiertes und beobachtetes Profil zueinander passen."""
    norm_a = float(np.linalg.norm(vector_a))
    norm_b = float(np.linalg.norm(vector_b))
    if norm_a == 0.0 or norm_b == 0.0:
        return 0.0
    return float(np.dot(vector_a, vector_b) / (norm_a * norm_b))


def fit_signature_shares(
    per_patient_profiles: dict[str, Counter[str]],
    signatures: dict[str, dict[str, float]],
) -> dict[str, dict[str, float]]:
    """
    Schaetzt Signaturanteile per Non-Negative Least Squares (NNLS).

    Gefittet wird gegen alle in der COSMIC-Matrix vorhandenen Signaturen.
    Fuer die Ausgabe werden anschliessend nur die gewuenschten Signaturen
    extrahiert; `SBS7` wird dabei als Summe von `SBS7a-d` berichtet.
    """
    signature_names = list(signatures)
    signature_matrix = np.array(
        [[signatures[name][channel] for name in signature_names] for channel in SBS96_CHANNELS],
        dtype=float,
    )

    patient_results: dict[str, dict[str, float]] = {}
    for patient_id, counts in per_patient_profiles.items():
        total_snv = sum(counts.values())
        observed_profile = np.array(
            [counts.get(channel, 0) / total_snv if total_snv else 0.0 for channel in SBS96_CHANNELS],
            dtype=float,
        )

        if observed_profile.sum() == 0.0:
            result = {signature: 0.0 for signature in REPORTED_SIGNATURES}
            result["SBS7"] = 0.0
            result["other_signatures"] = 0.0
            result["reconstruction_cosine_similarity"] = 0.0
            patient_results[patient_id] = result
            continue

        exposures, _ = nnls(signature_matrix, observed_profile)
        exposure_sum = float(exposures.sum())
        if exposure_sum == 0.0:
            normalized_exposures = np.zeros_like(exposures)
        else:
            normalized_exposures = exposures / exposure_sum

        exposure_map = {
            signature_name: float(normalized_exposures[index])
            for index, signature_name in enumerate(signature_names)
        }

        reconstructed_profile = signature_matrix @ normalized_exposures
        result = {
            signature: exposure_map.get(signature, 0.0)
            for signature in REPORTED_SIGNATURES
        }
        result["SBS7"] = sum(exposure_map.get(signature, 0.0) for signature in AGGREGATED_SIGNATURES["SBS7"])
        result["other_signatures"] = max(0.0, 1.0 - sum(result.values()))
        result["reconstruction_cosine_similarity"] = cosine_similarity(observed_profile, reconstructed_profile)
        patient_results[patient_id] = result

    return patient_results


def write_signature_shares(path: Path, signature_shares: dict[str, dict[str, float]]) -> None:
    """Schreibt die geschaetzten Anteile der gewuenschten Mutationssignaturen."""
    fieldnames = [
        "patient_id",
        "SBS1",
        "SBS2",
        "SBS5",
        "SBS7",
        "SBS13",
        "SBS38",
        "other_signatures",
        "reconstruction_cosine_similarity",
    ]

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for patient_id in sorted(signature_shares):
            writer.writerow({"patient_id": patient_id, **signature_shares[patient_id]})


def compute_signature_statistics(signature_shares: dict[str, dict[str, float]]) -> dict[str, dict[str, float]]:
    """Berechnet Mittelwert und Median der gewuenschten Signaturanteile."""
    statistics: dict[str, dict[str, float]] = {}
    if not signature_shares:
        return statistics

    for signature in SUMMARY_SIGNATURES:
        values = np.array([row[signature] for row in signature_shares.values()], dtype=float)
        statistics[signature] = {
            "mean": float(np.mean(values)),
            "median": float(np.median(values)),
        }
    return statistics


def write_signature_statistics(path: Path, signature_statistics: dict[str, dict[str, float]]) -> None:
    """Schreibt Mittelwert und Median der Signaturanteile in eine kompakte CSV."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["signature", "mean", "median"])
        writer.writeheader()
        for signature in SUMMARY_SIGNATURES:
            writer.writerow(
                {
                    "signature": signature,
                    "mean": signature_statistics[signature]["mean"],
                    "median": signature_statistics[signature]["median"],
                }
            )


def write_signature_share_plot(
    path: Path,
    signature_shares: dict[str, dict[str, float]],
) -> None:
    """
    Erstellt eine gut lesbare Abbildung der Signaturanteile.

    Ein horizontaler, gestapelter Balkenplot ist hier passend, weil viele
    Patienten gleichzeitig dargestellt werden und die Zusammensetzung pro
    Patient direkt vergleichbar bleibt.
    """
    if not signature_shares:
        return

    sorted_rows = sorted(
        signature_shares.items(),
        key=lambda item: (
            item[1]["SBS7"],
            item[1]["SBS38"],
            item[1]["SBS13"],
            item[1]["SBS1"],
        ),
        reverse=True,
    )

    patient_ids = [patient_id for patient_id, _ in sorted_rows]
    y_positions = np.arange(len(patient_ids), dtype=float)
    values_by_signature = {
        signature: np.array([row[signature] for _, row in sorted_rows], dtype=float)
        for signature in PLOT_SIGNATURE_ORDER
    }
    figure_width = 11.2
    figure_height = 7.6
    fig, ax_patients = plt.subplots(
        1,
        1,
        figsize=(figure_width, figure_height),
    )
    fig.subplots_adjust(top=0.88, left=0.10, right=0.98, bottom=0.10)

    cumulative_left = np.zeros(len(patient_ids), dtype=float)
    for signature in PLOT_SIGNATURE_ORDER:
        ax_patients.barh(
            y_positions,
            values_by_signature[signature],
            left=cumulative_left,
            height=0.95,
            color=PLOT_COLORS[signature],
            edgecolor="white",
            linewidth=0.15,
            label=signature,
        )
        cumulative_left += values_by_signature[signature]

    ax_patients.set_xlim(0, 1)
    ax_patients.xaxis.set_major_formatter(PercentFormatter(xmax=1, decimals=0))
    ax_patients.grid(axis="x", linestyle=":", linewidth=0.7, alpha=0.5)
    ax_patients.set_axisbelow(True)
    ax_patients.spines["top"].set_visible(False)
    ax_patients.spines["right"].set_visible(False)
    ax_patients.set_xlabel("Anteil", fontsize=10)
    ax_patients.set_ylabel(f"Patienten (n={len(patient_ids)})", fontsize=10)
    ax_patients.set_yticks([])
    ax_patients.tick_params(axis="y", left=False, labelleft=False)
    ax_patients.tick_params(axis="x", labelsize=9)
    ax_patients.margins(y=0)
    ax_patients.invert_yaxis()
    fig.suptitle("Anteile ausgewählter SBS-Signaturen pro Patient", fontsize=14, weight="bold", y=0.975)

    handles, labels = ax_patients.get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        ncol=len(PLOT_SIGNATURE_ORDER),
        loc="upper center",
        bbox_to_anchor=(0.5, 0.945),
        frameon=False,
        fontsize=8,
    )

    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=300, bbox_inches="tight")
    fig.savefig(path.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def write_cohort_mean_plot(path: Path, signature_shares: dict[str, dict[str, float]]) -> None:
    """Erstellt ein eigenes Diagramm fuer das Kohortenmittel der Signaturanteile."""
    if not signature_shares:
        return

    cohort_mean = {
        signature: float(np.mean([row[signature] for row in signature_shares.values()]))
        for signature in PLOT_SIGNATURE_ORDER
    }

    fig, ax = plt.subplots(1, 1, figsize=(10, 2.8), constrained_layout=True)

    left = 0.0
    for signature in PLOT_SIGNATURE_ORDER:
        ax.barh(
            ["Kohortenmittel"],
            [cohort_mean[signature]],
            left=left,
            height=0.95,
            color=PLOT_COLORS[signature],
            edgecolor="white",
            linewidth=0.7,
            label=signature,
        )
        left += cohort_mean[signature]

    ax.set_xlim(0, 1)
    ax.xaxis.set_major_formatter(PercentFormatter(xmax=1, decimals=0))
    ax.grid(axis="x", linestyle=":", linewidth=0.7, alpha=0.5)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_title("Kohortenmittel der ausgewählten SBS-Signaturen", fontsize=14, weight="bold", pad=10)
    ax.set_xlabel("Anteil")
    ax.tick_params(axis="y", labelsize=10)
    ax.margins(y=0)

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(
        handles,
        labels,
        ncol=len(PLOT_SIGNATURE_ORDER),
        loc="upper center",
        bbox_to_anchor=(0.5, 1.35),
        frameon=False,
        fontsize=9,
    )

    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=300, bbox_inches="tight")
    fig.savefig(path.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def write_summary(
    path: Path,
    maf_files: list[Path],
    total_rows: int,
    snv_rows: int,
    merged_records: list[MergedSnvRecord],
    per_patient_profiles: dict[str, Counter[str]],
    signature_shares: dict[str, dict[str, float]],
    signature_statistics: dict[str, dict[str, float]],
    empty_maf_files: list[Path],
    maf_files_without_valid_snvs: list[Path],
    cosmic_file: Path,
) -> None:
    """Erstellt eine kurze, gut lesbare Textzusammenfassung."""
    substitution_counts = Counter(record.sbs_substitution for record in merged_records)
    hg38_matches = sum(1 for record in merged_records if record.ref_matches_hg38)
    hg38_mismatches = len(merged_records) - hg38_matches

    with path.open("w") as handle:
        handle.write("Krebspatienten SNV/SBS96 Zusammenfassung\n")
        handle.write("=======================================\n\n")
        handle.write(f"Verarbeitete MAF-Dateien: {len(maf_files)}\n")
        handle.write(f"Alle eingelesenen MAF-Zeilen: {total_rows}\n")
        handle.write(f"MAF-Zeilen mit Variant_Type == SNP: {snv_rows}\n")
        handle.write(f"Gueltige SNVs mit hg38-Kontext: {len(merged_records)}\n")
        handle.write(f"Patienten mit SBS96-Profil: {len(per_patient_profiles)}\n")
        handle.write(f"Patienten mit Signatur-Fit: {len(signature_shares)}\n")
        handle.write(f"Leere MAF-Dateien: {len(empty_maf_files)}\n")
        handle.write(f"MAF-Dateien ohne gueltige SNVs: {len(maf_files_without_valid_snvs)}\n")
        handle.write(f"COSMIC-Referenzdatei: {cosmic_file}\n")
        handle.write(f"hg38-Referenz stimmt mit MAF ueberein: {hg38_matches}\n")
        handle.write(f"hg38-Referenz weicht von MAF ab: {hg38_mismatches}\n\n")
        handle.write("Substitutionsspektrum (kanonische Pyrimidin-Sicht):\n")
        total_substitutions = sum(substitution_counts.values())
        for substitution in SBS_SUBSTITUTIONS:
            count = substitution_counts.get(substitution, 0)
            fraction = count / total_substitutions if total_substitutions else 0.0
            handle.write(f"- {substitution}: {count} ({fraction:.2%})\n")
        if signature_statistics:
            handle.write("\nGeschaetzte Signaturanteile ueber alle Patienten:\n")
            for signature in SUMMARY_SIGNATURES:
                handle.write(
                    f"- {signature}: Mittelwert {signature_statistics[signature]['mean']:.2%}, "
                    f"Median {signature_statistics[signature]['median']:.2%}\n"
                )
        if empty_maf_files:
            handle.write("\nLeere MAF-Dateien:\n")
            for maf_path in empty_maf_files:
                handle.write(f"- {maf_path}\n")


def main() -> None:
    args = parse_args()
    maf_files = iter_maf_files(args.input_dir)
    args.outdir.mkdir(parents=True, exist_ok=True)

    merged_records: list[MergedSnvRecord] = []
    per_patient_profiles: dict[str, Counter[str]] = defaultdict(Counter)
    signature_shares: dict[str, dict[str, float]] = {}
    signature_statistics: dict[str, dict[str, float]] = {}
    total_rows = 0
    snv_rows = 0
    empty_maf_files: list[Path] = []
    maf_files_without_valid_snvs: list[Path] = []

    fasta = IndexedFasta(args.hg38_fasta, args.hg38_fai)
    try:
        for maf_path in maf_files:
            rows_in_file = 0
            valid_snvs_before = len(merged_records)
            for row in iter_maf_rows(maf_path):
                rows_in_file += 1
                total_rows += 1
                if (row.get("Variant_Type") or "").strip() == "SNP":
                    snv_rows += 1

                record = record_from_row(row, maf_path, fasta)
                if record is None:
                    continue

                merged_records.append(record)
                per_patient_profiles[record.patient_id][record.sbs96_channel] += 1

            valid_snvs_in_file = len(merged_records) - valid_snvs_before
            if rows_in_file == 0:
                empty_maf_files.append(maf_path)
            elif valid_snvs_in_file == 0:
                maf_files_without_valid_snvs.append(maf_path)
    finally:
        fasta.close()

    if not args.cosmic_file.exists():
        if not args.download_cosmic_if_missing:
            raise FileNotFoundError(
                "Fuer die Berechnung der Signaturanteile fehlt die COSMIC-Referenzdatei. "
                "Lege sie unter "
                f"{args.cosmic_file} ab oder starte das Skript mit --download-cosmic-if-missing."
            )
        if not maybe_download_cosmic(args.cosmic_file):
            raise RuntimeError("Die COSMIC-Referenzdatei konnte nicht heruntergeladen werden.")

    cosmic_signatures = load_cosmic_signatures(args.cosmic_file)
    signature_shares = fit_signature_shares(per_patient_profiles, cosmic_signatures)
    signature_statistics = compute_signature_statistics(signature_shares)

    write_merged_dataset(args.outdir / "merged_snvs_with_context.csv", merged_records)
    write_sbs96_per_patient(args.outdir / "sbs96_per_patient.csv", per_patient_profiles)
    write_signature_shares(args.outdir / "signature_shares_per_patient.csv", signature_shares)
    write_signature_statistics(args.outdir / "signature_share_summary.csv", signature_statistics)
    write_cohort_mean_plot(args.outdir / "signature_share_cohort_mean.svg", signature_shares)
    write_signature_share_plot(args.outdir / "signature_shares_stacked.svg", signature_shares)
    write_summary(
        args.outdir / "summary.txt",
        maf_files=maf_files,
        total_rows=total_rows,
        snv_rows=snv_rows,
        merged_records=merged_records,
        per_patient_profiles=per_patient_profiles,
        signature_shares=signature_shares,
        signature_statistics=signature_statistics,
        empty_maf_files=empty_maf_files,
        maf_files_without_valid_snvs=maf_files_without_valid_snvs,
        cosmic_file=args.cosmic_file,
    )

    print(f"Fertig. Ergebnisse liegen unter: {args.outdir}")
    print(f"Zusammengefuehrte SNVs: {len(merged_records)}")
    print(f"Patienten mit SBS96-Profil: {len(per_patient_profiles)}")
    print(f"Patienten mit Signaturanteilen: {len(signature_shares)}")


if __name__ == "__main__":
    main()
