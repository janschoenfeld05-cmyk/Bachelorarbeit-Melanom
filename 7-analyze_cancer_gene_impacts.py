#!/usr/bin/env python3
"""
Analysiert Mutationen in Protoonkogenen und Tumorsuppressorgenen.

Das Skript kombiniert drei Informationsquellen:
1. MAF-Dateien mit den genomischen Positionen und funktionellen Varianteneffekten
2. das Referenzgenom hg38 zur Validierung der Referenzallele
3. eine kuratierte Genrollen-Liste aus OncoKB

Ziel ist es, fuer alle Mutationen in relevanten Krebsgenen abzuschaetzen,
ob die beobachteten Veraenderungen eher zu einem Protoonkogen-typischen
aktivierenden Muster oder zu einem Tumorsuppressor-typischen
inaktivierenden Muster passen.

Ausgaben:
- `cancer_gene_mutations_detailed.csv`: detaillierte Mutationsliste
- `cancer_gene_summary.csv`: Zusammenfassung pro Gen
- `proto_oncogene_mutations.csv`
- `tumor_suppressor_mutations.csv`
- `context_dependent_gene_mutations.csv`
- `proto_oncogene_summary.csv`
- `tumor_suppressor_summary.csv`
- `context_dependent_gene_summary.csv`
- `report.txt`
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import urllib.request
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, TextIO

BASE_DIR = Path(__file__).resolve().parent.parent


ONCOKB_CANDIDATE_URL = "https://www.oncokb.org/api/v1/utils/cancerGeneList"
ROLE_PROTO = "Proto-Onkogen"
ROLE_TSG = "Tumorsuppressorgen"
ROLE_DUAL = "Kontextabhaengig"

LOSS_OF_FUNCTION_CLASSES = {
    "Frame_Shift_Del",
    "Frame_Shift_Ins",
    "Nonsense_Mutation",
    "Splice_Site",
    "Translation_Start_Site",
    "Nonstop_Mutation",
    "Start_Codon_SNP",
    "Start_Codon_Del",
    "Start_Codon_Ins",
    "De_novo_Start_OutOfFrame",
    "Stop_Codon_Del",
    "Stop_Codon_Ins",
}

ACTIVATING_LIKE_CLASSES = {
    "Missense_Mutation",
    "In_Frame_Del",
    "In_Frame_Ins",
    "De_novo_Start_InFrame",
}

NEUTRAL_OR_NONCODING_CLASSES = {
    "Silent",
    "Intron",
    "IGR",
    "3'UTR",
    "5'UTR",
    "3'Flank",
    "5'Flank",
    "RNA",
    "lincRNA",
    "Targeted_Region",
    "Splice_Region",
}

IMPACT_ORDER = {"HIGH": 3, "MODERATE": 2, "LOW": 1, "MODIFIER": 0, "NA": -1, "": -1}


@dataclass(frozen=True)
class CancerGeneMutation:
    """Repraesentiert eine einzelne Mutation in einem relevanten Krebsgen."""

    gene: str
    gene_role: str
    oncokb_gene_type: str
    patient_id: str
    sample_id: str
    case_id: str
    chromosome: str
    start_position: int
    end_position: int
    variant_type: str
    variant_classification: str
    impact: str
    consequence: str
    one_consequence: str
    hgvsp_short: str
    protein_position: str
    ref_allele_maf: str
    alt_allele_maf: str
    hg38_ref_segment: str
    ref_matches_hg38: str
    t_ref_count: str
    t_alt_count: str
    tumor_vaf: str
    role_consistency: str
    impact_interpretation: str
    maf_file: str


class IndexedFasta:
    """Kleiner Random-Access-Reader fuer FASTA-Dateien mit `.fai`-Index."""

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

        for contig in self.index:
            aliases[contig] = contig
            aliases[contig.upper()] = contig

        for chrom, accession in refseq_map.items():
            if accession in self.index:
                aliases[chrom] = accession
                aliases[f"chr{chrom}"] = accession
                aliases[f"CHR{chrom}"] = accession
        return aliases

    def resolve_chrom(self, chrom: str) -> str | None:
        return self.alias_map.get(chrom) or self.alias_map.get(chrom.upper())

    def fetch_base(self, chrom: str, pos1: int) -> str | None:
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

    def fetch_sequence(self, chrom: str, start1: int, end1: int) -> str | None:
        if end1 < start1:
            return None
        bases: list[str] = []
        for position in range(start1, end1 + 1):
            base = self.fetch_base(chrom, position)
            if base is None:
                return None
            bases.append(base)
        return "".join(bases)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Analysiert Mutationen in Protoonkogenen und Tumorsuppressorgenen "
            "auf Basis von MAF-Dateien, hg38 und einer OncoKB-Genrollenliste."
        )
    )
    parser.add_argument(
        "--maf-root",
        type=Path,
        default=BASE_DIR / "data" / "TCGA-SKCM_open_maf" / "files",
        help="Ordner mit rekursiv abgelegten `.maf` oder `.maf.gz` Dateien.",
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
        help="Pfad zur zugehoerigen `.fai`-Datei.",
    )
    parser.add_argument(
        "--oncokb-gene-list",
        type=Path,
        default=BASE_DIR / "reference" / "oncokb_cancer_gene_list.json",
        help="Lokale JSON-Datei mit der OncoKB-Krebsgene-Liste.",
    )
    parser.add_argument(
        "--download-oncokb-if-missing",
        action="store_true",
        help="Lade die OncoKB-Krebsgene-Liste herunter, falls sie lokal fehlt.",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=BASE_DIR / "results" / "cancer_gene_impact_analysis",
        help="Ausgabeordner fuer die Analyse.",
    )
    return parser.parse_args()


def maybe_download_oncokb_gene_list(target_path: Path) -> bool:
    """Lädt die öffentliche OncoKB-Genliste herunter."""
    target_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        with urllib.request.urlopen(ONCOKB_CANDIDATE_URL, timeout=60) as response:
            data = json.load(response)
        target_path.write_text(json.dumps(data, indent=2))
        return True
    except Exception:
        return False


def load_gene_roles(path: Path) -> dict[str, tuple[str, str]]:
    """
    Liefert ein Mapping `GENE -> (deutsche Rolle, OncoKB geneType)`.

    Nur Gene mit klarer Krebsrollen-Annotation werden uebernommen.
    """
    data = json.loads(path.read_text())
    roles: dict[str, tuple[str, str]] = {}
    for item in data:
        gene = (item.get("hugoSymbol") or "").strip().upper()
        gene_type = (item.get("geneType") or "").strip()
        if not gene:
            continue
        if gene_type == "ONCOGENE":
            roles[gene] = (ROLE_PROTO, gene_type)
        elif gene_type == "TSG":
            roles[gene] = (ROLE_TSG, gene_type)
        elif gene_type == "ONCOGENE_AND_TSG":
            roles[gene] = (ROLE_DUAL, gene_type)
    return roles


def iter_maf_files(maf_root: Path) -> list[Path]:
    if maf_root.is_file():
        return [maf_root]
    if not maf_root.exists():
        raise FileNotFoundError(f"MAF-Pfad nicht gefunden: {maf_root}")
    maf_files = sorted(maf_root.rglob("*.maf")) + sorted(maf_root.rglob("*.maf.gz"))
    if not maf_files:
        raise FileNotFoundError(f"Keine MAF-Dateien unter {maf_root} gefunden.")
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


def choose_alt_allele(row: dict[str, str]) -> str:
    ref = (row.get("Reference_Allele") or "").strip().upper()
    allele2 = (row.get("Tumor_Seq_Allele2") or "").strip().upper()
    allele1 = (row.get("Tumor_Seq_Allele1") or "").strip().upper()

    if allele2 and allele2 != ref:
        return allele2
    if allele1 and allele1 != ref:
        return allele1
    return allele2


def compute_vaf(t_ref_count: str, t_alt_count: str) -> str:
    """Berechnet die Tumor-VAF, falls Zaehlwerte vorliegen."""
    try:
        ref = float(t_ref_count) if t_ref_count else 0.0
        alt = float(t_alt_count) if t_alt_count else 0.0
    except ValueError:
        return ""
    total = ref + alt
    if total <= 0:
        return ""
    return f"{alt / total:.6f}"


def validate_ref_against_hg38(
    fasta: IndexedFasta,
    chrom: str,
    start_position: int,
    ref_allele_maf: str,
) -> tuple[str, str]:
    """
    Vergleicht das Referenzallel aus der MAF-Datei mit hg38.

    Rueckgabe:
    - hg38_ref_segment
    - ref_matches_hg38 als Text (`yes`, `no`, `not_checked`)
    """
    ref = ref_allele_maf.upper()
    if not ref or any(base not in "ACGT" for base in ref):
        return "", "not_checked"

    hg38_ref = fasta.fetch_sequence(chrom, start_position, start_position + len(ref) - 1)
    if hg38_ref is None:
        return "", "not_checked"
    return hg38_ref, "yes" if hg38_ref == ref else "no"


def classify_role_effect(gene_role: str, variant_classification: str) -> tuple[str, str]:
    """
    Interpretiert den Variantentyp im Kontext der Genrolle.

    `role_consistency` ist eine grobe maschinenlesbare Kategorie.
    `impact_interpretation` ist die sprechendere Beschreibung fuer den Bericht.
    """
    variant_classification = variant_classification.strip()

    if variant_classification in NEUTRAL_OR_NONCODING_CLASSES:
        return (
            "neutral_or_noncoding",
            "wahrscheinlich neutral oder nicht-kodierend",
        )

    if gene_role == ROLE_PROTO:
        if variant_classification in ACTIVATING_LIKE_CLASSES:
            return (
                "role_consistent",
                "mit Proto-Onkogen vereinbar (aktivierende/proteinveraendernde Aenderung)",
            )
        if variant_classification in LOSS_OF_FUNCTION_CLASSES:
            return (
                "role_inconsistent",
                "eher loss-of-function in Proto-Onkogen",
            )
        return (
            "unclear",
            "weitere proteinaendernde Aenderung im Proto-Onkogen",
        )

    if gene_role == ROLE_TSG:
        if variant_classification in LOSS_OF_FUNCTION_CLASSES:
            return (
                "role_consistent",
                "mit Inaktivierung eines Tumorsuppressorgens vereinbar",
            )
        if variant_classification in ACTIVATING_LIKE_CLASSES:
            return (
                "unclear",
                "missense/in-frame Aenderung im Tumorsuppressorgen",
            )
        return (
            "unclear",
            "weitere proteinaendernde Aenderung im Tumorsuppressorgen",
        )

    if variant_classification in LOSS_OF_FUNCTION_CLASSES:
        return (
            "context_dependent",
            "loss-of-function in kontextabhaengigem Krebsgen",
        )
    if variant_classification in ACTIVATING_LIKE_CLASSES:
        return (
            "context_dependent",
            "aktivierende/proteinveraendernde Aenderung in kontextabhaengigem Krebsgen",
        )
    return (
        "context_dependent",
        "weitere Aenderung in kontextabhaengigem Krebsgen",
    )


def build_record(
    row: dict[str, str],
    maf_path: Path,
    fasta: IndexedFasta,
    gene_roles: dict[str, tuple[str, str]],
) -> CancerGeneMutation | None:
    """Erzeugt eine annotierte Mutation in einem Krebsgen oder `None`."""
    gene = (row.get("Hugo_Symbol") or "").strip().upper()
    if not gene or gene not in gene_roles:
        return None

    chrom = (row.get("Chromosome") or "").strip()
    start_raw = (row.get("Start_Position") or "").strip()
    end_raw = (row.get("End_Position") or "").strip()
    try:
        start_position = int(float(start_raw))
        end_position = int(float(end_raw)) if end_raw else start_position
    except ValueError:
        return None

    gene_role, oncokb_gene_type = gene_roles[gene]
    variant_classification = (row.get("Variant_Classification") or "").strip() or "NA"
    role_consistency, interpretation = classify_role_effect(gene_role, variant_classification)
    ref_allele_maf = (row.get("Reference_Allele") or "").strip().upper()
    alt_allele_maf = choose_alt_allele(row)
    hg38_ref_segment, ref_matches_hg38 = validate_ref_against_hg38(
        fasta=fasta,
        chrom=chrom,
        start_position=start_position,
        ref_allele_maf=ref_allele_maf,
    )

    sample_id = (row.get("Tumor_Sample_Barcode") or "").strip() or "UNKNOWN_SAMPLE"
    case_id = (row.get("case_id") or "").strip() or maf_path.parent.name
    patient_id = case_id
    t_ref_count = (row.get("t_ref_count") or "").strip()
    t_alt_count = (row.get("t_alt_count") or "").strip()

    return CancerGeneMutation(
        gene=gene,
        gene_role=gene_role,
        oncokb_gene_type=oncokb_gene_type,
        patient_id=patient_id,
        sample_id=sample_id,
        case_id=case_id,
        chromosome=chrom,
        start_position=start_position,
        end_position=end_position,
        variant_type=(row.get("Variant_Type") or "").strip() or "NA",
        variant_classification=variant_classification,
        impact=(row.get("IMPACT") or "").strip() or "NA",
        consequence=(row.get("Consequence") or "").strip(),
        one_consequence=(row.get("One_Consequence") or "").strip(),
        hgvsp_short=(row.get("HGVSp_Short") or "").strip(),
        protein_position=(row.get("Protein_position") or "").strip(),
        ref_allele_maf=ref_allele_maf,
        alt_allele_maf=alt_allele_maf,
        hg38_ref_segment=hg38_ref_segment,
        ref_matches_hg38=ref_matches_hg38,
        t_ref_count=t_ref_count,
        t_alt_count=t_alt_count,
        tumor_vaf=compute_vaf(t_ref_count, t_alt_count),
        role_consistency=role_consistency,
        impact_interpretation=interpretation,
        maf_file=maf_path.name,
    )


def analyze_cancer_gene_mutations(
    maf_root: Path,
    fasta: IndexedFasta,
    gene_roles: dict[str, tuple[str, str]],
) -> tuple[list[CancerGeneMutation], int]:
    """Liest alle MAF-Dateien und extrahiert Mutationen in annotierten Krebsgenen."""
    maf_files = iter_maf_files(maf_root)
    records: list[CancerGeneMutation] = []
    total_rows = 0

    for maf_path in maf_files:
        for row in iter_maf_rows(maf_path):
            total_rows += 1
            record = build_record(row, maf_path, fasta, gene_roles)
            if record is not None:
                records.append(record)

    return records, total_rows


def top_labels(counter: Counter[str], limit: int = 3) -> str:
    """Formatiert die haeufigsten Elemente eines Counters fuer Report-Tabellen."""
    if not counter:
        return ""
    return "; ".join(f"{label}:{count}" for label, count in counter.most_common(limit) if label)


def summarize_records(records: list[CancerGeneMutation]) -> list[dict[str, str]]:
    """Verdichtet Mutationen zu einer Gen-Zusammenfassung."""
    grouped: dict[str, dict[str, object]] = defaultdict(
        lambda: {
            "gene_role": "",
            "oncokb_gene_type": "",
            "patients": set(),
            "samples": set(),
            "mutation_count": 0,
            "role_consistency_counts": Counter(),
            "impact_counts": Counter(),
            "variant_classification_counts": Counter(),
            "protein_changes": Counter(),
            "ref_yes": 0,
            "ref_no": 0,
            "ref_not_checked": 0,
        }
    )

    for record in records:
        entry = grouped[record.gene]
        entry["gene_role"] = record.gene_role
        entry["oncokb_gene_type"] = record.oncokb_gene_type
        entry["patients"].add(record.patient_id)
        entry["samples"].add(record.sample_id)
        entry["mutation_count"] += 1
        entry["role_consistency_counts"][record.role_consistency] += 1
        entry["impact_counts"][record.impact] += 1
        entry["variant_classification_counts"][record.variant_classification] += 1
        if (
            record.hgvsp_short
            and record.hgvsp_short != "p.?"
            and not record.hgvsp_short.endswith("=")
            and record.role_consistency != "neutral_or_noncoding"
        ):
            entry["protein_changes"][record.hgvsp_short] += 1
        if record.ref_matches_hg38 == "yes":
            entry["ref_yes"] += 1
        elif record.ref_matches_hg38 == "no":
            entry["ref_no"] += 1
        else:
            entry["ref_not_checked"] += 1

    rows: list[dict[str, str]] = []
    for gene, entry in grouped.items():
        role_consistency_counts: Counter[str] = entry["role_consistency_counts"]  # type: ignore[assignment]
        impact_counts: Counter[str] = entry["impact_counts"]  # type: ignore[assignment]
        variant_classification_counts: Counter[str] = entry["variant_classification_counts"]  # type: ignore[assignment]
        protein_changes: Counter[str] = entry["protein_changes"]  # type: ignore[assignment]
        rows.append(
            {
                "gene": gene,
                "gene_role": str(entry["gene_role"]),
                "oncokb_gene_type": str(entry["oncokb_gene_type"]),
                "mutation_count": str(entry["mutation_count"]),
                "patient_count": str(len(entry["patients"])),  # type: ignore[arg-type]
                "sample_count": str(len(entry["samples"])),  # type: ignore[arg-type]
                "role_consistent_count": str(role_consistency_counts.get("role_consistent", 0)),
                "unclear_count": str(role_consistency_counts.get("unclear", 0)),
                "role_inconsistent_count": str(role_consistency_counts.get("role_inconsistent", 0)),
                "context_dependent_count": str(role_consistency_counts.get("context_dependent", 0)),
                "neutral_or_noncoding_count": str(role_consistency_counts.get("neutral_or_noncoding", 0)),
                "high_impact_count": str(impact_counts.get("HIGH", 0)),
                "moderate_impact_count": str(impact_counts.get("MODERATE", 0)),
                "low_impact_count": str(impact_counts.get("LOW", 0)),
                "modifier_impact_count": str(impact_counts.get("MODIFIER", 0)),
                "top_variant_classifications": top_labels(variant_classification_counts, limit=4),
                "top_protein_changes": top_labels(protein_changes, limit=4),
                "hg38_ref_match_yes": str(entry["ref_yes"]),
                "hg38_ref_match_no": str(entry["ref_no"]),
                "hg38_ref_not_checked": str(entry["ref_not_checked"]),
            }
        )

    rows.sort(
        key=lambda row: (
            int(row["role_consistent_count"]),
            int(row["mutation_count"]),
            int(row["patient_count"]),
        ),
        reverse=True,
    )
    return rows


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def records_to_rows(records: list[CancerGeneMutation]) -> list[dict[str, str]]:
    """Wandelt Dataclass-Objekte in CSV-kompatible Zeilen um."""
    rows = [record.__dict__.copy() for record in records]
    rows.sort(
        key=lambda row: (
            row["gene_role"],
            row["gene"],
            IMPACT_ORDER.get(row["impact"], -1),
            row["chromosome"],
            row["start_position"],
        ),
        reverse=False,
    )
    return rows


def filter_by_role(records: list[CancerGeneMutation], role: str) -> list[CancerGeneMutation]:
    return [record for record in records if record.gene_role == role]


def write_report(
    path: Path,
    total_rows: int,
    records: list[CancerGeneMutation],
    summary_rows: list[dict[str, str]],
) -> None:
    """Schreibt einen kompakten Textbericht mit den wichtigsten Ergebnissen."""
    role_counts = Counter(record.gene_role for record in records)
    consistency_counts = Counter(record.role_consistency for record in records)
    proto_summary = [row for row in summary_rows if row["gene_role"] == ROLE_PROTO][:15]
    tsg_summary = [row for row in summary_rows if row["gene_role"] == ROLE_TSG][:15]
    dual_summary = [row for row in summary_rows if row["gene_role"] == ROLE_DUAL][:10]

    with path.open("w") as handle:
        handle.write("Analyse der Auswirkungen von Mutationen auf Krebsrelevante Gene\n")
        handle.write("=============================================================\n\n")
        handle.write(f"Alle eingelesenen MAF-Zeilen: {total_rows}\n")
        handle.write(f"Mutationen in OncoKB-annotierten Krebsgenen: {len(records)}\n")
        handle.write(f"Proto-Onkogen-Hits: {role_counts.get(ROLE_PROTO, 0)}\n")
        handle.write(f"Tumorsuppressor-Hits: {role_counts.get(ROLE_TSG, 0)}\n")
        handle.write(f"Kontextabhaengige Gen-Hits: {role_counts.get(ROLE_DUAL, 0)}\n\n")

        handle.write("Interpretation nach Genrolle:\n")
        handle.write(
            f"- Rolle-konsistente Treffer: {consistency_counts.get('role_consistent', 0)} "
            "(z. B. aktivierende Missense-Aenderungen in Proto-Onkogenen oder truncierende "
            "Aenderungen in Tumorsuppressorgenen)\n"
        )
        handle.write(f"- Unklare proteinveraendernde Treffer: {consistency_counts.get('unclear', 0)}\n")
        handle.write(f"- Rolle-inkonsistente Treffer: {consistency_counts.get('role_inconsistent', 0)}\n")
        handle.write(f"- Kontextabhaengige Treffer: {consistency_counts.get('context_dependent', 0)}\n")
        handle.write(f"- Wahrscheinlich neutrale/nicht-kodierende Treffer: {consistency_counts.get('neutral_or_noncoding', 0)}\n\n")

        handle.write("Top Proto-Onkogene nach rollen-konsistenten Treffern:\n")
        for row in proto_summary:
            handle.write(
                f"- {row['gene']}: konsistent={row['role_consistent_count']}, "
                f"gesamt={row['mutation_count']}, Patienten={row['patient_count']}, "
                f"Top-Proteinwechsel={row['top_protein_changes'] or 'n/a'}\n"
            )

        handle.write("\nTop Tumorsuppressorgene nach rollen-konsistenten Treffern:\n")
        for row in tsg_summary:
            handle.write(
                f"- {row['gene']}: konsistent={row['role_consistent_count']}, "
                f"gesamt={row['mutation_count']}, Patienten={row['patient_count']}, "
                f"Top-Proteinwechsel={row['top_protein_changes'] or 'n/a'}\n"
            )

        if dual_summary:
            handle.write("\nTop kontextabhaengige Gene:\n")
            for row in dual_summary:
                handle.write(
                    f"- {row['gene']}: kontextabhaengig={row['context_dependent_count']}, "
                    f"gesamt={row['mutation_count']}, Patienten={row['patient_count']}\n"
                )


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    if not args.oncokb_gene_list.exists():
        if not args.download_oncokb_if_missing:
            raise FileNotFoundError(
                "OncoKB-Genliste fehlt. Lege die Datei unter "
                f"{args.oncokb_gene_list} ab oder nutze --download-oncokb-if-missing."
            )
        if not maybe_download_oncokb_gene_list(args.oncokb_gene_list):
            raise RuntimeError("Die OncoKB-Genliste konnte nicht heruntergeladen werden.")

    gene_roles = load_gene_roles(args.oncokb_gene_list)
    fasta = IndexedFasta(args.hg38_fasta, args.hg38_fai)
    try:
        records, total_rows = analyze_cancer_gene_mutations(args.maf_root, fasta, gene_roles)
    finally:
        fasta.close()

    summary_rows = summarize_records(records)
    proto_records = filter_by_role(records, ROLE_PROTO)
    tsg_records = filter_by_role(records, ROLE_TSG)
    dual_records = filter_by_role(records, ROLE_DUAL)

    write_csv(args.outdir / "cancer_gene_mutations_detailed.csv", records_to_rows(records))
    write_csv(args.outdir / "cancer_gene_summary.csv", summary_rows)
    write_csv(args.outdir / "proto_oncogene_mutations.csv", records_to_rows(proto_records))
    write_csv(args.outdir / "tumor_suppressor_mutations.csv", records_to_rows(tsg_records))
    write_csv(args.outdir / "context_dependent_gene_mutations.csv", records_to_rows(dual_records))
    write_csv(
        args.outdir / "proto_oncogene_summary.csv",
        [row for row in summary_rows if row["gene_role"] == ROLE_PROTO],
    )
    write_csv(
        args.outdir / "tumor_suppressor_summary.csv",
        [row for row in summary_rows if row["gene_role"] == ROLE_TSG],
    )
    write_csv(
        args.outdir / "context_dependent_gene_summary.csv",
        [row for row in summary_rows if row["gene_role"] == ROLE_DUAL],
    )
    write_report(args.outdir / "report.txt", total_rows, records, summary_rows)

    print(f"Fertig. Ergebnisse liegen unter: {args.outdir}")
    print(f"Mutationen in relevanten Krebsgenen: {len(records)}")
    print(f"Proto-Onkogen-Hits: {len(proto_records)}")
    print(f"Tumorsuppressor-Hits: {len(tsg_records)}")
    print(f"Kontextabhaengige Gen-Hits: {len(dual_records)}")


if __name__ == "__main__":
    main()
