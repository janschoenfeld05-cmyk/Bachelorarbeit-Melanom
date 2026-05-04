#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import math
import re
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable

BASE_DIR = Path(__file__).resolve().parent.parent


SBS_SUBSTITUTIONS = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
SBS96_CHANNELS = [
    f"{left}[{sub}]{right}"
    for sub in SBS_SUBSTITUTIONS
    for left in "ACGT"
    for right in "ACGT"
]
SBS96_SET = set(SBS96_CHANNELS)

COMPLEMENT = str.maketrans("ACGT", "TGCA")

# COSMIC URLs: Falls die erste URL veraltet ist, wird die naechste probiert.
COSMIC_CANDIDATE_URLS = [
    "https://cancer.sanger.ac.uk/signatures/documents/596/COSMIC_v3.4_SBS_GRCh38.txt",
    "https://cancer.sanger.ac.uk/signatures/documents/297/COSMIC_v3.3.1_SBS_GRCh38.txt",
    "https://cancer.sanger.ac.uk/signatures/documents/383/COSMIC_v3.2_SBS_GRCh38.txt",
]


def reverse_complement(seq: str) -> str:
    """Reverse-Complement fuer DNA-Sequenzen (A<->T, C<->G)."""
    return seq.translate(COMPLEMENT)[::-1]


def normalize_sbs96(ref: str, alt: str, tri_context: str) -> tuple[str, str] | tuple[None, None]:
    """
    Bringt eine SNV in das kanonische SBS96-Format.

    Ergebnis:
    - channel: z. B. T[C>T]C
    - substitution: z. B. C>T
    """
    ref = ref.upper()
    alt = alt.upper()
    tri_context = tri_context.upper()

    if len(ref) != 1 or len(alt) != 1:
        return None, None
    if ref not in "ACGT" or alt not in "ACGT" or ref == alt:
        return None, None
    if len(tri_context) != 3 or any(base not in "ACGT" for base in tri_context):
        return None, None
    if tri_context[1] != ref:
        return None, None

    if ref in "AG":
        # Purin-zentrierte Ereignisse werden auf die Pyrimidin-Perspektive gespiegelt.
        tri_context = reverse_complement(tri_context)
        ref = tri_context[1]
        alt = alt.translate(COMPLEMENT)

    if ref not in "CT":
        return None, None

    substitution = f"{ref}>{alt}"
    channel = f"{tri_context[0]}[{substitution}]{tri_context[2]}"
    if channel not in SBS96_SET:
        return None, None
    return channel, substitution


def is_valid_base(base: str) -> bool:
    return len(base) == 1 and base in "ACGT"


class IndexedFasta:
    """Leichter FASTA-Reader auf Basis von .fai fuer schnellen Random-Access."""

    def __init__(self, fasta_path: Path, fai_path: Path | None = None):
        self.fasta_path = fasta_path
        self.fai_path = fai_path or Path(str(fasta_path) + ".fai")
        self.index = self._read_fai(self.fai_path)
        self.handle = self.fasta_path.open("rb")
        self.alias_map = self._build_alias_map()

    def close(self) -> None:
        self.handle.close()

    @staticmethod
    def _read_fai(path: Path) -> dict[str, tuple[int, int, int, int]]:
        if not path.exists():
            raise FileNotFoundError(f"FAI-Index nicht gefunden: {path}")
        index: dict[str, tuple[int, int, int, int]] = {}
        with path.open() as fh:
            for line in fh:
                if not line.strip():
                    continue
                name, length, offset, line_bases, line_width = line.rstrip("\n").split("\t")[:5]
                index[name] = (int(length), int(offset), int(line_bases), int(line_width))
        return index

    def _build_alias_map(self) -> dict[str, str]:
        aliases: dict[str, str] = {}
        # Direkte Namen
        for name in self.index:
            aliases[name] = name
            aliases[name.upper()] = name

        # Mapping chr1..chr22, chrX, chrY, chrM -> RefSeq Accession in GRCh38
        refseq_map = {
            **{str(i): f"NC_{i:06d}.{11 if i == 1 else 12 if i in {2, 3, 4, 6, 7, 12, 13, 14, 15, 16, 17, 18, 20} else 10 if i in {5, 8, 11} else 11 if i in {9, 10} else 9 if i in {19, 22} else 8 if i in {21} else 0}" for i in range(1, 23)},
            "X": "NC_000023.11",
            "Y": "NC_000024.10",
            "M": "NC_012920.1",
            "MT": "NC_012920.1",
        }

        # Korrekturen fuer einige Versionen (explizit)
        refseq_map.update(
            {
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
            }
        )

        for chrom, accession in refseq_map.items():
            if accession in self.index:
                aliases[chrom] = accession
                aliases[f"chr{chrom}"] = accession
                aliases[f"CHR{chrom}"] = accession
        return aliases

    def resolve_chrom(self, chrom: str) -> str | None:
        chrom = chrom.strip()
        if not chrom:
            return None
        if chrom in self.alias_map:
            return self.alias_map[chrom]
        key = chrom.upper()
        return self.alias_map.get(key)

    def fetch_base(self, chrom: str, pos1: int) -> str | None:
        # Byte-Offset-Berechnung gemaess .fai-Spezifikation.
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
        tri = left + center + right
        return tri if re.fullmatch(r"[ACGT]{3}", tri) else None


def iter_maf_files(maf_input: Path) -> list[Path]:
    """Sammelt alle MAF-Dateien rekursiv oder nutzt eine einzelne Eingabedatei."""
    if maf_input.is_file():
        return [maf_input]
    if not maf_input.exists():
        raise FileNotFoundError(f"MAF-Pfad nicht gefunden: {maf_input}")
    files = sorted(maf_input.rglob("*.maf")) + sorted(maf_input.rglob("*.maf.gz"))
    if not files:
        raise FileNotFoundError(f"Keine .maf/.maf.gz Dateien unter {maf_input} gefunden.")
    return files


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open("rt")


def detect_delimiter(header_line: str) -> str:
    return "\t" if header_line.count("\t") >= header_line.count(",") else ","


def normalize_vector(values: dict[str, float]) -> dict[str, float]:
    """Normiert einen Vektor auf Summe 1 (falls Summe > 0)."""
    total = sum(values.values())
    if total <= 0:
        return {k: 0.0 for k in values}
    return {k: v / total for k, v in values.items()}


def cosine_similarity(vec_a: dict[str, float], vec_b: dict[str, float]) -> float:
    """Aehnlichkeit zweier SBS96-Profile (1 = identisch, 0 = orthogonal)."""
    dot = 0.0
    norm_a = 0.0
    norm_b = 0.0
    for channel in SBS96_CHANNELS:
        a = vec_a.get(channel, 0.0)
        b = vec_b.get(channel, 0.0)
        dot += a * b
        norm_a += a * a
        norm_b += b * b
    if norm_a <= 0 or norm_b <= 0:
        return 0.0
    return dot / math.sqrt(norm_a * norm_b)


def maybe_download_cosmic(target_path: Path) -> bool:
    """Versucht mehrere bekannte COSMIC-URLs, bis eine erreichbar ist."""
    target_path.parent.mkdir(parents=True, exist_ok=True)
    for url in COSMIC_CANDIDATE_URLS:
        try:
            with urllib.request.urlopen(url, timeout=30) as response:
                content = response.read()
            target_path.write_bytes(content)
            return True
        except Exception:
            continue
    return False


def load_cosmic_signatures(cosmic_path: Path) -> dict[str, dict[str, float]]:
    """
    Laedt COSMIC-Signaturen robust fuer zwei gaengige Matrix-Layouts:
    - Kanaele in Zeilen, Signaturen in Spalten
    - Signaturen in Zeilen, Kanaele in Spalten
    """
    with cosmic_path.open() as fh:
        first_line = fh.readline()
    delim = detect_delimiter(first_line)

    with cosmic_path.open() as fh:
        reader = list(csv.reader(fh, delimiter=delim))

    if not reader or len(reader) < 2:
        raise ValueError(f"COSMIC-Datei leer oder ungueltig: {cosmic_path}")

    header = [cell.strip() for cell in reader[0]]
    row_first_col = [row[0].strip() for row in reader[1:] if row]
    row_is_channels = sum(1 for x in row_first_col if x in SBS96_SET) >= 80
    col_is_channels = sum(1 for x in header[1:] if x in SBS96_SET) >= 80

    signatures: dict[str, dict[str, float]] = {}

    if row_is_channels:
        sig_names = header[1:]
        for sig in sig_names:
            signatures[sig] = {channel: 0.0 for channel in SBS96_CHANNELS}
        for row in reader[1:]:
            if not row:
                continue
            channel = row[0].strip()
            if channel not in SBS96_SET:
                continue
            for idx, sig in enumerate(sig_names, start=1):
                try:
                    val = float(row[idx]) if idx < len(row) and row[idx] else 0.0
                except ValueError:
                    val = 0.0
                signatures[sig][channel] = val
    elif col_is_channels:
        channel_cols = [c.strip() for c in header[1:]]
        for row in reader[1:]:
            if not row:
                continue
            sig = row[0].strip()
            if not sig:
                continue
            vec = {channel: 0.0 for channel in SBS96_CHANNELS}
            for idx, channel in enumerate(channel_cols, start=1):
                if channel not in SBS96_SET:
                    continue
                try:
                    val = float(row[idx]) if idx < len(row) and row[idx] else 0.0
                except ValueError:
                    val = 0.0
                vec[channel] = val
            signatures[sig] = vec
    else:
        raise ValueError(
            "COSMIC-Format nicht erkannt. Erwartet wird eine SBS96-Matrix "
            "mit Kanaelen als Zeilen oder Spalten."
        )

    return {sig: normalize_vector(vec) for sig, vec in signatures.items()}


def compare_with_cosmic(
    profile_counts: dict[str, Counter[str]],
    signatures: dict[str, dict[str, float]],
    top_n: int = 5,
) -> dict[str, list[tuple[str, float]]]:
    """Vergleicht Profile mit allen COSMIC-Signaturen via Cosine Similarity."""
    result: dict[str, list[tuple[str, float]]] = {}
    for entity, counts in profile_counts.items():
        vec = normalize_vector({ch: float(counts.get(ch, 0)) for ch in SBS96_CHANNELS})
        sims = [(sig, cosine_similarity(vec, sig_vec)) for sig, sig_vec in signatures.items()]
        sims.sort(key=lambda x: x[1], reverse=True)
        result[entity] = sims[:top_n]
    return result


def write_sbs96_table(path: Path, data: dict[str, Counter[str]], include_total: bool = True) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        header = ["id"]
        if include_total:
            header.append("total_sbs")
        header.extend(SBS96_CHANNELS)
        writer.writerow(header)
        for entity in sorted(data):
            counts = data[entity]
            row = [entity]
            if include_total:
                row.append(sum(counts.values()))
            row.extend(counts.get(ch, 0) for ch in SBS96_CHANNELS)
            writer.writerow(row)


def write_cosmic_matches(path: Path, matches: dict[str, list[tuple[str, float]]]) -> None:
    with path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["id", "rank", "signature", "cosine_similarity"])
        for entity in sorted(matches):
            for rank, (sig, sim) in enumerate(matches[entity], start=1):
                writer.writerow([entity, rank, sig, f"{sim:.6f}"])


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Liest MAF-Daten ein, berechnet Mutationsarten, C>T-Anteil, SBS96-Profil, "
            "gleicht Referenzbasen gegen hg38 ab und vergleicht das Profil mit COSMIC-Signaturen."
        )
    )
    parser.add_argument(
        "--maf-input",
        type=Path,
        default=BASE_DIR / "data" / "TCGA-SKCM_open_maf" / "files",
        help="Datei (.maf/.maf.gz) oder Ordner mit MAF-Dateien.",
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
        help="Pfad zur FASTA-Indexdatei (.fai).",
    )
    parser.add_argument(
        "--cosmic-file",
        type=Path,
        default=BASE_DIR / "reference" / "COSMIC_v3.4_SBS_GRCh38.txt",
        help="Lokaler Pfad zur COSMIC-SBS96-Signaturdatei.",
    )
    parser.add_argument(
        "--download-cosmic-if-missing",
        action="store_true",
        help="Falls --cosmic-file fehlt, versuche automatischen Download.",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=BASE_DIR / "results" / "mutation_pipeline",
        help="Ausgabeordner fuer Reports und Tabellen.",
    )
    parser.add_argument(
        "--top-cosmic",
        type=int,
        default=5,
        help="Anzahl der Top-COSMIC-Matches pro Probe.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    maf_files = iter_maf_files(args.maf_input)
    fasta = IndexedFasta(args.hg38_fasta, args.hg38_fai)

    mutation_type_counts = Counter()
    variant_class_counts = Counter()
    substitution_counts = Counter()
    per_sample_sbs96: dict[str, Counter[str]] = defaultdict(Counter)
    cohort_sbs96: Counter[str] = Counter()

    total_rows = 0
    snp_rows = 0
    snp_with_context = 0
    hg38_ref_match = 0
    hg38_ref_mismatch = 0
    skipped_invalid = 0
    skipped_unknown_chrom = 0
    skipped_no_context = 0

    for maf in maf_files:
        with open_text(maf) as fh:
            # Kommentare verwerfen, Header erkennen
            data_lines: Iterable[str] = (line for line in fh if not line.startswith("#"))
            reader = csv.DictReader(data_lines, delimiter="\t")
            for row in reader:
                total_rows += 1
                mutation_type = (row.get("Variant_Type") or "").strip() or "NA"
                variant_class = (row.get("Variant_Classification") or "").strip() or "NA"
                mutation_type_counts[mutation_type] += 1
                variant_class_counts[variant_class] += 1

                if mutation_type != "SNP":
                    continue
                snp_rows += 1

                ref = (row.get("Reference_Allele") or "").strip().upper()
                alt = (row.get("Tumor_Seq_Allele2") or "").strip().upper()
                chrom = (row.get("Chromosome") or "").strip()
                sample = (row.get("Tumor_Sample_Barcode") or "").strip() or "UNKNOWN_SAMPLE"
                pos_raw = (row.get("Start_Position") or "").strip()

                if not (is_valid_base(ref) and is_valid_base(alt)):
                    skipped_invalid += 1
                    continue
                try:
                    pos = int(float(pos_raw))
                except ValueError:
                    skipped_invalid += 1
                    continue

                resolved = fasta.resolve_chrom(chrom)
                if resolved is None:
                    # Chromosom konnte nicht auf hg38-Accession aufgeloest werden.
                    skipped_unknown_chrom += 1
                    continue

                genome_ref = fasta.fetch_base(resolved, pos)
                if not genome_ref or genome_ref == "N":
                    # Keine nutzbare Referenzbase (N oder ausserhalb der Bounds).
                    skipped_no_context += 1
                    continue
                if genome_ref != ref:
                    hg38_ref_mismatch += 1
                    # Nutze die FASTA-Referenz fuer den Kontext trotzdem weiter.
                    ref_for_profile = genome_ref
                else:
                    hg38_ref_match += 1
                    ref_for_profile = ref

                tri = fasta.fetch_trinucleotide(resolved, pos)
                if tri is None:
                    skipped_no_context += 1
                    continue

                channel, substitution = normalize_sbs96(ref_for_profile, alt, tri)
                if channel is None or substitution is None:
                    skipped_invalid += 1
                    continue

                # Gueltige, hg38-kontextualisierte SNV fliesst in Spektren und 96er-Profil ein.
                snp_with_context += 1
                substitution_counts[substitution] += 1
                per_sample_sbs96[sample][channel] += 1
                cohort_sbs96[channel] += 1

    fasta.close()

    c_to_t_count = substitution_counts.get("C>T", 0)
    total_substitutions = sum(substitution_counts.values())
    c_to_t_fraction = (c_to_t_count / total_substitutions) if total_substitutions else 0.0

    # Tabellen schreiben
    with (outdir / "mutation_type_counts.csv").open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["variant_type", "count", "fraction"])
        for key, value in mutation_type_counts.most_common():
            writer.writerow([key, value, f"{value / total_rows if total_rows else 0:.8f}"])

    with (outdir / "variant_classification_counts.csv").open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["variant_classification", "count", "fraction"])
        for key, value in variant_class_counts.most_common():
            writer.writerow([key, value, f"{value / total_rows if total_rows else 0:.8f}"])

    with (outdir / "substitution_counts.csv").open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["substitution", "count", "fraction"])
        for sub in SBS_SUBSTITUTIONS:
            count = substitution_counts.get(sub, 0)
            writer.writerow([sub, count, f"{count / total_substitutions if total_substitutions else 0:.8f}"])

    write_sbs96_table(outdir / "sbs96_per_sample.csv", per_sample_sbs96, include_total=True)
    write_sbs96_table(outdir / "sbs96_cohort.csv", {"cohort": cohort_sbs96}, include_total=True)

    # COSMIC laden / optional herunterladen
    cosmic_status = "not_run"
    if not args.cosmic_file.exists() and args.download_cosmic_if_missing:
        downloaded = maybe_download_cosmic(args.cosmic_file)
        cosmic_status = "downloaded" if downloaded else "download_failed"
    elif args.cosmic_file.exists():
        cosmic_status = "found_local"
    else:
        cosmic_status = "missing"

    if args.cosmic_file.exists():
        # Nur wenn COSMIC vorliegt, werden Signatur-Matches geschrieben.
        signatures = load_cosmic_signatures(args.cosmic_file)
        cohort_matches = compare_with_cosmic({"cohort": cohort_sbs96}, signatures, top_n=args.top_cosmic)
        sample_matches = compare_with_cosmic(per_sample_sbs96, signatures, top_n=args.top_cosmic)
        write_cosmic_matches(outdir / "cosmic_matches_cohort.csv", cohort_matches)
        write_cosmic_matches(outdir / "cosmic_matches_per_sample.csv", sample_matches)
        cosmic_comp_done = True
    else:
        cosmic_comp_done = False

    with (outdir / "summary.txt").open("w") as fh:
        fh.write("Mutation Pipeline Summary\n")
        fh.write("=========================\n\n")
        fh.write(f"MAF files processed: {len(maf_files)}\n")
        fh.write(f"Total rows: {total_rows}\n")
        fh.write(f"SNP rows: {snp_rows}\n")
        fh.write(f"SNP rows with valid hg38 context: {snp_with_context}\n")
        fh.write(f"hg38 reference matches: {hg38_ref_match}\n")
        fh.write(f"hg38 reference mismatches: {hg38_ref_mismatch}\n")
        fh.write(f"Skipped invalid alleles/positions: {skipped_invalid}\n")
        fh.write(f"Skipped unknown chromosomes: {skipped_unknown_chrom}\n")
        fh.write(f"Skipped missing context (N/out-of-range): {skipped_no_context}\n\n")
        fh.write("Substitution spectrum (canonical pyrimidine view):\n")
        for sub in SBS_SUBSTITUTIONS:
            count = substitution_counts.get(sub, 0)
            fh.write(f"{sub}: {count} ({count / total_substitutions if total_substitutions else 0:.2%})\n")
        fh.write("\n")
        fh.write(f"C>T Anteil: {c_to_t_count}/{total_substitutions} ({c_to_t_fraction:.2%})\n\n")
        fh.write(f"COSMIC status: {cosmic_status}\n")
        fh.write(f"COSMIC comparison completed: {cosmic_comp_done}\n")
        if cosmic_comp_done:
            top = compare_with_cosmic({"cohort": cohort_sbs96}, signatures, top_n=min(5, args.top_cosmic))["cohort"]
            fh.write("Top COSMIC matches (cohort):\n")
            for rank, (sig, sim) in enumerate(top, start=1):
                fh.write(f"{rank}. {sig}: cosine={sim:.4f}\n")

    print(f"Fertig. Ergebnisse unter: {outdir}")
    print(f"C>T Anteil: {c_to_t_fraction:.2%} ({c_to_t_count}/{total_substitutions})")
    if cosmic_comp_done:
        print(f"COSMIC Vergleich erfolgreich: {outdir / 'cosmic_matches_cohort.csv'}")
    else:
        print("COSMIC Vergleich nicht durchgefuehrt (Datei fehlt oder Download nicht moeglich).")


if __name__ == "__main__":
    main()
