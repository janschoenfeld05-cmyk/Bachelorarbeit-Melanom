#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
from collections import Counter, defaultdict
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent


# Kleine eingebaute Referenz fuer haeufig relevante Krebsgene.
# Fuer eine umfassendere Annotation kann per --gene-role-file eine eigene TSV-Datei
# mit den Spalten gene und role uebergeben werden.
DEFAULT_GENE_ROLES = {
    "ABL1": "Proto-Onkogen",
    "AKT1": "Proto-Onkogen",
    "ALK": "Proto-Onkogen",
    "APC": "Tumorsuppressorgen",
    "ARID1A": "Tumorsuppressorgen",
    "ARID2": "Tumorsuppressorgen",
    "ATM": "Tumorsuppressorgen",
    "ATRX": "Tumorsuppressorgen",
    "BAP1": "Tumorsuppressorgen",
    "BRAF": "Proto-Onkogen",
    "BRCA1": "Tumorsuppressorgen",
    "BRCA2": "Tumorsuppressorgen",
    "CCND1": "Proto-Onkogen",
    "CDH1": "Tumorsuppressorgen",
    "CDK4": "Proto-Onkogen",
    "CDKN2A": "Tumorsuppressorgen",
    "CTNNB1": "Proto-Onkogen",
    "EGFR": "Proto-Onkogen",
    "ERBB2": "Proto-Onkogen",
    "FBXW7": "Tumorsuppressorgen",
    "FGFR2": "Proto-Onkogen",
    "FGFR3": "Proto-Onkogen",
    "GNA11": "Proto-Onkogen",
    "GNAQ": "Proto-Onkogen",
    "HRAS": "Proto-Onkogen",
    "IDH1": "Proto-Onkogen",
    "JAK2": "Proto-Onkogen",
    "KIT": "Proto-Onkogen",
    "KRAS": "Proto-Onkogen",
    "MAP2K1": "Proto-Onkogen",
    "MAP2K2": "Proto-Onkogen",
    "MET": "Proto-Onkogen",
    "MYC": "Proto-Onkogen",
    "NF1": "Tumorsuppressorgen",
    "NOTCH1": "Kontextabhaengig",
    "NRAS": "Proto-Onkogen",
    "NTRK1": "Proto-Onkogen",
    "PIK3CA": "Proto-Onkogen",
    "PTCH1": "Tumorsuppressorgen",
    "PTEN": "Tumorsuppressorgen",
    "PTPN11": "Proto-Onkogen",
    "RAC1": "Proto-Onkogen",
    "RB1": "Tumorsuppressorgen",
    "RET": "Proto-Onkogen",
    "SMAD4": "Tumorsuppressorgen",
    "SMARCA4": "Tumorsuppressorgen",
    "SMARCB1": "Tumorsuppressorgen",
    "STK11": "Tumorsuppressorgen",
    "TERT": "Proto-Onkogen",
    "TP53": "Tumorsuppressorgen",
    "TSC1": "Tumorsuppressorgen",
    "TSC2": "Tumorsuppressorgen",
    "VHL": "Tumorsuppressorgen",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Analysiert MAF-Dateien, zaehlt haeufig mutierte Gene und annotiert "
            "sie als Proto-Onkogen, Tumorsuppressorgen oder unbekannt."
        )
    )
    parser.add_argument(
        "--maf-root",
        type=Path,
        default=BASE_DIR / "data" / "TCGA-SKCM_open_maf" / "files",
        help="Wurzelverzeichnis mit Unterordnern, die *.maf.gz-Dateien enthalten.",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=BASE_DIR / "results" / "gene_analysis",
        help="Zielverzeichnis fuer Ausgabedateien.",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=50,
        help="Wie viele Top-Gene in den Textbericht aufgenommen werden sollen.",
    )
    parser.add_argument(
        "--gene-role-file",
        type=Path,
        default=None,
        help=(
            "Optionale TSV-Datei mit mindestens den Spalten gene und role, "
            "um die eingebaute Genrollen-Referenz zu erweitern oder zu ueberschreiben."
        ),
    )
    parser.add_argument(
        "--include-silent",
        action="store_true",
        help="Wenn gesetzt, werden auch synonyme/stille Varianten mitgezaehlt.",
    )
    return parser.parse_args()


def load_gene_roles(path: Path | None) -> dict[str, str]:
    roles = dict(DEFAULT_GENE_ROLES)
    if path is None:
        return roles
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if "gene" not in reader.fieldnames or "role" not in reader.fieldnames:
            raise ValueError(
                f"{path} muss die TSV-Spalten 'gene' und 'role' enthalten."
            )
        for row in reader:
            gene = (row.get("gene") or "").strip().upper()
            role = (row.get("role") or "").strip()
            if gene and role:
                roles[gene] = role
    return roles


def iter_maf_rows(maf_path: Path):
    with gzip.open(maf_path, "rt") as handle:
        reader = csv.DictReader(
            (line for line in handle if not line.startswith("#")),
            delimiter="\t",
        )
        for row in reader:
            yield row


def is_counted_variant(row: dict[str, str], include_silent: bool) -> bool:
    hugo = (row.get("Hugo_Symbol") or "").strip()
    if not hugo or hugo == "Unknown":
        return False
    variant_type = (row.get("Variant_Type") or "").strip()
    if variant_type and variant_type not in {"SNP", "DNP", "TNP", "ONP", "INS", "DEL"}:
        return False
    variant_class = (row.get("Variant_Classification") or "").strip()
    if not include_silent and variant_class in {"Silent", "Intron", "IGR", "3'UTR", "5'UTR", "RNA", "lincRNA"}:
        return False
    return True


def analyze_mafs(
    maf_root: Path, include_silent: bool
) -> tuple[Counter[str], dict[str, set[str]], dict[str, set[str]], int, int]:
    maf_files = sorted(maf_root.glob("*/*.maf.gz"))
    if not maf_files:
        raise FileNotFoundError(f"Keine MAF-Dateien unter {maf_root} gefunden.")

    mutation_counts: Counter[str] = Counter()
    gene_to_patients: dict[str, set[str]] = defaultdict(set)
    gene_to_samples: dict[str, set[str]] = defaultdict(set)
    total_rows = 0
    counted_rows = 0

    for maf_path in maf_files:
        for row in iter_maf_rows(maf_path):
            total_rows += 1
            if not is_counted_variant(row, include_silent):
                continue
            counted_rows += 1
            gene = row["Hugo_Symbol"].strip().upper()
            patient = (row.get("case_id") or "").strip() or "UNKNOWN_PATIENT"
            sample = (row.get("Tumor_Sample_Barcode") or "").strip() or "UNKNOWN_SAMPLE"
            mutation_counts[gene] += 1
            gene_to_patients[gene].add(patient)
            gene_to_samples[gene].add(sample)

    return mutation_counts, gene_to_patients, gene_to_samples, total_rows, counted_rows


def write_outputs(
    outdir: Path,
    mutation_counts: Counter[str],
    gene_to_patients: dict[str, set[str]],
    gene_to_samples: dict[str, set[str]],
    gene_roles: dict[str, str],
    top_n: int,
    total_rows: int,
    counted_rows: int,
) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    rows = []
    total_mutations = sum(mutation_counts.values())
    for gene, mutation_count in mutation_counts.most_common():
        patient_count = len(gene_to_patients[gene])
        sample_count = len(gene_to_samples[gene])
        role = gene_roles.get(gene, "Nicht annotiert")
        rows.append(
            {
                "gene": gene,
                "mutation_count": mutation_count,
                "patient_count": patient_count,
                "sample_count": sample_count,
                "mutation_fraction": f"{(mutation_count / total_mutations if total_mutations else 0):.8f}",
                "gene_role": role,
            }
        )

    table_path = outdir / "gene_mutation_summary.csv"
    with table_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "gene",
                "mutation_count",
                "patient_count",
                "sample_count",
                "mutation_fraction",
                "gene_role",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    proto_path = outdir / "proto_oncogenes.csv"
    tsg_path = outdir / "tumor_suppressor_genes.csv"
    other_path = outdir / "other_or_unannotated_genes.csv"
    category_to_path = {
        "Proto-Onkogen": proto_path,
        "Tumorsuppressorgen": tsg_path,
    }
    categorized = {
        "Proto-Onkogen": [],
        "Tumorsuppressorgen": [],
        "other": [],
    }
    for row in rows:
        if row["gene_role"] in categorized:
            categorized[row["gene_role"]].append(row)
        else:
            categorized["other"].append(row)

    for role, path in category_to_path.items():
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(categorized[role])
    with other_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(categorized["other"])

    report_path = outdir / "report.txt"
    with report_path.open("w") as handle:
        handle.write("Analyse haeufig mutierter Gene\n")
        handle.write("============================\n\n")
        handle.write(f"Alle eingelesenen MAF-Zeilen: {total_rows}\n")
        handle.write(f"Fuer die Genanalyse gezaehlte Varianten: {counted_rows}\n")
        handle.write(f"Unterschiedliche mutierte Gene: {len(rows)}\n")
        handle.write(f"Gesamtzahl gezaehlter Mutationen: {total_mutations}\n\n")

        role_counts = Counter(row["gene_role"] for row in rows)
        handle.write("Anzahl Gene pro Rolle:\n")
        for role, count in role_counts.most_common():
            handle.write(f"{role}\t{count}\n")

        handle.write(f"\nTop {top_n} Gene nach Mutationshaeufigkeit:\n")
        for row in rows[:top_n]:
            handle.write(
                f"{row['gene']}\tMutationen={row['mutation_count']}\t"
                f"Patienten={row['patient_count']}\tSamples={row['sample_count']}\t"
                f"Rolle={row['gene_role']}\n"
            )

    print(f"Wrote {table_path}")
    print(f"Wrote {proto_path}")
    print(f"Wrote {tsg_path}")
    print(f"Wrote {other_path}")
    print(f"Wrote {report_path}")


def main() -> None:
    args = parse_args()
    gene_roles = load_gene_roles(args.gene_role_file)
    (
        mutation_counts,
        gene_to_patients,
        gene_to_samples,
        total_rows,
        counted_rows,
    ) = analyze_mafs(args.maf_root, args.include_silent)
    write_outputs(
        outdir=args.outdir,
        mutation_counts=mutation_counts,
        gene_to_patients=gene_to_patients,
        gene_to_samples=gene_to_samples,
        gene_roles=gene_roles,
        top_n=args.top_n,
        total_rows=total_rows,
        counted_rows=counted_rows,
    )


if __name__ == "__main__":
    main()
