#!/usr/bin/env python3
"""
Lädt offene TCGA-SKCM-MAF-Dateien für Tumorproben aus dem GDC herunter.

Die Abfrage ist auf folgende Kriterien eingeschränkt:
- Projekt: TCGA-SKCM
- Dateityp: Masked Somatic Mutation
- Zugriff: open
- Proben: Tumorproben (Primary Tumor, Metastatic, Recurrent Tumor, Additional - New Primary)

Erzeugt:
- manifest.tsv mit Dateimetadaten
- query_summary.json mit der exakten GDC-Abfrage
- skcm_open_masked_somatic_mutation.tar.gz als Rohdownload
- entpackten Dateien im Zielordner
"""

from __future__ import annotations

import argparse
import csv
import json
import tarfile
import urllib.request
from pathlib import Path
from urllib.parse import urlencode

BASE_DIR = Path(__file__).resolve().parent.parent


GDC_FILES_ENDPOINT = "https://api.gdc.cancer.gov/files"
GDC_DATA_ENDPOINT = "https://api.gdc.cancer.gov/data"
TUMOR_SAMPLE_TYPES = [
    "Primary Tumor",
    "Metastatic",
    "Recurrent Tumor",
    "Additional - New Primary",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Lädt offene TCGA-SKCM-Tumor-MAF-Dateien aus dem GDC herunter."
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=BASE_DIR / "data" / "TCGA-SKCM_open_maf",
        help="Zielordner für Manifest, Archiv und entpackte MAF-Dateien.",
    )
    parser.add_argument(
        "--project-id",
        default="TCGA-SKCM",
        help="GDC-Projekt-ID. Standard: TCGA-SKCM",
    )
    return parser.parse_args()


def build_filters(project_id: str) -> dict[str, object]:
    return {
        "op": "and",
        "content": [
            {"op": "in", "content": {"field": "cases.project.project_id", "value": [project_id]}},
            {"op": "in", "content": {"field": "data_type", "value": ["Masked Somatic Mutation"]}},
            {"op": "in", "content": {"field": "access", "value": ["open"]}},
            {"op": "in", "content": {"field": "cases.samples.sample_type", "value": TUMOR_SAMPLE_TYPES}},
        ],
    }


def query_files(project_id: str) -> dict[str, object]:
    filters = build_filters(project_id)
    fields = [
        "file_id",
        "file_name",
        "md5sum",
        "file_size",
        "data_type",
        "access",
        "cases.case_id",
        "cases.submitter_id",
        "cases.samples.sample_type",
        "cases.samples.submitter_id",
    ]
    params = {
        "filters": json.dumps(filters),
        "fields": ",".join(fields),
        "format": "JSON",
        "size": "5000",
    }
    url = GDC_FILES_ENDPOINT + "?" + urlencode(params)
    with urllib.request.urlopen(url, timeout=120) as response:
        return json.load(response)


def write_manifest(path: Path, hits: list[dict[str, object]]) -> None:
    """Schreibt ein tab-separiertes Manifest mit den wichtigsten Dateimetadaten."""
    fieldnames = [
        "id",
        "filename",
        "md5",
        "size",
        "case_ids",
        "case_submitter_ids",
        "sample_types",
        "sample_submitter_ids",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for hit in hits:
            cases = hit.get("cases") or []
            case_ids = sorted({case.get("case_id", "") for case in cases if case.get("case_id")})
            case_submitter_ids = sorted({case.get("submitter_id", "") for case in cases if case.get("submitter_id")})
            sample_types = sorted(
                {
                    sample.get("sample_type", "")
                    for case in cases
                    for sample in (case.get("samples") or [])
                    if sample.get("sample_type")
                }
            )
            sample_submitter_ids = sorted(
                {
                    sample.get("submitter_id", "")
                    for case in cases
                    for sample in (case.get("samples") or [])
                    if sample.get("submitter_id")
                }
            )
            writer.writerow(
                {
                    "id": hit.get("file_id", ""),
                    "filename": hit.get("file_name", ""),
                    "md5": hit.get("md5sum", ""),
                    "size": hit.get("file_size", ""),
                    "case_ids": ";".join(case_ids),
                    "case_submitter_ids": ";".join(case_submitter_ids),
                    "sample_types": ";".join(sample_types),
                    "sample_submitter_ids": ";".join(sample_submitter_ids),
                }
            )


def download_tarball(file_ids: list[str], output_path: Path) -> None:
    """Lädt alle Dateien gesammelt als tar.gz vom GDC-Data-Endpunkt herunter."""
    payload = json.dumps({"ids": file_ids}).encode("utf-8")
    request = urllib.request.Request(
        GDC_DATA_ENDPOINT,
        data=payload,
        headers={"Content-Type": "application/json"},
        method="POST",
    )
    with urllib.request.urlopen(request, timeout=3600) as response, output_path.open("wb") as handle:
        while True:
            chunk = response.read(1024 * 1024)
            if not chunk:
                break
            handle.write(chunk)


def extract_tarball(archive_path: Path, extract_dir: Path) -> None:
    """Entpackt das heruntergeladene Archiv in einen Unterordner."""
    with tarfile.open(archive_path, "r:gz") as archive:
        archive.extractall(extract_dir, filter="data")


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    response = query_files(args.project_id)
    hits = response["data"]["hits"]
    total = response["data"]["pagination"]["total"]
    file_ids = [hit["file_id"] for hit in hits]

    manifest_path = args.outdir / "manifest.tsv"
    summary_path = args.outdir / "query_summary.json"
    archive_path = args.outdir / "skcm_open_masked_somatic_mutation.tar.gz"
    extract_dir = args.outdir / "files"

    write_manifest(manifest_path, hits)
    summary_path.write_text(
        json.dumps(
            {
                "project_id": args.project_id,
                "expected_file_count": total,
                "downloaded_file_count": len(file_ids),
                "tumor_sample_types": TUMOR_SAMPLE_TYPES,
                "filters": build_filters(args.project_id),
            },
            indent=2,
        )
    )

    if len(file_ids) != total:
        raise RuntimeError(f"Abfrage inkonsistent: erwartet {total}, erhalten {len(file_ids)} Dateieinträge.")

    download_tarball(file_ids, archive_path)
    extract_dir.mkdir(parents=True, exist_ok=True)
    extract_tarball(archive_path, extract_dir)

    print(f"Fertig. Dateien heruntergeladen: {len(file_ids)}")
    print(f"Manifest: {manifest_path}")
    print(f"Archiv: {archive_path}")
    print(f"Entpackt unter: {extract_dir}")


if __name__ == "__main__":
    main()
