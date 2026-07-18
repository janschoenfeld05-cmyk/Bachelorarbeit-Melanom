---

editor_options: 
  markdown: 
    wrap: 72
---

# Bachelorarbeit Pipeline

Diese Skripte wurden so umgestellt, dass ihre Standardpfade relativ zum Projektverzeichnis aufgelöst werden:

``` python
BASE_DIR = Path(__file__).resolve().parent.parent
```

Dadurch sind keine benutzerspezifischen `/Users/...`-Pfade mehr nötig.

## Voraussetzungen

Die Standardpfade sind jetzt reproduzierbar, aber die Skripte benötigen auf einem frischen Gerät trotzdem die passenden Daten, Referenzen und ggf. Python-Pakete.

### Python-Pakete

Die benötigten Pakete stehen in [requirements.txt]:

``` bash
python -m pip install -r requirements.txt
```

Verwendet werden:

- `matplotlib`
- `numpy`
- `pandas`
- `scipy`
- `seaborn`

### Referenz- und Eingabedateien

Vor einem Lauf müssen je nach Skript folgende Dateien bzw. Ordner vorhanden sein:

- `data/TCGA-SKCM_open_maf/files` MAF-Dateien
- `GCF_000001405.40_GRCh38.p14_genomic.fna` hg38-Referenzgenom
- `GCF_000001405.40_GRCh38.p14_genomic.fna.fai` Indexdatei zum Referenzgenom
- `reference/COSMIC_v3.4_SBS_GRCh38.txt` COSMIC-Signaturmatrix
- `reference/oncokb_cancer_gene_list.json` OncoKB-Krebsgene-Liste

Hinweise:

- `scripts/1-download_tcga_skcm_open_maf.py` kann die TCGA-SKCM-MAF-Dateien herunterladen.
- `scripts/3-build_krebspatienten_sbs96_dataset.py` kann COSMIC bei Bedarf mit `--download-cosmic-if-missing` nachladen.
- `scripts/7-analyze_cancer_gene_impacts.py` kann OncoKB bei Bedarf mit `--download-oncokb-if-missing` nachladen.
- Das hg38-Referenzgenom und die `.fai`-Datei müssen lokal vorhanden sein.

## Projektstruktur

- `data/TCGA-SKCM_open_maf/files` Eingabedaten: rekursiv abgelegte `.maf` / `.maf.gz` Dateien
- `reference/` Referenzdateien wie COSMIC und OncoKB
- `results/` Standard-Ausgabeordner für alle Analysen
- `scripts/` Analyse- und Plot-Skripte

## Ausführungsreihenfolge

Die Skripte sind numeriert, um die Ausführungsreihenfolge zu visualisieren.

1.  `scripts/1-download_tcga_skcm_open_maf.py` Optional. Lädt die offenen TCGA-SKCM-MAF-Dateien nach `data/TCGA-SKCM_open_maf/`.
2.  `scripts/2-mutation_uv_cosmic_pipeline.py` Allgemeine UV-/COSMIC-Übersicht direkt aus den MAF-Dateien.
3.  `scripts/3-build_krebspatienten_sbs96_dataset.py` Zentrale SBS96-Pipeline: SNV-Filter, hg38-Kontext, SBS96-Profile, Signaturanteile.
4.  `scripts/4-generate_sbs96_extended_reports.py` Zusatzreports auf Basis von Schritt 3: Heatmap, Cluster, PCA, Mittelprofile.
5.  `scripts/5-plot_sbs96_heatmap.py` Optionale separate Heatmap-Erzeugung.
6.  `scripts/6-analyze_mutated_genes.py` Allgemeine Genanalyse mit kleiner eingebauter Genrollenliste.
7.  `scripts/7-analyze_cancer_gene_impacts.py` OncoKB-basierte Krebsgenanalyse.
8.  `scripts/8-plot_cancer_gene_impacts.py` Visualisierung der Ergebnisse aus Schritt 7.
9.  `scripts/9-analyze_target_gene_variant_classes.py` Zielgenanalyse für `BRAF`, `NRAS`, `TP53`, `CDKN2A`, `NF1`, `PTEN`, `RB1`, `RAC1`, `MAP2K1`, `MAP2K2`.

## Abhängigkeiten zwischen den Skripten

Nicht alle Skripte sind direkt auf jedem Gerät startbar. Einige Plot-Skripte erwarten Ergebnisse aus vorherigen Schritten:

- `scripts/4-generate_sbs96_extended_reports.py` benötigt die Ausgaben aus `scripts/3-build_krebspatienten_sbs96_dataset.py`
- `scripts/5-plot_sbs96_heatmap.py` benötigt `results/krebspatienten_sbs96/sbs96_per_patient.csv`
- `scripts/8-plot_cancer_gene_impacts.py` benötigt die Ausgaben aus `scripts/7-analyze_cancer_gene_impacts.py`

Direkt auf Rohdaten laufen:

- `scripts/2-mutation_uv_cosmic_pipeline.py`
- `scripts/3-build_krebspatienten_sbs96_dataset.py`
- `scripts/6-analyze_mutated_genes.py`
- `scripts/7-analyze_cancer_gene_impacts.py`
- `scripts/9-analyze_target_gene_variant_classes.py`

### Eingaben

- MAF-Dateien: `data/TCGA-SKCM_open_maf/files`
- hg38-FASTA: `GCF_000001405.40_GRCh38.p14_genomic.fna`
- hg38-FAI: `GCF_000001405.40_GRCh38.p14_genomic.fna.fai`
- COSMIC: `reference/COSMIC_v3.4_SBS_GRCh38.txt`
- OncoKB: `reference/oncokb_cancer_gene_list.json`

### Ausgaben

- UV/COSMIC-Pipeline: `results/mutation_pipeline`
- SBS96-Hauptpipeline: `results/krebspatienten_sbs96`
- Erweiterte SBS96-Reports: `results/`
- Allgemeine Genanalyse: `results/gene_analysis`
- Krebsgen-Impact-Analyse: `results/cancer_gene_impact_analysis`
- Zielgenanalyse: `results/target_gene_variant_analysis`
