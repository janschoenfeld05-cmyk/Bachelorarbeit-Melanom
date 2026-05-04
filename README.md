# Bachelorarbeit Pipeline

Diese Skripte wurden so umgestellt, dass ihre Standardpfade relativ zum
Projektverzeichnis aufgeloest werden:

```python
BASE_DIR = Path(__file__).resolve().parent.parent
```

Dadurch sind keine benutzerspezifischen `/Users/...`-Pfade mehr noetig.

## Projektstruktur

- `data/TCGA-SKCM_open_maf/files`
  Eingabedaten: rekursiv abgelegte `.maf` / `.maf.gz` Dateien
- `reference/`
  Referenzdateien wie COSMIC und OncoKB
- `results/`
  Standard-Ausgabeordner fuer alle Analysen
- `scripts/`
  Analyse- und Plot-Skripte

## Ausfuehrungsreihenfolge

1. `scripts/1-download_tcga_skcm_open_maf.py`
   Optional. Laedt die offenen TCGA-SKCM-MAF-Dateien nach `data/TCGA-SKCM_open_maf/`.
2. `scripts/2-mutation_uv_cosmic_pipeline.py`
   Allgemeine UV-/COSMIC-Uebersicht direkt aus den MAF-Dateien.
3. `scripts/3-build_krebspatienten_sbs96_dataset.py`
   Zentrale SBS96-Pipeline: SNV-Filter, hg38-Kontext, SBS96-Profile, Signaturanteile.
4. `scripts/4-generate_sbs96_extended_reports.py`
   Zusatzreports auf Basis von Schritt 3: Heatmap, Cluster, PCA, Mittelprofile.
5. `scripts/5-plot_sbs96_heatmap.py`
   Optionale separate Heatmap-Erzeugung.
6. `scripts/6-analyze_mutated_genes.py`
   Allgemeine Genanalyse mit kleiner eingebauter Genrollenliste.
7. `scripts/7-analyze_cancer_gene_impacts.py`
   OncoKB-basierte Krebsgenanalyse.
8. `scripts/8-plot_cancer_gene_impacts.py`
   Visualisierung der Ergebnisse aus Schritt 7.
9. `scripts/9-analyze_target_gene_variant_classes.py`
   Zielgenanalyse fuer `BRAF`, `NRAS`, `TP53`, `CDKN2A`, `NF1`, `PTEN`, `RB1`, `RAC1`, `MAP2K1`, `MAP2K2`.

### Eingaben

- MAF-Dateien:
  `data/TCGA-SKCM_open_maf/files`
- hg38-FASTA:
  `GCF_000001405.40_GRCh38.p14_genomic.fna`
- hg38-FAI:
  `GCF_000001405.40_GRCh38.p14_genomic.fna.fai`
- COSMIC:
  `reference/COSMIC_v3.4_SBS_GRCh38.txt`
- OncoKB:
  `reference/oncokb_cancer_gene_list.json`

### Ausgaben

- UV/COSMIC-Pipeline:
  `results/mutation_pipeline`
- SBS96-Hauptpipeline:
  `results/krebspatienten_sbs96`
- Erweiterte SBS96-Reports:
  `results/`
- Allgemeine Genanalyse:
  `results/gene_analysis`
- Krebsgen-Impact-Analyse:
  `results/cancer_gene_impact_analysis`
- Zielgenanalyse:
  `results/target_gene_variant_analysis`
