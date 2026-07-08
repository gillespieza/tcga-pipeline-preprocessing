# tcga-pipeline-preprocessing

A Python CLI tool for downloading and preprocessing multi-omics data from [cBioPortal](https://www.cbioportal.org/) for any TCGA study.

## Features

- **Interactive mode** — fuzzy-search studies, checkbox data-type selection
- **CLI args mode** — fully scriptable for batch use
- **Supports**: Clinical, RNA-seq, CNA, Somatic Mutations, RPPA
- **Mutations**: written as both long (MAF) and wide (binary gene matrix)
- **Clinical Cleaning**:
  - Automatically deduplicates clinical records to ensure 1 sample per patient (keeps first sample for statistical independence).
  - Parses survival endpoints (`OS`, `DSS`, `DFS`, `PFS`) into binary event variables (`0.0` for censored, `1.0` for event).
  - Filters out records with missing/invalid Overall Survival data.
- **RNA-seq Preprocessing**:
  - Deduplicates RNA-seq samples to ensure 1 sample per patient (keeps the first sample).
  - Automatically removes genes (columns) containing any missing (`NaN`) values.
- **Caching**: re-running skips the download if data is already present
- **Outputs**: Organized into distinct `raw` and `processed` directories.

## Setup

```bash
cd tcga-pipeline-preprocessing
python -m venv .venv
.venv\Scripts\activate        # Windows
# source .venv/bin/activate   # macOS/Linux
pip install -r requirements.txt
```

## Usage

### Interactive (recommended for first use)

```bash
python tcga_fetch.py
```

You will be prompted to:
1. Search and select a cBioPortal study
2. Choose which data types to fetch
3. Set an output directory
4. Optionally apply log2(x+1) to RNA-seq values

### Command-line arguments (for scripting)

```bash
python tcga_fetch.py \
  --study kirc_tcga_pan_can_atlas_2018 \
  --data clinical rnaseq cna mutations \
  --output ./output

# With log2 transformation on RNA-seq:
python tcga_fetch.py --study kirc_tcga_pan_can_atlas_2018 --data clinical rnaseq --log2

# List all available studies:
python tcga_fetch.py --list-studies

# Force re-download (ignore cache):
python tcga_fetch.py --study kirc_tcga_pan_can_atlas_2018 --data clinical --force

# Write individual files only, skip merged.csv:
python tcga_fetch.py --study kirc_tcga_pan_can_atlas_2018 --data clinical rnaseq --no-merge
```

## Output Structure

Outputs are written into two subdirectories under the study folder:

```
output/
└── kirc_tcga_pan_can_atlas_2018/
    ├── raw/
    │   ├── clinical.csv          # Unaltered, raw clinical attributes (pivoted & merged)
    │   ├── rnaseq.csv            # Unaltered RNA-seq expression matrix
    │   ├── cna.csv               # Unaltered CNA matrix
    │   ├── mutations_long.csv    # Unaltered variant list (MAF format)
    │   └── mutations_wide.csv    # Unaltered binary mutation matrix
    └── processed/
        ├── clinical_cleaned.csv  # Deduplicated, binarized survival endpoints, filtered OS
        ├── rna_clean.csv         # Deduplicated RNA-seq expression matrix with missing values removed
        ├── merged.csv            # All processed/cleaned layers joined on SAMPLE_ID
        └── manifest.json         # Run metadata, sample/column counts
```

## Data Type Notes

| Type | Raw File (`raw/`) | Processed File (`processed/`) | Description & Format |
|---|---|---|---|
| Clinical | `clinical.csv` | `clinical_cleaned.csv` | Raw vs. cleaned (deduplicated, binary survival events, OS-filtered) |
| RNA-seq | `rnaseq.csv` | `rna_clean.csv` | Raw vs. preprocessed (deduplicated, genes with missing values removed) |
| CNA | `cna.csv` | `cna.csv` | Samples × genes (-2 deep del → +2 amp) |
| Mutations (long) | `mutations_long.csv` | `mutations_long.csv` | One row per variant call (MAF) |
| Mutations (wide) | `mutations_wide.csv` | `mutations_wide.csv` | Samples × genes (0/1 binary) |
| RPPA | `rppa.csv` | `rppa.csv` | Samples × proteins (z-scores) |
| Merged | - | `merged.csv` | All processed/cleaned layers joined on `SAMPLE_ID` |

## Data Source

Data is downloaded from [cBioPortal](https://www.cbioportal.org/datasets) via the public REST API and study data packages. No authentication is required for public TCGA studies.

