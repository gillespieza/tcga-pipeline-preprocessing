# tcga-pipeline-preprocessing

A Python CLI tool for downloading and preprocessing multi-omics data from [cBioPortal](https://www.cbioportal.org/) for any TCGA study.

## Features

- **Interactive mode** — fuzzy-search studies, checkbox data-type selection
- **CLI args mode** — fully scriptable for batch use
- **Supports**: Clinical, RNA-seq, CNA, Somatic Mutations, RPPA
- **Mutations**: written as both long (MAF) and wide (binary gene matrix)
- **Caching**: re-running skips the download if data is already present
- **Outputs**: individual layer CSVs + a merged wide-format CSV + `manifest.json`

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

```
output/
└── kirc_tcga_pan_can_atlas_2018/
    ├── clinical.csv          # Patient + sample metadata merged
    ├── rnaseq.csv            # Samples × genes expression matrix
    ├── cna.csv               # Samples × genes CNA matrix (-2 to +2)
    ├── mutations_long.csv    # One row per variant (MAF format)
    ├── mutations_wide.csv    # Samples × genes binary mutation matrix
    ├── rppa.csv              # Samples × proteins RPPA z-scores
    ├── merged.csv            # All selected layers joined on SAMPLE_ID
    └── manifest.json         # Study ID, fetch date, row/column counts
```

## Data Type Notes

| Type | File | Format |
|---|---|---|
| Clinical | `clinical.csv` | One row per sample |
| RNA-seq | `rnaseq.csv` | Samples × genes (RSEM or log2+1) |
| CNA | `cna.csv` | Samples × genes (-2 deep del → +2 amp) |
| Mutations (long) | `mutations_long.csv` | One row per variant call (MAF) |
| Mutations (wide) | `mutations_wide.csv` | Samples × genes (0/1 binary) |
| RPPA | `rppa.csv` | Samples × proteins (z-scores) |
| Merged | `merged.csv` | All layers joined on SAMPLE_ID |

## Data Source

Data is downloaded from [cBioPortal](https://www.cbioportal.org/datasets) via the public REST API and study data packages. No authentication is required for public TCGA studies.
