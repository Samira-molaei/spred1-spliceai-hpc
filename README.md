# SPRED1 SpliceAI HPC Workflow

A reproducible Slurm + Singularity workflow for batch **SpliceAI** scoring of **SPRED1** variants defined on the MANE Select transcript.

## Overview

This repository provides a transcript-specific workflow for scoring SPRED1 variants on an HPC cluster using the Broad SpliceAI Lookup local API. It converts transcript-level HGVS variants to GRCh38 genomic coordinates, runs SpliceAI locally, retains the SPRED1 MANE Select transcript, and exports standardized splice-effect scores for downstream variant curation.

### Target transcript

- **Gene:** `SPRED1`
- **RefSeq:** `NM_152594.3`
- **MANE Select Ensembl transcript:** `ENST00000299084.9`
- **Genome build:** GRCh38 / hg38

Because one genomic variant may overlap multiple SPRED1 isoforms, the workflow intentionally retains only `ENST00000299084.9`, the transcript corresponding to `NM_152594.3`.

## Analysis parameters

| Parameter | Setting |
|---|---|
| Genome | GRCh38 / hg38 |
| GENCODE gene set | `basic` |
| Maximum distance | 500 bp |
| Masking | `mask=0` (unmasked) |
| Target transcript | `ENST00000299084.9` |
| Output precision | 2 decimal places |

SpliceAI maxima are determined from the **raw API values** before scores are formatted for display.

Two-decimal output uses conventional half-up rounding:

```text
0.005 -> 0.01
0.014 -> 0.01
0.211 -> 0.21
```

## Workflow

```text
NM_152594.3 HGVS variants
        |
        v
GeneBe HGVS -> GRCh38 conversion
        |
        v
Local Broad SpliceAI Lookup API
(hg38 / GENCODE basic / distance 500 / mask 0)
        |
        v
Retain ENST00000299084.9 only
        |
        v
DS_AG / DS_AL / DS_DG / DS_DL
        |
        v
Determine maximum from raw scores
        |
        v
Format final scores to 2 decimals
        |
        v
CSV output
```

## Repository structure

```text
.
├── README.md
├── LICENSE
├── requirements.txt
├── .gitignore
├── examples/
│   └── spred1_hgvs_variants.example.txt
└── scripts/
    ├── run_spred1_spliceai_batch.py
    └── run_spred1_spliceai_full.sbatch
```

## Requirements

- Slurm workload manager
- Singularity or Apptainer
- Python 3
- `requests`
- `genebe`
- Broad `spliceai-38` container image

Install the Python dependencies:

```bash
python3 -m pip install -r requirements.txt
```

Create the SpliceAI container image if needed:

```bash
singularity pull spliceai-38.sif docker://weisburd/spliceai-38:latest
```

## Input

Provide one transcript HGVS variant per line:

```text
NM_152594.3:c.6C>G
NM_152594.3:c.7G>C
NM_152594.3:c.27C>A
```

An example input file is provided at:

```text
examples/spred1_hgvs_variants.example.txt
```

By default, the Slurm script expects a working input file named:

```text
spred1_hgvs_variants.txt
```

## Running the workflow

From the working directory:

```bash
sbatch scripts/run_spred1_spliceai_full.sbatch
```

If a non-default Python executable is required:

```bash
PYTHON_BIN=/path/to/python3 \
sbatch scripts/run_spred1_spliceai_full.sbatch
```

Paths can be overridden at submission time:

```bash
SPRED1_WORKDIR=/path/to/workdir \
SPLICEAI_SIF=/path/to/spliceai-38.sif \
SPRED1_INPUT=/path/to/variants.txt \
SPRED1_OUTPUT=/path/to/results.csv \
PYTHON_BIN=/path/to/python3 \
sbatch scripts/run_spred1_spliceai_full.sbatch
```

The Slurm launcher starts a local SpliceAI API instance on the assigned compute node and passes its address to the Python scoring script through `SPLICEAI_API`.

## Batch script behavior

`run_spred1_spliceai_batch.py` performs the following steps:

1. Reads transcript HGVS variants.
2. Resolves each variant to GRCh38 genomic coordinates with GeneBe.
3. Queries the local SpliceAI API using fixed analysis parameters.
4. Selects only `ENST00000299084.9` for SPRED1.
5. Records the four delta scores:
   - `DS_AG` — acceptor gain
   - `DS_AL` — acceptor loss
   - `DS_DG` — donor gain
   - `DS_DL` — donor loss
6. Identifies the maximum delta score from the raw values.
7. Formats reported scores to two decimal places.
8. Writes results incrementally and supports `--resume` for incomplete runs.

## Output

The final CSV contains:

```text
HGVS
resolved_hg38
transcript
transcript_priority
gene
DS_AG
DS_AL
DS_DG
DS_DL
SpliceAI_max
max_type
status
```

Example output values:

```text
HGVS: NM_152594.3:c.6C>G
resolved_hg38: 15-38253191-C-G
transcript: ENST00000299084.9
transcript_priority: MS
gene: SPRED1
DS_AG: 0.01
DS_AL: 0.00
DS_DG: 0.21
DS_DL: 0.00
SpliceAI_max: 0.21
max_type: DS_DG
status: OK
```

## Reproducibility and resume support

Results are written incrementally. When `--resume` is used, the script preserves completed `OK` rows and processes variants that are missing or previously returned an error.

For reproducibility, an existing output file should only be resumed when it was generated with the same genome build, transcript, annotation set, distance, masking, and rounding settings.

## Quality control

Check the number of output rows:

```bash
wc -l spred1_spliceai_results_2dp.csv
```

For example, 539 input variants should produce 540 CSV lines: one header plus 539 data rows.

Check status counts:

```bash
python3 - <<'PY'
import csv
from collections import Counter

with open("spred1_spliceai_results_2dp.csv") as f:
    print(Counter(row["status"] for row in csv.DictReader(f)))
PY
```

## Data scope

This repository contains the reusable computational workflow only. Project-specific curation spreadsheets, unpublished annotations, full working variant lists, generated container images, and analysis result files are excluded.

## Citation

If this workflow is used in research, please cite the original SpliceAI publication and the annotation/resources used in the analysis:

> Jaganathan K, et al. Predicting Splicing from Primary Sequence with Deep Learning. *Cell*. 2019;176(3):535-548.e24.

The workflow also uses GeneBe for HGVS-to-genomic coordinate conversion and the Broad SpliceAI Lookup local API implementation.

## License

Released under the MIT License. See [LICENSE](LICENSE).
