# IonTorrent-Varlociraptor Classifier

A reproducible pipeline for classifying variants detected from **Ion Torrent amplicon sequencing data** using [Varlociraptor](https://varlociraptor.github.io/).

This workflow estimates the probability of each variant being **somatic**, **germline**, **an artifact**, **an FFPE artifact**, or **absent**, taking into account the tumor cell fraction and sequencing characteristics of Ion Torrent data.

---

## Overview

The pipeline:

1. Normalizes and filters candidate variants from a provided VCF file.  
2. Estimates alignment properties from the corresponding BAM file.  
3. Splits and preprocesses candidates in parallel using `Varlociraptor`.  
4. Calls variant probabilities under an **FFPE scenario** (with or without contamination modeling).  
5. Postprocesses the output with Python to:
   - Merge relevant fields into a clean TSV file.
   - Compute posterior probabilities and predicted classes.
   - Generate a diagnostic scatter plot of **estimated vs. observed AF**.
   - Print the summarized table to `stdout` (for piping in larger workflows).

Outputs are organized in a structured format:

```
output/
└── <sample_id>/
    ├── data/ # Intermediate data files (alignment properties, candidates, BCFs)
    ├── logs/ # Process logs for each stage
    └── results/ # Final TSV table and scatter plot
```

---

## Requirements

* Install Python dependencies using poetry
* Install Varlociraptor

---

## Usage

The repository provides a wrapper script (`main.sh`) that shows a complete example.

```
bash varlociraptor.sh \
    --id test \
    --tumor_cell_quantity 0.2 \
    --bam /path/to/sample.bam \
    --vcf /path/to/sample.vcf \
    --reference /path/to/reference/hg19.fa \
    --n_cores 32 \
    --max_depth 20000
```

### Arguments

| Option | Description | Required | Default |
|:--------|:-------------|:----------:|:---------|
| `--id` | Sample identifier | ✅ | — |
| `--reference` | Reference genome (FASTA) | ✅ | — |
| `--bam` | Input BAM file | ✅ | — |
| `--vcf` | Input VCF file (candidate variants) | ✅ | — |
| `--n_cores` | Number of parallel cores | ❌ | 1 |
| `--max_depth` | Max coverage (downsamples if higher) | ❌ | 20,000 |
| `--tumor_cell_quantity` | Tumor cell fraction (0–1) | ❌ | 1.0 |
| `--use_mapq` | Use MAPQ for likelihood estimation | ❌ | Off |
| `--scenario` | FFPE model: `ffpe` (uses contamination) or `ffpe_no_cont` | ❌ | ffpe |
| `--output_folder` | Base output directory | ❌ | `./output/` |

---

## Output Description

### Final Outputs (in `results/`)

| File | Description |
|:------|:-------------|
| `varlociraptor_probs.tsv` | Tab-separated file with key columns: variant coordinates, class predictions, probabilities, observed/estimated AF, etc. |
| `estimated_vs_observed_AF.png` | Scatter plot showing **observed** vs **estimated** AF (adjusted for contamination). |

The TSV is also **printed to `stdout`**, enabling integration into pipelines like:

```bash
bash varlociraptor.sh ... | grep "SOMATIC" > somatic_calls.tsv
```

---

## Interpretation Notes

- `PROB_SOMATIC`, `PROB_GERMLINE`, etc. are posterior probabilities output by Varlociraptor.
- The column `pred` indicates the most likely class. If both `GERMLINE` and `SOMATIC` have similar high probabilities, it outputs `GERMLINE/SOMATIC` as the class.
- `AF_OBS_VARL` represents the expected allele frequency after adjusting the estimated value (`AF_VARL`) by the tumor cell fraction.
- The scatter plot helps visualize whether the estimated allele frequencies align with observed AF values.

A warning will be printed to `stderr` if the tumor cell content prior appears underestimated, based on high-AF variants.

---

## Repository Structure

```
.
├── main.sh              # Example run script
├── varlociraptor.sh     # Main workflow script
├── preprocess.py        # Postprocessing and plotting step
├── scenarios/           # YAML scenario templates for Varlociraptor
├── README.md            # This file
├── pyproject.toml       # Poetry file with dependencies
├── poetry.lock          # Poetry lock file for installation
└── .gitignore           # gitignore file
```

---

## Citation

If you use this workflow in your research, please cite:

> Köster, J., Dijkstra, L. J., Marschall, T., & Schönhuth, A. (2020). Varlociraptor: enhancing sensitivity and controlling false discovery rate in somatic indel discovery. *Genome biology, 21(1)*, 98. https://doi.org/10.1186/s13059-020-01993-6

---

## Author

Developed by **Javier Maroto / Pathology department @ Universitätsspital Basel**,  
to facilitate robust probabilistic variant classification for Ion Torrent amplicon sequencing data.

---

## License

This project is released under the **GNU General Public License v3.0**.  
See `LICENSE` for details.