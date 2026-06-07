# GRNTWAS

Gene Regulatory Network-guided Transcriptome-Wide Association Study.

This repository contains two connected workflows:

1. Train GRN-guided gene expression prediction models from genotype, expression, gene annotation, eQTL, and GRN data.
2. Run TWAS-style association tests using trained GRNTWAS weights, GWAS summary statistics, and optionally LD reference data.

The default training model is DPR. On systems where DPR cannot run, such as the tested macOS environment here, use `--model elasticnet`, `--model lasso`, `--model ridge`, or `--model auto`.

## Repository Structure

```text
.
├── README.md
├── LICENSE
├── .gitignore
├── code/
│   ├── main.py                              # training CLI entry point
│   ├── config.py                            # default training paths and parameters
│   ├── GRN_guided_adaptive_selection.py     # GRN-guided predictor training workflow
│   ├── GRNutils.py                          # VCF, expression, genotype, and utility functions
│   ├── Regular_subgraph_build.py            # GRN subgraph and TF influence selection
│   ├── GRNet_build.py                       # GRN graph construction helpers
│   ├── pipeline.py                          # legacy/helper pipeline code
│   ├── requirements.txt                     # Python dependencies
│   ├── model/
│   │   ├── DPR                              # provided DPR executable; default model when runnable
│   │   ├── Group_spares_lasso.py            # DPR / GroupLasso + ElasticNet/Lasso/Ridge training
│   │   ├── DAG.png
│   │   └── Y_compose.png
│   └── Association_GWAS/
│       ├── Association_GRNTWAS.py           # association with optional LD reference covariance
│       └── utils.py                         # association I/O, tabix, LD, and SNP utilities
├── example_data/
│   ├── gene_example_55.bed                  # example gene annotation
│   ├── GRN_example.tsv
│   ├── GRN_example.gexf
│   ├── sample_ids.txt
│   ├── tpm_normalized_gene_peer_vcf.csv
│   ├── TF_eQTL_example.txt
│   ├── vcf/                                 # bgzipped/indexed example VCF
│   ├── gwas/                                # synthetic GWAS Z-score file for association demo
│   └── ld/                                  # synthetic LD covariance file for association demo
├── scripts/
│   ├── run_training_example.sh              # train example weights on macOS using ElasticNet
│   ├── run_association_example.sh           # run LD and no-LD association examples
│   ├── run_full_example.sh                  # run training then association
│   └── generate_example_association_data.py # generate synthetic GWAS and LD files
```

Generated outputs are written under `result/`, which is ignored by git.

## Installation

Create and activate a Python environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
pip install -r code/requirements.txt
```

The training and association examples also require `bgzip` and `tabix`.

On macOS:

```bash
brew install htslib
```

Check the environment:

```bash
python -c "import numpy, pandas, scipy, sklearn, networkx, group_lasso; print('Python environment OK')"
tabix --version
bgzip --version
```

## Training Models

Use `--model` to choose the expression prediction model:

- `dpr`: default GRNTWAS model. Requires a runnable DPR executable, usually `code/model/DPR` or a path supplied with `--dpr-path`.
- `elasticnet`: GroupLasso feature selection followed by `ElasticNetCV` with sklearn defaults.
- `lasso`: GroupLasso feature selection followed by `LassoCV` with sklearn defaults.
- `ridge`: GroupLasso feature selection followed by `RidgeCV` with sklearn defaults.
- `auto`: compare available models and keep the one with the best training R2. DPR is included only when the executable is available.

The macOS example script uses `--model elasticnet` because DPR is not runnable on the tested Mac environment. For Linux/HPC runs, keep the default DPR model or pass it explicitly:

```bash
python code/main.py --model dpr --dpr-path code/model/DPR
```

## Quick Start With Example Data

Run training only:

```bash
bash scripts/run_training_example.sh
```

This writes:

```text
result/example_weight/weight_GRN.csv
result/example_info/info_GRN.csv
```

Run association using the trained model, synthetic GWAS data, and synthetic LD data:

```bash
bash scripts/run_association_example.sh
```

This writes:

```text
result/example_association_ld/GRNTWAS_results_AD_2019_with_LD.tsv
result/example_association_no_ld/GRNTWAS_results_AD_2019_no_LD.tsv
```

Run the full example pipeline:

```bash
bash scripts/run_full_example.sh
```

## Training Usage

The main training command is:

```bash
python code/main.py \
  --bed example_data/gene_example_55.bed \
  --grn example_data/GRN_example.gexf \
  --grn-tsv example_data/GRN_example.tsv \
  --geno example_data/vcf/ \
  --exp example_data/tpm_normalized_gene_peer_vcf.csv \
  --sample example_data/sample_ids.txt \
  --eqtl example_data/TF_eQTL_example.txt \
  --out-weight result/example_weight/ \
  --out-info result/example_info/ \
  --threads 1 \
  --method LTM \
  --mode lasso \
  --windows 100000 \
  --tf-numbers 10 \
  --no-cv-r2 \
  --model elasticnet
```

Important training inputs:

- `--bed`: gene BED file. Expected columns are `CHROM`, `GeneStart`, `GeneEnd`, `Strand`, `TargetID`, `GeneName`, and `GeneType`.
- `--grn`: GRN graph in GEXF format.
- `--grn-tsv`: GRN edge list in TSV format.
- `--geno`: directory containing bgzipped and tabix-indexed VCF files.
- `--exp`: expression matrix with gene metadata columns followed by sample columns.
- `--sample`: sample IDs used for genotype and expression alignment.
- `--eqtl`: eQTL file for filtering/selecting SNPs.
- `--model`: `dpr`, `elasticnet`, `lasso`, `ridge`, or `auto`.

Training outputs:

- `weight_GRN.csv`: trained SNP weights. Main columns include `CHROM`, `POS`, `snpID`, `TargetID`, `GeneID`, `MAF`, `p_HWE`, and `ES`.
- `info_GRN.csv`: gene-level training summary. It includes `BestModel`, so downstream users can see whether DPR, ElasticNet, Lasso, or Ridge produced the saved model.

## Association Usage

### With LD

Use `Association_GRNTWAS.py` with `--LD_pattern` when LD reference covariance files are available:

```bash
python code/Association_GWAS/Association_GRNTWAS.py \
  --gene_anno result/example_info/info_GRN.csv \
  --weight result/example_weight/weight_GRN.csv \
  --Zscore example_data/gwas/AD_2019_tigar_GWAS.hg38.sorted.indexed.tsv.gz \
  --LD_pattern example_data/ld/CHR{chrom}_reference_cov.nochr.txt.gz \
  --window 100000 \
  --weight_threshold 0 \
  --thread 1 \
  --out_dir result/example_association_ld \
  --out_twas_file GRNTWAS_results_AD_2019_with_LD.tsv \
  --gtf example_data/gene_example_55.bed
```

### Without LD

Use the same `Association_GRNTWAS.py` command without `--LD_pattern` when LD reference data are unavailable:

```bash
python code/Association_GWAS/Association_GRNTWAS.py \
  --gene_anno result/example_info/info_GRN.csv \
  --weight result/example_weight/weight_GRN.csv \
  --Zscore example_data/gwas/AD_2019_tigar_GWAS.hg38.sorted.indexed.tsv.gz \
  --window 100000 \
  --weight_threshold 0 \
  --thread 1 \
  --out_dir result/example_association_no_ld \
  --out_twas_file GRNTWAS_results_AD_2019_no_LD.tsv \
  --gtf example_data/gene_example_55.bed
```

Association inputs:

- GWAS file: bgzipped and tabix-indexed TSV with columns `CHROM`, `POS`, `REF`, `ALT`, and `Zscore`.
- LD file: bgzipped and tabix-indexed TSV with columns `row`, `CHROM`, `POS`, `snpID`, and `COV`.
- LD chromosome naming: chromosome values are numeric without a `chr` prefix, and the path pattern must contain `{chrom}`, for example `CHR{chrom}_reference_cov.nochr.txt.gz`.
- Weight file: output from GRNTWAS training.
- Gene annotation: usually `info_GRN.csv` from training; `--gtf` points to the BED gene annotation.

The example synthetic GWAS and LD files were generated to mimic this layout:

```text
example_data/gwas/AD_2019_tigar_GWAS.hg38.sorted.indexed.tsv.gz
example_data/gwas/AD_2019_tigar_GWAS.hg38.sorted.indexed.tsv.gz.tbi
example_data/ld/CHR1_reference_cov.nochr.txt.gz
example_data/ld/CHR1_reference_cov.nochr.txt.gz.tbi
```

To regenerate them from trained example weights:

```bash
python scripts/generate_example_association_data.py \
  --weight result/example_weight/weight_GRN.csv \
  --gwas-dir example_data/gwas \
  --ld-dir example_data/ld
```

## Real Data Notes

For real GWAS and LD data, keep the same file formats as the examples:

```text
GWAS: CHROM POS REF ALT Zscore
LD:   row CHROM POS snpID COV
```

The LD pattern should look like:

```text
/path/to/no_chr_prefix/CHR{chrom}_reference_cov.nochr.txt.gz
```

The GWAS file must be sorted, bgzipped, and tabix-indexed by `CHROM` and `POS`.

## License

See `LICENSE`.
