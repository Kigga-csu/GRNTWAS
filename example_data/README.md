# Example Data

This directory contains synthetic, public-safe example data for the GRNTWAS training and association workflows.

Files:

- `gene_example_55.bed`: 55 example genes.
- `GRN_example.tsv`: synthetic GRN edge list.
- `GRN_example.gexf`: synthetic GRN graph.
- `sample_ids.txt`: example sample IDs.
- `tpm_normalized_gene_peer_vcf.csv`: synthetic expression matrix.
- `TF_eQTL_example.txt`: synthetic eQTL input.
- `vcf/`: bgzipped and tabix-indexed example VCF.
- `gwas/`: synthetic GWAS Z-score file for association testing.
- `ld/`: synthetic LD covariance file using numeric chromosomes without a `chr` prefix.

The GWAS and LD files can be regenerated after training:

```bash
python scripts/generate_example_association_data.py \
  --weight result/example_weight/weight_GRN.csv \
  --gwas-dir example_data/gwas \
  --ld-dir example_data/ld
```

See the repository-level `README.md` for runnable examples.
