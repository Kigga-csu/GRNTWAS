# Code

This directory contains the GRNTWAS training code and the GWAS association scripts.

Use the repository-level `README.md` for installation, example commands, model selection, and input/output formats.

Main entry points:

- `main.py`: train GRNTWAS expression prediction weights.
- `model/Group_spares_lasso.py`: DPR and GroupLasso + ElasticNet/Lasso/Ridge model training.
- `Association_GWAS/Association_GRNTWAS.py`: run association with LD when `--LD_pattern` is supplied, or without LD when it is omitted.
