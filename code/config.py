#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
GRNTWAS Configuration Module

All paths and parameters are configured in this file.
Users should modify paths according to their environment before use.
"""

# ==============================================================================
# Input Data Path Configuration (must be modified before use)
# ==============================================================================

# Gene annotation file path (BED format)
BED_PATH = '/path/to/gene.bed'

# Gene regulatory network file path (GEXF format)
GRN_NETWORK_PATH = '/path/to/network.gexf'

# Gene regulatory network raw file path (TSV format)
NET_RAW_FILE_PATH = '/path/to/network.tsv'

# Genotype data directory path (VCF format)
GENO_PATH = '/path/to/genotype/'

# Gene expression data file path
GENE_EXPRESSION_PATH = '/path/to/expression.csv'

# Sample ID file path
SAMPLE_PATH = '/path/to/sample_ids.txt'

# eQTL data file path
EQTL_PATH = '/path/to/eQTL.txt'

# DPR executable path
DPR_PATH = '/path/to/DPR'

# Training model.
# Options: 'dpr' (default/recommended), 'elasticnet', 'lasso', 'ridge', 'auto'
MODEL = 'dpr'

# ==============================================================================
# Output Path Configuration
# ==============================================================================

# Model weight output directory
OUT_WEIGHT_PATH = 'result/weight/'

# Training info output directory
OUT_INFO_PATH = 'result/info/'

# ==============================================================================
# Algorithm Parameters
# ==============================================================================

# Gene window size (bp)
WINDOWS = 100000

# Number of parallel threads
THREAD = 5

# Whether to perform cross-validation R2 calculation
CV_R2 = True

# Cross-validation R2 threshold
CVR2_THRESHOLD = 0.001

# Influence propagation algorithm
# Options: 'RRW' (Restart Random Walk), 'LTM' (Linear Threshold Model), 
#          'KATZ' (Katz Centrality), 'PATH' (Path Influence)
INFLUENCE_METHOD = 'LTM'

# Number of TFs to select
TF_NUMBERS = 10

# Prediction mode
# Options: 'eqtl' (use predictor4vcf_GRN_eQTL), 'lasso' (use predictor4vcf_GRN_lasso)
PREDICTION_MODE = 'eqtl'

# ==============================================================================
# Data Filtering Parameters
# ==============================================================================

# Minor allele frequency threshold
MAF_THRESHOLD = 0.01

# Missing rate threshold
MISSING_RATE = 0.2

# Hardy-Weinberg equilibrium test p-value threshold
HWE_THRESHOLD = 0.001

# ==============================================================================
# Output File Column Definitions
# ==============================================================================

# Weight file columns
WEIGHT_COLS = ['CHROM', 'POS', 'snpID', 'TargetID', 'GeneID', 'MAF', 'p_HWE', 'ES']

# Info file columns
INFO_COLS = [
    'CHROM', 'GeneStart', 'GeneEnd', 'TargetID', 'GeneName',
    'sample_size', 'n_snp', 'n_effect_snp', 'CVR2', 'TrainPVALUE',
    'TrainR2', 'k-fold', 'alpha', 'Lambda', 'cvm', 'CVR2_threshold',
    'BestModel'
]
