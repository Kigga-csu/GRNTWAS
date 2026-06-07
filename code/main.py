#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os
import warnings
import multiprocessing

import numpy as np
import pandas as pd
import networkx as nx
from sklearn.model_selection import KFold

import config
import GRNutils
import GRN_guided_adaptive_selection


def parse_args():
    parser = argparse.ArgumentParser(description='GRNTWAS main')

    parser.add_argument('--bed', type=str, default=config.BED_PATH)
    parser.add_argument('--grn', type=str, default=config.GRN_NETWORK_PATH)
    parser.add_argument('--grn-tsv', type=str, default=config.NET_RAW_FILE_PATH)
    parser.add_argument('--geno', type=str, default=config.GENO_PATH)
    parser.add_argument('--exp', type=str, default=config.GENE_EXPRESSION_PATH)
    parser.add_argument('--sample', type=str, default=config.SAMPLE_PATH)
    parser.add_argument('--eqtl', type=str, default=config.EQTL_PATH)
    parser.add_argument('--dpr-path', type=str, default=config.DPR_PATH)
    parser.add_argument(
        '--model',
        type=str,
        default=config.MODEL,
        choices=['dpr', 'elasticnet', 'lasso', 'ridge', 'auto'],
        help='Training model. Default is dpr; use elasticnet/lasso/ridge for non-DPR systems or auto to compare available models.'
    )

    parser.add_argument('--out-weight', type=str, default=config.OUT_WEIGHT_PATH)
    parser.add_argument('--out-info', type=str, default=config.OUT_INFO_PATH)

    parser.add_argument('--threads', type=int, default=config.THREAD)
    parser.add_argument('--windows', type=int, default=config.WINDOWS)
    parser.add_argument('--method', type=str, default=config.INFLUENCE_METHOD)
    parser.add_argument('--tf-numbers', type=int, default=config.TF_NUMBERS)
    parser.add_argument('--mode', type=str, default='lasso', choices=['lasso','eqtl'])

    parser.add_argument('--maf', type=float, default=config.MAF_THRESHOLD)
    parser.add_argument('--missing-rate', type=float, default=config.MISSING_RATE)
    parser.add_argument('--hwe', type=float, default=config.HWE_THRESHOLD)
    parser.add_argument('--cv-r2', action='store_true', default=config.CV_R2)
    parser.add_argument('--no-cv-r2', action='store_false', dest='cv_r2')
    parser.add_argument('--cvr2-threshold', type=float, default=config.CVR2_THRESHOLD)

    # pre-A/A state kept these two switches for TF selection behavior
    parser.add_argument('--tf-select-mode', choices=['all', 'fold'], default='fold')
    parser.add_argument('--max-tf', type=int, default=8)

    return parser.parse_args()


def load_grn_network(net_raw_file_path):
    grn_df = pd.read_csv(net_raw_file_path, sep='\t')
    target_gene = grn_df.iloc[:, 1].unique().tolist()
    tf_gene = grn_df.iloc[:, 0].unique().tolist()
    return target_gene, tf_gene


def load_gene_expression(gene_expression_path):
    gene_exp = pd.read_csv(gene_expression_path, sep='\t')
    gene_names = gene_exp['GeneName'].tolist()
    return gene_exp, gene_names


def filter_genes_by_expression(gene_names, target_genes, tf_genes):
    expression_tf_gene = np.array(list(set(gene_names) & set(tf_genes)))
    expression_tg_gene = list(set(gene_names) & set(target_genes))
    return expression_tf_gene, expression_tg_gene


def match_samples(sample_path, geno_path, gene_exp):
    gene_exp_header = gene_exp.columns.tolist()
    gene_exp_sample = gene_exp_header[5:]

    geno_sample = GRNutils.sampleid_vcf(geno_path)
    sample_id_list = GRN_guided_adaptive_selection.sample_id_build(sample_path)

    sample_id = list(set(sample_id_list) & set(geno_sample) & set(gene_exp_sample))
    return np.array(sample_id)


def setup_cross_validation(sample_id, n_splits=5, random_state=42):
    kf = KFold(n_splits=n_splits, shuffle=True, random_state=random_state)
    kf_splits = [(sample_id[x], sample_id[y]) for x, y in kf.split(sample_id)]
    cv_train_id, cv_test_id = zip(*kf_splits)
    return cv_train_id, cv_test_id


def setup_output_dirs(out_weight_path, out_info_path):
    for path in [out_weight_path, out_info_path]:
        if not os.path.exists(path):
            os.makedirs(path)

    pd.DataFrame(columns=config.WEIGHT_COLS).to_csv(
        out_weight_path + 'weight_GRN.csv', sep='\t', index=False, header=True, mode='w'
    )
    pd.DataFrame(columns=config.INFO_COLS).to_csv(
        out_info_path + 'info_GRN.csv', sep='\t', index=False, header=True, mode='w'
    )


def _worker_process(num, shared):
    row = shared['gene_exp'].iloc[num]
    target_id = row.iloc[3]
    target_name = row.iloc[4]
    chrom = row.iloc[0]
    target_expr = pd.DataFrame([shared['gene_exp'].iloc[num]], columns=shared['gene_exp'].columns)

    print('#' * 65)
    print(f'Processing gene: {target_name} ({target_id}) - Chromosome {chrom}')
    print('#' * 65)

    GRN_guided_adaptive_selection.predictor4vcf_GRN_lasso(
        target_id,
        shared['args'].geno,
        shared['sample_id'],
        target_expr,
        shared['cv_train_id'],
        shared['cv_test_id'],
        shared['args'].windows,
        shared['args'].out_weight,
        shared['args'].out_info,
        shared['args'].method,
        shared['args'].tf_numbers,
        shared['args'].bed,
        shared['args'].sample,
        shared['eqtl_all'],
        shared['expression_tf_gene'],
        shared['gene_exp'],
        shared['graph'],
        shared['args'].exp,
        maf=shared['args'].maf,
        missing_rate=shared['args'].missing_rate,
        hwe=shared['args'].hwe,
        cv_r2=shared['args'].cv_r2,
        cvr2_threshold=shared['args'].cvr2_threshold,
        dpr_path=shared['args'].dpr_path,
        model=shared['args'].model,
        tf_select_mode=shared['args'].tf_select_mode,
        max_tf=shared['args'].max_tf
    )


def main():
    warnings.filterwarnings('ignore')
    args = parse_args()

    target_genes, tf_genes = load_grn_network(args.grn_tsv)
    gene_exp, gene_names = load_gene_expression(args.exp)
    expression_tf_gene, _ = filter_genes_by_expression(gene_names, target_genes, tf_genes)

    sample_id = match_samples(args.sample, args.geno, gene_exp)
    eqtl_all = pd.read_csv(args.eqtl, sep='\t')

    gene_exp = gene_exp.iloc[:, :5].join(gene_exp[sample_id], how='inner')
    graph = nx.read_gexf(args.grn)

    if args.cv_r2:
        cv_train_id, cv_test_id = setup_cross_validation(sample_id)
    else:
        cv_train_id, cv_test_id = None, None

    setup_output_dirs(args.out_weight, args.out_info)

    shared = {
        'args': args,
        'gene_exp': gene_exp,
        'sample_id': sample_id,
        'cv_train_id': cv_train_id,
        'cv_test_id': cv_test_id,
        'eqtl_all': eqtl_all,
        'expression_tf_gene': expression_tf_gene,
        'graph': graph,
    }

    gene_size = len(gene_exp)
    if args.threads <= 1:
        for num in range(gene_size):
            _worker_process(num, shared)
    else:
        with multiprocessing.Pool(args.threads) as pool:
            pool.starmap(_worker_process, [(num, shared) for num in range(gene_size)])


if __name__ == "__main__":
    main()
