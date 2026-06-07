#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Graph-guided adaptive TF selection (Lasso mode) pipeline wrapper.
Only exposes predictor4vcf_GRN_lasso.
"""

import os
import glob
import operator
import numpy as np
import pandas as pd

# Local imports
import GRNutils as tg
import Regular_subgraph_build as rsb
from model.Group_spares_lasso import compare_lasso_enet_cv_revise


def sample_id_build(human_class_file):
    sample_list = []
    with open(human_class_file, 'r') as file:
        for line in file:
            sample_list.append(line.strip())
    return sample_list


def gene_info_4_id(bed_file_path, gene_ids, windows):
    df = pd.read_csv(bed_file_path, sep='\t', header=None,
                     names=['chrom', 'start', 'end', 'strand', 'gene_id', 'gene_name', 'gene_type'])
    matches = df[df['gene_id'].isin(gene_ids)]
    if matches.empty:
        raise ValueError(f"No gene annotation found for gene IDs: {gene_ids}")

    chrom = matches['chrom'].values
    gene_names = matches['gene_name'].values
    first_match = matches.iloc[0]
    start = max(int(first_match['start']) - windows, 0)
    end = int(first_match['end']) + windows
    return gene_names, start, end, chrom


def _resolve_vcf_path(geno_path, chrom):
    if geno_path and os.path.isfile(geno_path):
        return geno_path
    if geno_path and os.path.isdir(geno_path):
        pattern = os.path.join(geno_path, f"*{chrom}*.vcf.gz")
        candidates = sorted(glob.glob(pattern))
        if candidates:
            return candidates[0]
        any_vcfs = sorted(glob.glob(os.path.join(geno_path, "*.vcf.gz")))
        if any_vcfs:
            return any_vcfs[0]
    return f'{geno_path}NIA_JG_1898_samples_GRM_WGS_b37_JointAnalysis01_2017-12-08_{chrom}.recalibrated_variants.Mayo.hg38.lifted.vcf.gz'


def extract_genotype_vcf(genofile_type="vcf", data_format="GT", bed_file=None,
                         geno_path=None, gene_ID=None, gene_exp_path=None,
                         sampleid_path=None, windows=5000):
    target_name, start, end, chrom = gene_info_4_id(bed_file, [gene_ID], windows)
    chrom_str = str(chrom[0])
    try:
        chrom_int = int(chrom_str)
        if chrom_int not in range(1, 23):
            print(f"{target_name} does not have genotype data")
            return False, False
    except ValueError:
        print(f"{target_name} does not have genotype data")
        return False, False

    genofile_path = _resolve_vcf_path(geno_path, chrom_str)
    ret = tg.sampleid_startup(
        chrm=chrom_str,
        genofile_type=genofile_type,
        data_format=data_format,
        geno_path=genofile_path,
        sampleid_path=sampleid_path,
        geneexp_path=gene_exp_path
    )
    if len(ret) == 4:
        sampleID, sample_size, geno_info, exp_info = ret
    else:
        sampleID, sample_size, geno_info = ret
        exp_info = None
    if type(sampleID) is bool:
        print(f"Sample ID extraction failed")

    df_Geno = tg.read_tabix(str(start), str(end), sampleID, **geno_info)
    return df_Geno, sampleID


def read_tf_geno_vcf(bed_path, geno_path, gene_ids, sample_path, windows, gene_exp_path):
    df_list = []
    for gene_id in gene_ids:
        df, sample_ID = extract_genotype_vcf(
            bed_file=bed_path,
            geno_path=geno_path,
            gene_ID=gene_id,
            sampleid_path=sample_path,
            windows=windows,
            gene_exp_path=gene_exp_path
        )
        if isinstance(df, pd.DataFrame):
            df['GeneID'] = gene_id
            df_list.append(df)
    return df_list


def extract_eqtl(eQTL, TF_ids):
    TF_eQTL_list = []
    for i in TF_ids:
        eQTL_TF = eQTL[eQTL['Gene'] == i]
        if len(eQTL_TF) > 0:
            TF_eQTL_list.append(eQTL_TF)
    if TF_eQTL_list:
        TF_eQTL = pd.concat(TF_eQTL_list, axis=0, ignore_index=True)
        return TF_eQTL
    return False


def Geno_filter_eQTL(TF_Geno_list, eQTL_TF):
    if eQTL_TF is False:
        return TF_Geno_list
    eQTL_POS_list = eQTL_TF['SNPPos'].tolist()
    filtered_list = []
    for TF_Geno in TF_Geno_list:
        filtered_TF_Geno = TF_Geno[TF_Geno['POS'].isin(eQTL_POS_list)]
        if len(filtered_TF_Geno) > 0:
            filtered_list.append(filtered_TF_Geno)
    return filtered_list


def _select_tfs_graph_guided(target_name, target_expr_df, sample_ids,
                             gene_exp, graph, bed_path, expression_TF_gene):
    return rsb.select_TFs_via_graph_guided_bayesian(
        gexf_path=graph,
        target_node_name=target_name,
        target_expr_values=target_expr_df[sample_ids].values.flatten(),
        gene_exp_df=gene_exp,
        sampleID=sample_ids,
        bed_path=bed_path,
        expression_TF_gene_names=expression_TF_gene
    )


def _build_geno_exp_and_meta(geno_list, sample_ids, target_expr,
                             missing_rate, maf, hwe):
    for i, df in enumerate(geno_list):
        if isinstance(df, pd.DataFrame):
            df['group'] = i
        else:
            print(f"Warning: Element at index {i} in Geno_list is not a DataFrame, skipping.")
            continue

    try:
        Geno = pd.concat(geno_list, axis=0, ignore_index=True)
    except ValueError as e:
        print(f"Error: Cannot merge data - {e}")
        return None, None, None, None, None

    if Geno.empty:
        print("Merged DataFrame is empty")
        return None, None, None, None, None

    Geno.drop_duplicates(subset='snpID', keep='first', inplace=True)
    Geno = tg.handle_missing_wsx(Geno, sample_ids, missing_rate)
    if Geno[sample_ids].shape[0] == 0:
        print('This gene is not in genotype data')
        return None, None, None, None, None

    Geno = tg.calc_maf(Geno, sample_ids, maf, op=operator.ge)
    if Geno[sample_ids].shape[0] == 0:
        print('This gene is not in genotype data')
        return None, None, None, None, None

    Geno = tg.calc_p_hwe(Geno, sample_ids, hwe, op=operator.ge)
    if Geno[sample_ids].shape[0] == 0:
        print('This gene is not in genotype data')
        return None, None, None, None, None

    groups = Geno['group'].tolist()
    n_snp = Geno['snpID'].size
    Weight = Geno[['CHROM', 'POS', 'snpID', 'REF', 'ALT', 'p_HWE', 'MAF']].copy()
    Geno = Geno.drop('group', axis=1)
    Geno = tg.center(Geno, sample_ids)
    Geno_meta = Geno[['snpID', 'REF', 'ALT']].copy()

    target_expr = tg.center(target_expr, sample_ids)
    Geno_Exp = pd.concat([
        Geno.set_index(['snpID'])[sample_ids],
        target_expr.set_index(['TargetID'])[sample_ids]
    ]).T

    return Geno_Exp, Geno_meta, groups, Weight, n_snp


def predictor4vcf_GRN_lasso(target_ID, genofile_path, sampleID4exp, target_Expr,
                            CV_trainID, CV_testID, windows,
                            out_weight_path, out_info_path, Influence_M, numbers,
                            bed_path, sample_path, eQTL_file, expression_TF_gene,
                            gene_exp, graph, gene_expression_path,
                            maf, missing_rate, hwe, cv_r2, cvr2_threshold,
                            dpr_path=None, model="dpr", tf_select_mode="fold", max_tf=8):
    """
    Prediction function using Bayesian graph-guided TF selection.
    """
    model = (model or "dpr").lower()
    if dpr_path is None:
        dpr_path = "/data/lab/wangshixian/TIGAR-master/Model_Train_Pred/DPR"
    if model in {"dpr", "auto"} and dpr_path and (not os.path.exists(dpr_path)):
        print(f"Warning: DPR not found at {dpr_path}. Skipping DPR model.")
        dpr_path = None
    elif model not in {"dpr", "auto"} and dpr_path and (not os.path.exists(dpr_path)):
        dpr_path = None

    print('Starting GRN_TWAS calculation (Lasso mode)')
    target_name = target_Expr.iloc[0, 4]
    chrom = target_Expr.iloc[0, 0]
    print(f"Target gene: {target_name}")

    try:
        chrom_int = int(chrom)
        if chrom_int not in range(1, 23):
            print(f"{target_ID} does not have genotype data")
            return 0, 0
    except ValueError:
        print("Chromosome must be an integer")
        return 0, 0

    target_Geno, sampleID = extract_genotype_vcf(
        bed_file=bed_path,
        gene_ID=target_ID,
        geno_path=genofile_path,
        sampleid_path=sample_path,
        windows=windows
    )
    sample_size = sampleID.size

    if target_Geno is not None:
        target_Geno['GeneID'] = target_ID

    def build_geno_list(selected_tf_ids):
        if not selected_tf_ids:
            TF_Geno_list_local = []
        else:
            print(f"TFs selected by Bayesian method: {selected_tf_ids}")
            try:
                TF_Geno_list_local = read_tf_geno_vcf(
                    bed_path, genofile_path, selected_tf_ids,
                    sample_path, 100000, gene_expression_path
                )
            except Exception as e:
                print(f"Error occurred: {e}")
                TF_Geno_list_local = pd.DataFrame()

        if len(TF_Geno_list_local) == 0:
            if target_Geno is None:
                print('No available SNP data')
                return None
            Geno_list_local = [target_Geno]
        else:
            eQTL_TF = extract_eqtl(eQTL_file, selected_tf_ids)
            TF_Geno_list_local = Geno_filter_eQTL(TF_Geno_list_local, eQTL_TF)
            TF_Geno_list_local.insert(0, target_Geno)
            Geno_list_local = TF_Geno_list_local

        print("Checking and filtering empty elements in Geno_list...")
        valid_geno_list_local = [df for df in Geno_list_local if df is not None and not df.empty]

        if not valid_geno_list_local:
            print('No valid genotype DataFrame after filtering, skipping this target.')
            return None

        if len(valid_geno_list_local) < len(Geno_list_local):
            removed = len(Geno_list_local) - len(valid_geno_list_local)
            print(f"Removed {removed} empty/invalid elements from Geno_list.")

        return valid_geno_list_local

    # Cross-validation
    if cv_r2:
        print('Running 5-fold CV...')
        k_fold_R2 = []
        K_best_Models = []
        for i in range(5):
            train_ids = CV_trainID[i]
            test_ids = CV_testID[i]
            tf_ids_source = train_ids if tf_select_mode == "fold" else sampleID4exp
            selected_tf_ids = _select_tfs_graph_guided(
                target_name, target_Expr, tf_ids_source,
                gene_exp, graph, bed_path, expression_TF_gene
            )
            if max_tf and len(selected_tf_ids) > max_tf:
                selected_tf_ids = selected_tf_ids[:max_tf]
            geno_list_cv = build_geno_list(selected_tf_ids)
            if geno_list_cv is None:
                k_fold_R2.append(0)
                K_best_Models.append("NA")
                continue

            Geno_Exp_df_cv, Geno_meta_cv, groups_cv, _, _ = _build_geno_exp_and_meta(
                geno_list_cv, sampleID4exp, target_Expr, missing_rate, maf, hwe
            )
            if Geno_Exp_df_cv is None or Geno_meta_cv is None:
                k_fold_R2.append(0)
                K_best_Models.append("NA")
                continue

            cv_rsquared, best_Model = compare_lasso_enet_cv_revise(
                groups_cv,
                Geno_Exp_df_cv.loc[train_ids],
                Geno_Exp_df_cv.loc[test_ids],
                Geno_meta=Geno_meta_cv,
                dpr_path=dpr_path,
                used_model=model
            )
            k_fold_R2.append(cv_rsquared)
            K_best_Models.append(best_Model)

        k_fold_R2 = np.array(k_fold_R2, dtype=float)
        avg_r2_cv = np.mean(k_fold_R2)
        best_Model = max(set(K_best_Models), key=K_best_Models.count)
        print(f"Average R2 for 5-fold CV: {avg_r2_cv:.4f}")
        print('*' * 36)
        print(f"Best model: {best_Model}")
        if avg_r2_cv < cvr2_threshold:
            print(f'Average R2 < {cvr2_threshold}; Skipping training for TargetID: {target_ID}\n')
            return avg_r2_cv, 0
    else:
        avg_r2_cv = 0
        best_Model = None

    # Train final model
    selected_tf_ids_all = _select_tfs_graph_guided(
        target_name, target_Expr, sampleID4exp,
        gene_exp, graph, bed_path, expression_TF_gene
    )
    if max_tf and len(selected_tf_ids_all) > max_tf:
        selected_tf_ids_all = selected_tf_ids_all[:max_tf]
    geno_list_all = build_geno_list(selected_tf_ids_all)
    if geno_list_all is None:
        return None

    Geno_Exp_df_all, Geno_meta_all, groups_all, Weight_all, n_snp = _build_geno_exp_and_meta(
        geno_list_all, sampleID4exp, target_Expr, missing_rate, maf, hwe
    )
    if Geno_Exp_df_all is None or Geno_meta_all is None:
        print("No valid Geno-Expression data for training.")
        return None

    beta, Rsquared, Pvalue, Alpha, Lambda, cvm, final_model = compare_lasso_enet_cv_revise(
        groups_all, Geno_Exp_df_all, Geno_meta=Geno_meta_all, dpr_path=dpr_path, used_model=model
    )

    # Write weight file
    weight_cols = ['CHROM', 'POS', 'snpID', 'TargetID', 'GeneID', 'MAF', 'p_HWE', 'ES']
    out_weight = pd.DataFrame(columns=weight_cols)
    if beta is not None:
        out_weight['CHROM'] = Weight_all['CHROM']
        out_weight['POS'] = Weight_all['POS']
        out_weight['snpID'] = Weight_all['snpID']
        out_weight['TargetID'] = target_ID
        out_weight['GeneID'] = target_ID
        out_weight['MAF'] = Weight_all['MAF']
        out_weight['p_HWE'] = Weight_all['p_HWE']
        out_weight['ES'] = beta
        out_weight.to_csv(out_weight_path + 'weight_GRN.csv', header=None, index=None, sep='\t', mode='a')

    # Write info file (align with GRN_COP_GSL)
    n_effect_snp = int(np.sum(beta != 0)) if beta is not None else 0
    Info = target_Expr[['CHROM', 'GeneStart', 'GeneEnd', 'TargetID', 'GeneName']].copy()
    Info['sample_size'] = sample_size
    Info['n_snp'] = n_snp
    Info['n_effect_snp'] = n_effect_snp
    Info['CVR2'] = avg_r2_cv
    Info['TrainPVALUE'] = Pvalue if Pvalue is not None and not np.isnan(Pvalue) else 'NaN'
    Info['TrainR2'] = Rsquared if n_effect_snp else 0
    Info['k-fold'] = 5
    Info['alpha'] = Alpha
    Info['Lambda'] = Lambda
    Info['cvm'] = cvm
    Info['CVR2_threshold'] = cvr2_threshold if cv_r2 else 0
    Info['BestModel'] = final_model
    Info.to_csv(out_info_path + 'info_GRN.csv', header=None, index=None, sep='\t', mode='a')

    print('Target Elastic Net training completed.\n')
    return avg_r2_cv, Pvalue
