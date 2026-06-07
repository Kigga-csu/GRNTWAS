#!/usr/bin/env python

import argparse
import multiprocessing
import os
import sys
from collections import defaultdict
from time import time

import numpy as np
import pandas as pd
from scipy.stats import chi2

try:
    from natsort import natsorted
except ImportError:
    import re

    def natsorted(values):
        def key(value):
            return [int(part) if part.isdigit() else part for part in re.split(r"(\d+)", str(value))]
        return sorted(values, key=key)


def import_utils(tigar_dir):
    if tigar_dir:
        sys.path.append(tigar_dir)
    try:
        import TIGARutils as tg
    except ImportError:
        if __package__ is None:
            sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        import utils as tg
    return tg


def parse_args():
    parser = argparse.ArgumentParser(
        description="GRNTWAS association test with optional LD correction. "
                    "Provide --LD_pattern to use LD; omit it for no-LD association."
    )
    parser.add_argument("--gene_anno", type=str, dest="annot_path", required=True)
    parser.add_argument("--weight", type=str, dest="w_path", required=True)
    parser.add_argument("--Zscore", type=str, dest="z_path", required=True)
    parser.add_argument(
        "--LD_pattern",
        type=str,
        dest="ld_pattern",
        default=None,
        help="Optional LD file pattern, for example /path/to/CHR{chrom}_reference_cov.nochr.txt.gz.",
    )
    parser.add_argument("--window", type=int, default=1000000)
    parser.add_argument("--weight_threshold", type=float, default=0.0)
    parser.add_argument(
        "--test_stat",
        type=str,
        choices=["SPrediXcan"],
        default="SPrediXcan",
        help="Association statistic. Only SPrediXcan is supported.",
    )
    parser.add_argument("--thread", type=int, default=1)
    parser.add_argument("--out_dir", type=str, required=True)
    parser.add_argument("--out_twas_file", type=str, required=True)
    parser.add_argument("--gtf", type=str, required=True)
    parser.add_argument("--TIGAR_dir", type=str, default=None, help="Optional directory containing TIGARutils.py.")
    return parser.parse_args()


def get_pval(z):
    return np.format_float_scientific(chi2.sf(z ** 2, 1), precision=15, exp_digits=0)


def build_gene_regions(target_weights_df, gene_loc_dict, window, tg):
    regions_unmerged = []
    missing_geneids = []
    for gene_id in target_weights_df["GeneID"].unique():
        if gene_id in gene_loc_dict:
            loc = gene_loc_dict[gene_id]
            start = max(0, int(loc["GeneStart"]) - window)
            end = int(loc["GeneEnd"]) + window
            regions_unmerged.append((str(loc["CHROM"]), start, end))
        else:
            missing_geneids.append(gene_id)
    if missing_geneids:
        print(f"Missing {len(missing_geneids)} GeneID annotations: {', '.join(map(str, missing_geneids[:5]))}")
    return tg.merge_regions(regions_unmerged)


def match_zscores_to_weights(zscore_df, target_weights_df, tg):
    zscore_df = zscore_df.drop_duplicates(subset=["snpID"], keep="first").reset_index(drop=True)
    zscore_df["snpIDflip"] = tg.flip_snpIDs(zscore_df["snpID"].values)
    target_weight_snpids = set(target_weights_df["snpID"])
    matched_indices = zscore_df[
        zscore_df["snpID"].isin(target_weight_snpids) |
        zscore_df["snpIDflip"].isin(target_weight_snpids)
    ].index

    if matched_indices.empty:
        return pd.DataFrame()

    flip_multipliers = {}
    zscore_id_to_use = {}
    final_weight_snps_matched = set()
    for idx, row in zscore_df.loc[matched_indices].iterrows():
        if row["snpID"] in target_weight_snpids:
            flip_multipliers[idx] = 1
            zscore_id_to_use[idx] = row["snpID"]
            final_weight_snps_matched.add(row["snpID"])
        elif row["snpIDflip"] in target_weight_snpids:
            flip_multipliers[idx] = -1
            zscore_id_to_use[idx] = row["snpIDflip"]
            final_weight_snps_matched.add(row["snpIDflip"])

    if not zscore_id_to_use:
        return pd.DataFrame()

    zscore_matched = zscore_df.loc[list(zscore_id_to_use.keys())].copy()
    zscore_matched["snpID"] = zscore_matched.index.map(zscore_id_to_use.get)
    zscore_matched["Zscore"] = zscore_matched["Zscore"] * zscore_matched.index.map(flip_multipliers.get)
    zscore_aligned = zscore_matched[["snpID", "Zscore"]].drop_duplicates(subset=["snpID"], keep="first")

    weights = target_weights_df[target_weights_df["snpID"].isin(final_weight_snps_matched)].copy()
    zw = weights.merge(zscore_aligned, on="snpID", how="inner")
    if zw.empty:
        return zw
    snp_order = natsorted(zw["snpID"].unique())
    return zw.drop_duplicates(subset=["snpID"], keep="first").set_index("snpID").loc[snp_order].reset_index()


def compute_no_ld_stats(zw):
    results = {}
    weights = zw["ES"].astype(float).values
    zscores = zw["Zscore"].astype(float).values

    maf = zw["MAF"].astype(float).clip(0, 1).values
    snp_sd = np.sqrt(2 * maf * (1 - maf))
    denom = np.sqrt(np.dot(weights * snp_sd, weights * snp_sd))
    z_spred = np.sum(weights * snp_sd * zscores) / denom if denom > 0 else np.nan
    results["SPred_Z"] = z_spred
    results["SPred_PVAL"] = get_pval(z_spred) if not pd.isna(z_spred) else np.nan

    return results, zw


def expand_query_snps(snp_ids, tg):
    query_snps = set(snp_ids)
    for snp_id in snp_ids:
        try:
            query_snps.add(tg.flip_snpIDs([snp_id])[0])
        except Exception as exc:
            print(f"Could not generate flipped SNP ID for {snp_id}: {exc}")
    return natsorted(query_snps)


def group_snps_by_region(query_snps, merged_regions):
    region_to_snps = defaultdict(list)
    snp_pos_map = {}
    for snp_id in query_snps:
        parts = str(snp_id).split(":")
        if len(parts) < 2:
            continue
        try:
            snp_chrom = parts[0]
            snp_pos = int(parts[1])
        except ValueError:
            continue
        snp_pos_map[snp_id] = snp_pos
        for region_chrom, region_start, region_end in merged_regions:
            if snp_chrom == str(region_chrom) and int(region_start) <= snp_pos <= int(region_end):
                region_to_snps[f"{region_chrom}:{region_start}-{region_end}"].append(snp_id)
                break
    return region_to_snps, snp_pos_map


def load_ld_rows(args, merged_regions, snp_ids, tg):
    query_snps = expand_query_snps(snp_ids, tg)
    region_to_snps, snp_pos_map = group_snps_by_region(query_snps, merged_regions)
    if not region_to_snps:
        return pd.DataFrame(), snp_pos_map

    all_mcov_dfs = []
    found_ld_files = {}
    for region_key, snps_in_region in region_to_snps.items():
        if not snps_in_region:
            continue
        chrom = region_key.split(":", 1)[0]
        ld_file_path = args.ld_pattern.format(chrom=chrom)
        if chrom not in found_ld_files:
            if not os.path.exists(ld_file_path):
                print(f"Missing LD file for chromosome {chrom}: {ld_file_path}")
                continue
            if not os.path.exists(ld_file_path + ".tbi"):
                print(f"Missing LD index for chromosome {chrom}: {ld_file_path}.tbi")
                continue
            found_ld_files[chrom] = ld_file_path

        snps_sorted = sorted(snps_in_region, key=lambda snp: snp_pos_map.get(snp, float("inf")))
        try:
            mcov_df = tg.get_ld_data(found_ld_files[chrom], snps_sorted)
        except tg.NoTargetDataError:
            continue
        except Exception as exc:
            print(f"Failed reading LD for {region_key}: {exc}")
            continue
        if not mcov_df.empty:
            all_mcov_dfs.append(mcov_df)

    if not all_mcov_dfs:
        return pd.DataFrame(), snp_pos_map

    mcov_all = pd.concat(all_mcov_dfs)
    return (
        mcov_all.reset_index()
        .drop_duplicates(subset=["row", "snpID"], keep="first")
        .set_index("row")
        .drop_duplicates(subset=["snpID"], keep="first"),
        snp_pos_map,
    )


def align_ld_to_weights(args, merged_regions, zw, tg):
    snp_search_ids = natsorted(zw["snpID"].unique())
    mcov_all, _ = load_ld_rows(args, merged_regions, snp_search_ids, tg)
    if mcov_all.empty:
        return pd.DataFrame(), pd.DataFrame(), {}, {}

    mcov_all["snpIDflip"] = tg.flip_snpIDs(mcov_all["snpID"].values)
    mcov_all = mcov_all.dropna(subset=["snpID", "snpIDflip"])
    mcov_lookup = {
        "direct": pd.Series(mcov_all.index, index=mcov_all["snpID"]).to_dict(),
        "flip": pd.Series(mcov_all.index, index=mcov_all["snpIDflip"]).to_dict(),
    }

    target_snp_to_mcov_idx = {}
    target_snp_to_mcov_snp = {}
    used_mcov_indices = set()
    for target_snp in snp_search_ids:
        matched_idx = None
        matched_mcov_snp = None
        if target_snp in mcov_lookup["direct"]:
            matched_idx = mcov_lookup["direct"][target_snp]
            matched_mcov_snp = target_snp
        elif target_snp in mcov_lookup["flip"]:
            matched_idx = mcov_lookup["flip"][target_snp]
            matched_mcov_snp = mcov_all.loc[matched_idx, "snpID"]

        if matched_idx is not None and matched_idx not in used_mcov_indices:
            target_snp_to_mcov_idx[target_snp] = matched_idx
            target_snp_to_mcov_snp[target_snp] = matched_mcov_snp
            used_mcov_indices.add(matched_idx)

    snps_with_ld = natsorted(target_snp_to_mcov_idx.keys())
    if not snps_with_ld:
        return pd.DataFrame(), pd.DataFrame(), {}, {}

    zw_ld = zw[zw["snpID"].isin(snps_with_ld)].copy()
    if zw_ld.empty:
        return pd.DataFrame(), pd.DataFrame(), {}, {}
    zw_ld = zw_ld.set_index("snpID").loc[snps_with_ld].reset_index()
    mcov_final = mcov_all.loc[list(target_snp_to_mcov_idx.values())].copy().sort_index()
    return zw_ld, mcov_final, target_snp_to_mcov_snp, target_snp_to_mcov_idx


def build_ld_blocks(mcov_final, target_snp_to_mcov_snp, zw_ld, tg):
    mcov_final = mcov_final.copy()
    mcov_final["CHROM"] = mcov_final["snpID"].map(lambda snp: str(snp).split(":")[0])
    mcov_snp_to_target_snp = {mcov_snp: target_snp for target_snp, mcov_snp in target_snp_to_mcov_snp.items()}
    weight_by_target_snp = zw_ld.set_index("snpID")["ES"]

    blocks = []
    for chrom, mcov_block in mcov_final.groupby("CHROM"):
        mcov_block = mcov_block.sort_index()
        if mcov_block.empty:
            continue
        try:
            snp_sd, v_cov, _ = tg.get_ld_matrix(mcov_block, return_diag=True)
        except Exception as exc:
            print(f"Failed building LD matrix for chromosome {chrom}: {exc}")
            continue

        ld_snps = mcov_block["snpID"].tolist()
        target_snps = []
        for ld_snp in ld_snps:
            target_snp = mcov_snp_to_target_snp.get(ld_snp)
            if target_snp not in weight_by_target_snp.index:
                target_snps = []
                break
            target_snps.append(target_snp)
        if len(target_snps) != len(ld_snps):
            print(f"Skipping chromosome {chrom}; LD SNP order could not be aligned to weights.")
            continue

        blocks.append({
            "chrom": chrom,
            "target_snps": target_snps,
            "snp_sd": snp_sd,
            "v_cov": v_cov,
        })
    return blocks


def compute_ld_stats(args, merged_regions, zw, tg):
    zw_ld, mcov_final, target_snp_to_mcov_snp, _ = align_ld_to_weights(args, merged_regions, zw, tg)
    if zw_ld.empty or mcov_final.empty:
        return {}, pd.DataFrame()

    blocks = build_ld_blocks(mcov_final, target_snp_to_mcov_snp, zw_ld, tg)
    if not blocks:
        return {}, pd.DataFrame()

    results = {}
    weight_by_target_snp = zw_ld.set_index("snpID")["ES"]

    denominator_sq = 0.0
    snp_sd_map = {}
    for block in blocks:
        w_block = weight_by_target_snp.loc[block["target_snps"]].values
        denom_contrib = np.linalg.multi_dot([w_block, block["v_cov"], w_block])
        denominator_sq += max(float(denom_contrib), 0.0)
        for idx, target_snp in enumerate(block["target_snps"]):
            snp_sd_map[target_snp] = block["snp_sd"][idx]

    if set(zw_ld["snpID"]) == set(snp_sd_map):
        snp_sd_final = np.array([snp_sd_map[snp_id] for snp_id in zw_ld["snpID"]])
        numerator = np.sum(zw_ld["ES"].values * snp_sd_final * zw_ld["Zscore"].values)
        denominator = np.sqrt(denominator_sq) if denominator_sq > 1e-10 else 0
        z_spred = numerator / denominator if denominator > 0 else np.nan
    else:
        z_spred = np.nan
    results["SPred_Z"] = z_spred
    results["SPred_PVAL"] = get_pval(z_spred) if not pd.isna(z_spred) else np.nan

    return results, zw_ld


def read_inputs(args, tg):
    gene_annot_all = pd.read_csv(
        args.annot_path,
        sep="\t",
        usecols=["CHROM", "GeneStart", "GeneEnd", "TargetID", "GeneName"],
        dtype={"CHROM": object, "GeneStart": np.int64, "GeneEnd": np.int64, "TargetID": object, "GeneName": object},
    )
    gene_annot_all = tg.optimize_cols(gene_annot_all)
    target_ids = natsorted(gene_annot_all["TargetID"].unique())

    gene_gtf_all = pd.read_csv(
        args.gtf,
        sep="\t",
        header=None,
        names=["CHROM", "GeneStart", "GeneEnd", "_", "TargetID", "GeneName", "type"],
        usecols=[0, 1, 2, 4, 5],
        dtype={"CHROM": object, "GeneStart": np.int64, "GeneEnd": np.int64, "TargetID": object, "GeneName": object},
    )
    gene_loc_dict = gene_gtf_all.set_index("TargetID")[["CHROM", "GeneStart", "GeneEnd"]].to_dict("index")
    gene_info_dict = gene_gtf_all.set_index("TargetID")[["CHROM", "GeneStart", "GeneEnd", "GeneName"]].to_dict("index")

    weight_cols = ["CHROM", "POS", "snpID", "TargetID", "GeneID", "MAF", "p_HWE", "ES"]
    weight_df_all = pd.read_csv(args.w_path, sep="\t", usecols=weight_cols, low_memory=False)
    weight_df_all = tg.optimize_cols(weight_df_all)
    if args.weight_threshold > 0:
        weight_df_all = weight_df_all[np.abs(weight_df_all["ES"]) > args.weight_threshold].copy()

    weights_by_targetid = {
        target: group_df.reset_index(drop=True)
        for target, group_df in weight_df_all.groupby("TargetID")
        if not group_df.empty
    }
    target_ids = natsorted([target for target in target_ids if target in weights_by_targetid])

    zscore_header = tg.get_header(args.z_path, zipped=True)
    zscore_info = tg.get_cols_dtype(
        zscore_header,
        cols=["CHROM", "POS", "REF", "ALT", "Zscore"],
        ind_namekey=True,
    )
    zscore_info["path"] = args.z_path
    return target_ids, gene_loc_dict, gene_info_dict, weights_by_targetid, zscore_info


def output_columns():
    return ["CHROM", "GeneStart", "GeneEnd", "TargetID", "GeneName", "n_snps", "used_regions",
            "SPred_Z", "SPred_PVAL"]


def write_result(out_twas_path, target, gene_info_dict, merged_regions, zw_used, results):
    base_info = gene_info_dict.get(target, {"CHROM": "NA", "GeneStart": 0, "GeneEnd": 0, "GeneName": "NA"})
    output_df = pd.DataFrame({
        "CHROM": [base_info["CHROM"]],
        "GeneStart": [base_info["GeneStart"]],
        "GeneEnd": [base_info["GeneEnd"]],
        "TargetID": [target],
        "GeneName": [base_info["GeneName"]],
        "n_snps": [zw_used.shape[0]],
        "used_regions": [";".join([f"{r[0]}:{r[1]}-{r[2]}" for r in merged_regions])],
    })
    output_df["SPred_Z"] = results.get("SPred_Z", np.nan)
    output_df["SPred_PVAL"] = results.get("SPred_PVAL", np.nan)
    output_df.to_csv(out_twas_path, sep="\t", index=False, header=False, mode="a", na_rep="NA")


def main():
    start_time = time()
    args = parse_args()
    tg = import_utils(args.TIGAR_dir)
    out_twas_path = os.path.join(args.out_dir, args.out_twas_file)
    os.makedirs(args.out_dir, exist_ok=True)

    mode = "with LD" if args.ld_pattern else "without LD"
    print("###############################################################")
    print(f"GRNTWAS Association {mode}")
    print("###############################################################")
    print(f"Output: {out_twas_path}")

    target_ids, gene_loc_dict, gene_info_dict, weights_by_targetid, zscore_info = read_inputs(args, tg)
    pd.DataFrame(columns=output_columns()).to_csv(out_twas_path, sep="\t", index=False, mode="w")

    @tg.error_handler
    def thread_process(num):
        target = target_ids[num]
        print(f"\n--- {num + 1}/{len(target_ids)}: {target} ---")
        target_weights_df = weights_by_targetid[target]
        merged_regions = build_gene_regions(target_weights_df, gene_loc_dict, args.window, tg)
        if not merged_regions:
            print(f"{target}: no SNP search regions.")
            return None

        all_zscore_dfs = []
        for region_chrom, region_start, region_end in merged_regions:
            try:
                z_df_region = tg.read_tabix(
                    start=str(region_start),
                    end=str(region_end),
                    chrm=str(region_chrom),
                    sampleID=[],
                    raise_error=False,
                    **zscore_info,
                )
                if not z_df_region.empty:
                    all_zscore_dfs.append(z_df_region)
            except Exception as exc:
                print(f"{target}: failed reading GWAS region {region_chrom}:{region_start}-{region_end}: {exc}")

        if not all_zscore_dfs:
            print(f"{target}: no GWAS Z-scores in search regions.")
            return None

        zscore_df = pd.concat(all_zscore_dfs, ignore_index=True)
        zw = match_zscores_to_weights(zscore_df, target_weights_df, tg)
        if zw.empty:
            print(f"{target}: no matching GWAS SNPs.")
            return None

        if args.ld_pattern:
            results, zw_used = compute_ld_stats(args, merged_regions, zw, tg)
        else:
            results, zw_used = compute_no_ld_stats(zw)

        if not results or zw_used.empty:
            print(f"{target}: no association result.")
            return None

        write_result(out_twas_path, target, gene_info_dict, merged_regions, zw_used, results)
        print(f"{target}: association {mode} complete with {zw_used.shape[0]} SNPs.")
        return None

    if args.thread <= 1:
        for num in range(len(target_ids)):
            thread_process(num)
    else:
        with multiprocessing.Pool(args.thread) as pool:
            list(pool.imap_unordered(thread_process, range(len(target_ids))))

    print("\nElapsed: " + tg.format_elapsed_time(time() - start_time))


if __name__ == "__main__":
    main()
