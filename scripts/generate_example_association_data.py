#!/usr/bin/env python

import argparse
import os
import subprocess

import numpy as np
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description="Generate synthetic GWAS and LD files from GRNTWAS example weights.")
    parser.add_argument("--weight", default="result/example_weight/weight_GRN.csv")
    parser.add_argument("--gwas-dir", default="example_data/gwas")
    parser.add_argument("--ld-dir", default="example_data/ld")
    parser.add_argument("--trait", default="AD_2019")
    return parser.parse_args()


def run_command(command):
    subprocess.check_call(command)


def split_snp_id(snp_id):
    chrom, pos, ref, alt = str(snp_id).split(":")
    return chrom, int(pos), ref, alt


def main():
    args = parse_args()
    os.makedirs(args.gwas_dir, exist_ok=True)
    os.makedirs(args.ld_dir, exist_ok=True)

    weights = pd.read_csv(args.weight, sep="\t")
    weights = weights.drop_duplicates(subset=["snpID"]).copy()
    parts = weights["snpID"].map(split_snp_id)
    weights["CHROM"] = [item[0] for item in parts]
    weights["POS"] = [item[1] for item in parts]
    weights["REF"] = [item[2] for item in parts]
    weights["ALT"] = [item[3] for item in parts]
    weights = weights.sort_values(["CHROM", "POS", "REF", "ALT"]).reset_index(drop=True)

    gwas = weights[["CHROM", "POS", "REF", "ALT"]].copy()
    signal = weights["ES"].astype(float).fillna(0).to_numpy()
    gwas["Zscore"] = np.round(1.5 + signal * 2.0 + np.linspace(-0.3, 0.3, len(gwas)), 6)

    gwas_tsv = os.path.join(args.gwas_dir, f"{args.trait}_tigar_GWAS.hg38.sorted.indexed.tsv")
    gwas_gz = gwas_tsv + ".gz"
    gwas.to_csv(gwas_tsv, sep="\t", index=False)
    run_command(["bgzip", "-f", gwas_tsv])
    run_command(["tabix", "-f", "-S", "1", "-s", "1", "-b", "2", "-e", "2", gwas_gz])

    for chrom, chrom_df in weights.groupby("CHROM", sort=False):
        chrom_df = chrom_df.sort_values("POS").reset_index(drop=True)
        rows = []
        maf = chrom_df["MAF"].astype(float).clip(0.001, 0.999).to_numpy()
        variances = 2 * maf * (1 - maf)
        for idx, row in chrom_df.iterrows():
            cov_values = []
            for j in range(idx, len(chrom_df)):
                if j == idx:
                    cov_values.append(variances[idx])
                else:
                    distance = abs(int(chrom_df.loc[j, "POS"]) - int(row["POS"]))
                    cov_values.append(0.05 * np.exp(-distance / 50000.0))
            rows.append({
                "row": idx,
                "CHROM": row["CHROM"],
                "POS": int(row["POS"]),
                "snpID": row["snpID"],
                "COV": ",".join(f"{value:.8f}" for value in cov_values)
            })

        ld = pd.DataFrame(rows)
        ld_tsv = os.path.join(args.ld_dir, f"CHR{chrom}_reference_cov.nochr.txt")
        ld_gz = ld_tsv + ".gz"
        ld.to_csv(ld_tsv, sep="\t", index=False)
        run_command(["bgzip", "-f", ld_tsv])
        run_command(["tabix", "-f", "-S", "1", "-s", "2", "-b", "3", "-e", "3", ld_gz])

    print(f"Wrote {gwas_gz}")
    print(f"Wrote LD files under {args.ld_dir}")


if __name__ == "__main__":
    main()
