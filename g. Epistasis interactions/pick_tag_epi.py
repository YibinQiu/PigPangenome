#!/usr/bin/env python3
# Script: pick_tag_epi.py
# Date: 2024/12/31
# Purpose: Select representative epistatic effect pairs for each QTL (gwas < 0.001 locus to all).
# Author: Yibin Qiu

import os
import pandas as pd
import numpy as np
import argparse
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor

# Constants
VERSION = "1.0"
# have fun

VERSION = "1.1"
# update 
# 20241113: split 18*18 chromsome matrix and then process

VERSION = "1.2"
# update 
# 20241231: we define updownstream 1MB as one QTL

# Function to print script information
def print_script_info():
    print("*********************************************************************", flush=True)
    print("* pick cherry -- a script: Pick representative epistatic effect pairs", flush=True)
    print(f"* Version {VERSION}", flush=True)
    print("* Yibin Qiu (13422157044qyb@gmail.com)", flush=True)
    print("* South China Agricultural University", flush=True)
    print("*********************************************************************", flush=True)

print_script_info()
#now we define updownstream 1MB as one QTL
def calculate_bounds(lead_pos, chr_name, chr_len_dict, distance=1000000):
    lower_bound = max(lead_pos - distance, 0)
    upper_bound = min(lead_pos + distance, chr_len_dict[chr_name])
    return lower_bound, upper_bound

def process_chromosome_data(epi_chr, chr_len_dict):
    """Process data for a specific chromosome."""
    qtls = []
    groups = []
    processed = np.zeros(len(epi_chr), dtype=bool)  # Boolean mask for processed rows
    while not processed.all():
        lead_pair_idx = epi_chr.loc[~processed, "P"].idxmin()
        lead_pair = epi_chr.loc[lead_pair_idx]
        snp1lead_chr, snp1lead_pos = lead_pair["SNP1chr"], lead_pair["SNP1position"]
        snp2lead_chr, snp2lead_pos = lead_pair["SNP2chr"], lead_pair["SNP2position"]   
        snp1_bound = calculate_bounds(snp1lead_pos, snp1lead_chr, chr_len_dict)
        snp2_bound = calculate_bounds(snp2lead_pos, snp2lead_chr, chr_len_dict)
        candidate_mask = (
            (epi_chr["SNP1chr"] == snp1lead_chr) & (epi_chr["SNP2chr"] == snp2lead_chr) &
            (epi_chr["SNP1position"].between(*snp1_bound)) &
            (epi_chr["SNP2position"].between(*snp2_bound)) & ~processed
        )
        candidate = epi_chr[candidate_mask].copy()
        candidate.loc[:, "group"] = f"{lead_pair['SNP1']}_{lead_pair['SNP2']}_{lead_pair['BETA_INT']}_{lead_pair['P']}"
        groups.append(candidate)
        snp1_qtl_start = candidate[candidate["SNP1position"] <= snp1lead_pos]["SNP1position"].min() or snp1_bound[0]
        snp1_qtl_end = candidate[candidate["SNP1position"] >= snp1lead_pos]["SNP1position"].max() or snp1_bound[1]
        snp2_qtl_start = candidate[candidate["SNP2position"] <= snp2lead_pos]["SNP2position"].min() or snp2_bound[0]
        snp2_qtl_end = candidate[candidate["SNP2position"] >= snp2lead_pos]["SNP2position"].max() or snp2_bound[1]
        qtls.append(pd.DataFrame([{
            "SNP1CHR": snp1lead_chr, "SNP1QTL_start": snp1_qtl_start, "SNP1QTL_end": snp1_qtl_end,
            "SNP2CHR": snp2lead_chr, "SNP2QTL_start": snp2_qtl_start, "SNP2QTL_end": snp2_qtl_end,
            "lead_CHR1": snp1lead_chr, "lead_POS1": snp1lead_pos, "lead_SNP1": lead_pair["SNP1"],
            "lead_CHR2": snp2lead_chr, "lead_POS2": snp2lead_pos, "lead_SNP2": lead_pair["SNP2"],
            "lead_BETA_INT": lead_pair["BETA_INT"], "lead_P": lead_pair["P"]
        }]))
        processed |= candidate_mask
    qtls_df = pd.concat(qtls, ignore_index=True)
    groups_df = pd.concat(groups, ignore_index=True)
    return qtls_df, groups_df

#def process_group(grp, chr_len_dict):
#    return process_chromosome_data(grp, chr_len_dict)
def process_group(grp, chr_len_dict):
    try:
        qtls_df, groups_df = process_chromosome_data(grp, chr_len_dict)
        print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - Processing result for: {grp['SNP1chr'].iloc[0]} - {grp['SNP2chr'].iloc[0]}", flush=True)
        if not qtls_df.empty and not groups_df.empty:
            return qtls_df, groups_df
        else:
            print(f"Warning: empty result for: {grp['SNP1chr'].iloc[0]} - {grp['SNP2chr'].iloc[0]}", flush=True)
            return None
    except Exception as e:
        print(f"Error processing group {grp['SNP1chr'].iloc[0]} - {grp['SNP2chr'].iloc[0]}: {e}", flush=True)
        return None

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="This script was used to pick representative epistatic effect pairs.")
    parser.add_argument("--trait", help="Input trait name", type=str, default="carcass_length")
    parser.add_argument("--pvalue", help="Input epi pvalue", type=float, default=1e-14)
    parser.add_argument("--rsqure", help="Input R2", type=float, default=0.1)
    parser.add_argument("--chrLen", help="Chromosome Length file path", type=str, default="./")
    parser.add_argument("--epi", help="Epi result path", type=str, default="./")
    parser.add_argument("--output", help="Output path", type=str, default="./")
    args = parser.parse_args()

    # Load chromosome length file
    trait = args.trait
    chr_len_path = os.path.join(args.chrLen, "chr.txt")
    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}. INFO: loading chromosome length file for {trait}.....", flush=True)
    chr_len = pd.read_csv(chr_len_path, sep=" ", header=None, names=["chr", "length"])

    # Load and filter epi data
    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}. INFO: loading epi file for {trait}.....", flush=True)
    epi_path = os.path.join(args.epi, f"{trait}.epi_r2_p1df.txt")
    epi = pd.read_csv(epi_path, sep="\t", na_values="NA")
    epi.columns = ["SNP1", "SNP2", "BETA_INT", "P", "R2", "P1df"]
    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}. INFO: total {len(epi)} raw pairs for {trait}.....", flush=True)

    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}. INFO: filter epi file for {trait}.....", flush=True)
    epi = epi[(epi["P"] <= args.pvalue) & ((epi["R2"].isna()) | (epi["R2"] < args.rsqure))]
    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}. INFO: remaining {len(epi)} pairs for {trait}.....", flush=True)

    epi["SNP1chr"], epi["SNP1position"] = zip(*epi["SNP1"].str.split(":|_").map(lambda x: (x[0], int(x[1]))))
    epi["SNP2chr"], epi["SNP2position"] = zip(*epi["SNP2"].str.split(":|_").map(lambda x: (x[0], int(x[1]))))
    
    # SNP1和SNP2位置对调的pair和原pair是同一个pair，只保留其中的一个
    epi["SNP_combination"] = np.where(
        epi["SNP1position"] < epi["SNP2position"],
        epi["SNP1"] + "_" + epi["SNP2"],
        epi["SNP2"] + "_" + epi["SNP1"]
    )
    epi = epi.drop_duplicates(subset="SNP_combination").drop(columns=["SNP_combination"])
    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}. INFO: after remove dup pairs, remaining {len(epi)} pairs for {trait}.....", flush=True)
    epi.sort_values(by=["SNP1chr", "SNP1position", "SNP2chr", "SNP2position"], inplace=True)
    output_path_epi = os.path.join(args.output, f"{trait}_filterAndDedup_pair.txt")
    epi.to_csv(output_path_epi, sep=" ", index=False, na_rep="NA")

    epi.reset_index(drop=True, inplace=True)
    
    chr_len_dict = dict(zip(chr_len["chr"], chr_len["length"]))

    # Split the data by SNP1chr and process each chromosome in parallel
    epi_chrs = {chr_name: grp for chr_name, grp in epi.groupby(["SNP1chr", "SNP2chr"]) if not grp.empty}
    with ProcessPoolExecutor() as executor:
        results = executor.map(process_group, epi_chrs.values(), [chr_len_dict]*len(epi_chrs))
    
    results_list = list(results)
    # Combine results from each chromosome
    qtls_df = pd.concat([result[0] for result in results_list], ignore_index=True)
    groups_df = pd.concat([result[1] for result in results_list], ignore_index=True)

    # Sort and save results
    qtls_df.sort_values(by=["SNP1CHR", "SNP1QTL_start", "SNP1QTL_end", "SNP2CHR", "SNP2QTL_start", "SNP2QTL_end"], inplace=True)
    groups_df.sort_values(by=["SNP1chr", "SNP1position", "SNP2chr", "SNP2position"], inplace=True)

    output_path_grouped = os.path.join(args.output, f"{trait}_grouped_result.txt")
    output_path_qtls = os.path.join(args.output, f"{trait}_representative_pair.txt")

    groups_df.to_csv(output_path_grouped, sep=" ", index=False, na_rep="NA")
    qtls_df.to_csv(output_path_qtls, sep=" ", index=False, na_rep="NA")

    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}. INFO: all done for {trait}.....", flush=True)

if __name__ == "__main__":
    main()