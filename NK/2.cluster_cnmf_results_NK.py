import sys
import os
import re
import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from starcat import BuildConsensusReference, starCAT

root = "/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.3.cNMF_NK/cnmf_NK_out"

cnmf_paths = []
Ks = []

for dataset in os.listdir(root):
    dataset_dir = os.path.join(root, dataset)
    if not os.path.isdir(dataset_dir):
        continue
    for sub in os.listdir(dataset_dir):
        # 匹配 K*_final
        m = re.match(r"K(\d+)_final", sub)
        if m is None:
            continue
        K = int(m.group(1))
        k_final_dir = os.path.join(dataset_dir, sub)
        # 找 *_cNMF_K*_final 目录
        for leaf in os.listdir(k_final_dir):
            if leaf.endswith(f"_cNMF_K{K}_final"):
                cnmf_paths.append(os.path.join(k_final_dir, leaf))
                Ks.append(K)

dts = [0.01] * len(cnmf_paths)
outdir = '/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.3.cNMF_NK/Example_refbuilder_NK_v2'
prefix = 'starcat_ref'
refbuilder = BuildConsensusReference(cnmf_paths, ks=Ks, density_thresholds=dts,output_dir=outdir, prefix=prefix)

clus_df, spectra_tpm_grouped, spectra_scores_grouped, hvgs_union, top_genes = refbuilder.cluster_cnmf_results()


outdir = refbuilder.output_dir
topgene_dir = os.path.join(outdir, "top_genes_by_cGEP_NK")
os.makedirs(topgene_dir, exist_ok=True)

for cgep, df in top_genes.items():
    out_file = os.path.join(topgene_dir, f"{cgep}_top_genes.csv")
    df.to_csv(out_file)
    print("✅ saved:", out_file)

