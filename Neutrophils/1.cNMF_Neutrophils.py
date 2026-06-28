import os
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from cnmf import cNMF

# =============================
# 路径配置（Mono / Macro）
# =============================
h5ad_dir = "/sibcb1/bioinformatics/yangyue/project/immunotherapy/1.TISCH_data_Neutrophils/"
qc_dir   = "/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.5.cNMF_Neutrophils/h5ad_Neutrophils_QC/"
out_dir  = "/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.5.cNMF_Neutrophils/cnmf_Neutrophils_out/"

os.makedirs(qc_dir, exist_ok=True)
os.makedirs(out_dir, exist_ok=True)

# =============================
# cNMF 参数
# =============================
Ks = np.arange(5, 11)
n_iter = 20
final_n_iter = 200
seed = 14
density_threshold = 0.01
num_highvar_genes = 2000

# =============================
# 批量处理
# =============================
for h5ad_file in sorted(os.listdir(h5ad_dir)):

    if not h5ad_file.endswith(".h5ad"):
        continue

    dataset = h5ad_file.replace(".h5ad", "")
    print(f"\n========== {dataset} ==========")

    # -----------------------------
    # 1. 读取 & QC
    # -----------------------------
    adata = sc.read_h5ad(os.path.join(h5ad_dir, h5ad_file))

    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()

    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_cells(adata, min_counts=200)
    sc.pp.filter_genes(adata, min_cells=3)

    if adata.n_obs < 200:
        print("⚠️ Too few cells, skip dataset")
        continue

    qc_file = os.path.join(qc_dir, f"{dataset}_cnmfinput.h5ad")
    adata.write(qc_file)

    # -----------------------------
    # 2. exploratory cNMF（多 K）
    # -----------------------------
    ds_out = os.path.join(out_dir, dataset)
    os.makedirs(ds_out, exist_ok=True)

    cnmf_obj = cNMF(
        output_dir=ds_out,
        name=f"{dataset}_cNMF_explore"
    )

    cnmf_obj.prepare(
        counts_fn=qc_file,
        components=Ks,
        n_iter=n_iter,
        seed=seed,
        num_highvar_genes=num_highvar_genes
    )

    cnmf_obj.factorize(worker_i=0, total_workers=1)
    cnmf_obj.combine()

    # -----------------------------
    # 3. K* 选择（stability）
    # -----------------------------
    stability_scores = {}

    for k in Ks:
        try:
            cnmf_obj.consensus(
                k=k,
                density_threshold=density_threshold,
                show_clustering=False
            )

            usage, spectra_scores, spectra_tpm, top_genes = cnmf_obj.load_results(
                K=k,
                density_threshold=density_threshold
            )

            R = np.corrcoef(spectra_scores.T)
            iu = np.triu_indices_from(R, k=1)
            stability = np.nanmedian(R[iu])

            stability_scores[k] = stability
            print(f"K={k}, stability={stability:.4f}")

        except Exception as e:
            print(f"K={k} failed: {e}")

    if len(stability_scores) == 0:
        print("❌ No valid K, skip dataset")
        continue

    best_k = max(stability_scores, key=stability_scores.get)
    print(f"⭐ Selected K* = {best_k}")

    # -----------------------------
    # 4. 保存 K-selection 曲线
    # -----------------------------
    plt.figure(figsize=(5, 4))
    plt.plot(
        list(stability_scores.keys()),
        list(stability_scores.values()),
        marker="o"
    )
    plt.axvline(best_k, color="red", linestyle="--", label=f"K*={best_k}")
    plt.xlabel("K (number of programs)")
    plt.ylabel("Stability (median pairwise correlation)")
    plt.title(f"{dataset} cNMF K selection (Mono_Macro)")
    plt.legend()
    plt.tight_layout()
    plt.savefig(
        os.path.join(ds_out, f"{dataset}_K_selection_curve.png"),
        dpi=300
    )
    plt.close()

    # -----------------------------
    # 5. exploratory consensus
    # -----------------------------
    cnmf_obj.consensus(
        k=best_k,
        density_threshold=density_threshold,
        show_clustering=True,
        close_clustergram_fig=False
    )

    plt.savefig(
        os.path.join(ds_out, f"{dataset}_consensus_k{best_k}.png"),
        dpi=300,
        bbox_inches="tight"
    )
    plt.close()

    # -----------------------------
    # 6. Final cNMF（K*，200 iter）
    # -----------------------------
    final_out = os.path.join(ds_out, f"K{best_k}_final")
    os.makedirs(final_out, exist_ok=True)

    cnmf_final = cNMF(
        output_dir=final_out,
        name=f"{dataset}_cNMF_K{best_k}_final"
    )

    cnmf_final.prepare(
        counts_fn=qc_file,
        components=[best_k],
        n_iter=final_n_iter,
        seed=seed,
        num_highvar_genes=num_highvar_genes
    )

    cnmf_final.factorize(worker_i=0, total_workers=1)
    cnmf_final.combine()

    cnmf_final.consensus(
        k=best_k,
        density_threshold=density_threshold,
        show_clustering=False
    )

    # -----------------------------
    # 7. 保存 final 结果
    # -----------------------------
    usage_f, spectra_scores_f, spectra_tpm_f, top_genes_f = cnmf_final.load_results(
        K=best_k,
        density_threshold=density_threshold
    )

    pd.DataFrame(
        usage_f,
        index=adata.obs_names,
        columns=[f"GEP{i+1}" for i in range(usage_f.shape[1])]
    ).to_csv(os.path.join(final_out, f"{dataset}_usage_K{best_k}_final.csv"))

    spectra_scores_f.to_csv(
        os.path.join(final_out, f"{dataset}_spectra_scores_K{best_k}_final.csv")
    )
    spectra_tpm_f.to_csv(
        os.path.join(final_out, f"{dataset}_spectra_tpm_K{best_k}_final.csv")
    )
    top_genes_f.to_csv(
        os.path.join(final_out, f"{dataset}_top_genes_K{best_k}_final.csv")
    )

    print(f"🎯 Final cNMF finished for {dataset} (Mono_Macro, K*={best_k}, 200 iters)")


