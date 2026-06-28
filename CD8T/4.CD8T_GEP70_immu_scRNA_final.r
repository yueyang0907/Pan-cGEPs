suppressPackageStartupMessages({
  library(Seurat)
  library(CelliD)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(tidyr)
  library(tibble)
  library(tools)
  library(pROC)
})

# =============================
# 1. 参数设置
# =============================
PATHWAYS_PER_PAGE <- 20
GEP_NCOL          <- 4
HEIGHT_PER_GEP    <- 3.8
FIG_WIDTH         <- 16

# 定义分组名称 (请确保这与Seurat对象metadata中Response列的内容一致)
GROUP_LEVELS <- c("Nonresponder", "Responder")
# 定义在统计表中对应的标签 (用于提取数据)
LABEL_NR <- "Nonresponder"
LABEL_R  <- "Responder"

# =============================
# 2. 输入文件与路径
# =============================
# 样本文件路径
sample_files <- list.files(
  path = "/sibcb1/bioinformatics/yangyue/project/immunotherapy/1.new_data/",
  pattern = "_seurat_CD8T.rds$",
  full.names = TRUE
)

# 读取基因集
GEP70_top30 <- read.csv("/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.2.cNMF_CD8T/3.cGEP_topgene/ALL_cGEP_top_genes_long_wide_by_cGEP1_70.csv",stringsAsFactors = FALSE,check.names = FALSE)

# 转换基因集格式
cluster_genes <- lapply(as.list(GEP70_top30), function(x) {
  x <- as.character(x)
  x[!is.na(x) & x != ""]
})

# =============================
# 3. 输出目录设置
# =============================
outdir <- "/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.2.cNMF_CD8T/4.CD8T_GEP70_immu_scRNA"
cell_dir <- file.path(outdir, "cell_level")
patient_dir <- file.path(outdir, "patient_level")

if(!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
dir.create(cell_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(patient_dir, recursive = TRUE, showWarnings = FALSE)

# =============================
# 4. 辅助函数
# =============================
# 安全提取统计量，防止因某组缺失导致报错
get_stat_value <- function(df_stat, group_name, col_name) {
  val <- df_stat[df_stat$group == group_name, col_name]
  if (length(val) == 0) return(NA)
  return(as.numeric(val))
}

# =============================
# 5. 主循环处理
# =============================
stat_all <- list()

for (file in sample_files) {

  sample_name <- file_path_sans_ext(basename(file))
  message(">>> Processing: ", sample_name)

  # 读取与预处理
  obj <- readRDS(file)
  obj <- NormalizeData(obj)
  obj <- RunMCA(obj, assay = "RNA")

  # 运行 CelliD HGT
  HGT <- RunCellHGT(
    obj,
    pathways = cluster_genes,
    dims = 1:50,
    n.features = 200
  )

  # 整理 Metadata
  # 注意：这里确保 Response 列转换为因子，顺序对应 GROUP_LEVELS
  meta <- obj@meta.data %>%
    rownames_to_column("cell_id") %>%
    mutate(group = factor(Response, levels = GROUP_LEVELS))

  # 整理 HGT 结果 (Long Format)
  HGT_long <- as.data.frame(as.matrix(HGT)) %>%
    rownames_to_column("pathway") %>%
    pivot_longer(-pathway, names_to = "cell_id", values_to = "score") %>%
    left_join(meta, by = "cell_id")

  # -----------------------------
  # 计算病人水平均值
  # -----------------------------
  HGT_patient <- HGT_long %>%
    group_by(sample, pathway) %>%
    summarise(
      score_mean = mean(score, na.rm = TRUE),
      group = first(group),
      .groups = "drop"
    )

  # -----------------------------
  # 分页绘图与统计
  # -----------------------------
  all_pathways <- unique(HGT_long$pathway)
  n_pages <- ceiling(length(all_pathways) / PATHWAYS_PER_PAGE)

  for (page in seq_len(n_pages)) {

    pw_subset <- all_pathways[
      ((page - 1) * PATHWAYS_PER_PAGE + 1):
        min(page * PATHWAYS_PER_PAGE, length(all_pathways))
    ]

    plot_cell <- list()
    plot_patient <- list()

    for (pw in pw_subset) {

      # =========================
      # A. CELL LEVEL 分析
      # =========================
      dfc <- HGT_long %>%
        filter(pathway == pw, !is.na(group)) %>%
        droplevels()

      # --- 统计量计算 ---
      stats_c <- dfc %>%
        group_by(group) %>%
        summarise(
          n = n(),
          median = median(score, na.rm = TRUE),
          .groups = "drop"
        )

      # 提取 R (Responder) 和 NR (Nonresponder)
      cell_n_R   <- get_stat_value(stats_c, LABEL_R, "n")
      cell_n_NR  <- get_stat_value(stats_c, LABEL_NR, "n")
      cell_med_R <- get_stat_value(stats_c, LABEL_R, "median")
      cell_med_NR<- get_stat_value(stats_c, LABEL_NR, "median")

      # --- 差异检验与 ROC ---
      cell_p <- NA
      cell_auc <- NA
      roc_c <- NULL

      if (n_distinct(dfc$group) == 2) {
        # Wilcoxon 检验
        cell_p <- wilcox.test(score ~ group, data = dfc)$p.value
        # ROC 分析
        roc_c <- roc(dfc$group, dfc$score, levels = GROUP_LEVELS, quiet = TRUE)
        cell_auc <- as.numeric(auc(roc_c))
      }

      # --- 绘图 Cell ---
      p_cell <- ggplot(dfc, aes(group, score, fill = group)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(size = 0.2, alpha = 0.3, width = 0.2) +
        theme_bw(base_size = 8) +
        theme(legend.position = "none") +
        labs(title = pw, y = "Cell HGT")

      if (!is.na(cell_p)) {
        p_cell <- p_cell + annotate("text", x = 1.5,
                                    y = max(dfc$score, na.rm = TRUE),
                                    label = paste0("p=", signif(cell_p, 3)))
      }

      p_cell_roc <- ggplot() + theme_void()
      if (!is.na(cell_auc)) {
        p_cell_roc <- ggplot(data.frame(f = 1 - roc_c$specificities, t = roc_c$sensitivities), aes(f, t)) +
          geom_line(color = "red") + geom_abline(linetype = 2) +
          annotate("text", x = 0.6, y = 0.15, label = paste0("AUC=", round(cell_auc, 2))) +
          theme_bw(base_size = 8)
      }
      plot_cell[[pw]] <- p_cell + p_cell_roc


      # =========================
      # B. PATIENT LEVEL 分析
      # =========================
      dfp <- HGT_patient %>%
        filter(pathway == pw, !is.na(group)) %>%
        droplevels()

      # --- 统计量计算 ---
      stats_p <- dfp %>%
        group_by(group) %>%
        summarise(
          n = n(),
          median = median(score_mean, na.rm = TRUE),
          .groups = "drop"
        )

      pat_n_R   <- get_stat_value(stats_p, LABEL_R, "n")
      pat_n_NR  <- get_stat_value(stats_p, LABEL_NR, "n")
      pat_med_R <- get_stat_value(stats_p, LABEL_R, "median")
      pat_med_NR<- get_stat_value(stats_p, LABEL_NR, "median")

      # --- 差异检验与 ROC ---
      pat_p <- NA
      pat_auc <- NA
      roc_p <- NULL

      if (n_distinct(dfp$group) == 2) {
        pat_p <- wilcox.test(score_mean ~ group, data = dfp)$p.value
        roc_p <- roc(dfp$group, dfp$score_mean, levels = GROUP_LEVELS, quiet = TRUE)
        pat_auc <- as.numeric(auc(roc_p))
      }

      # --- 绘图 Patient ---
      p_pat <- ggplot(dfp, aes(group, score_mean, fill = group)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(width = 0.15) +
        theme_bw(base_size = 8) +
        theme(legend.position = "none") +
        labs(title = pw, y = "Patient mean HGT")

      if (!is.na(pat_p)) {
        p_pat <- p_pat + annotate("text", x = 1.5,
                                  y = max(dfp$score_mean, na.rm = TRUE),
                                  label = paste0("p=", signif(pat_p, 3)))
      }

      p_pat_roc <- ggplot() + theme_void()
      if (!is.na(pat_auc)) {
        p_pat_roc <- ggplot(data.frame(f = 1 - roc_p$specificities, t = roc_p$sensitivities), aes(f, t)) +
          geom_line(color = "blue") + geom_abline(linetype = 2) +
          annotate("text", x = 0.6, y = 0.15, label = paste0("AUC=", round(pat_auc, 2))) +
          theme_bw(base_size = 8)
      }
      plot_patient[[pw]] <- p_pat + p_pat_roc

      # =========================
      # C. 保存统计结果到列表
      # =========================
      stat_all[[length(stat_all) + 1]] <- data.frame(
        dataset = sample_name,
        pathway = pw,

        # Cell Level Stats
        cell_N_R = cell_n_R,          # R组 细胞数
        cell_N_NR = cell_n_NR,        # NR组 细胞数
        cell_Median_R = cell_med_R,   # R组 评分中位数
        cell_Median_NR = cell_med_NR, # NR组 评分中位数
        cell_pvalue = cell_p,
        cell_auc = cell_auc,

        # Patient Level Stats
        patient_R = pat_n_R,           # R组 病人数
        patient_NR = pat_n_NR,         # NR组 病人数
        patient_Median_R = pat_med_R,    # R组 评分中位数
        patient_Median_NR = pat_med_NR,  # NR组 评分中位数
        patient_pvalue = pat_p,
        patient_auc = pat_auc
      )
    } # end pathway loop

    # 保存图片
    ggsave(
      file.path(cell_dir, paste0(sample_name, "_cell_page", page, ".png")),
      wrap_plots(plot_cell, ncol = GEP_NCOL),
      width = FIG_WIDTH,
      height = HEIGHT_PER_GEP * ceiling(length(pw_subset) / GEP_NCOL),
      dpi = 300, limitsize = FALSE
    )

    ggsave(
      file.path(patient_dir, paste0(sample_name, "_patient_page", page, ".png")),
      wrap_plots(plot_patient, ncol = GEP_NCOL),
      width = FIG_WIDTH,
      height = HEIGHT_PER_GEP * ceiling(length(pw_subset) / GEP_NCOL),
      dpi = 300, limitsize = FALSE
    )

  } # end page loop

  message("<<< DONE: ", sample_name)
}

# =============================
# 6. 输出最终统计表
# =============================
stat_df <- bind_rows(stat_all)

write.csv(
  stat_df,
  file.path(outdir, "CD8T_GEP70_cell_patient_pvalue_auc.csv"),
  row.names = FALSE
)

message("ALL FINISHED. Statistics saved to GEP70_cell_patient_stats_full.csv")