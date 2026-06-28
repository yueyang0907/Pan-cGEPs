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

GROUP_LEVELS <- c("Nonresponder", "Responder")
LABEL_NR     <- "Nonresponder"
LABEL_R      <- "Responder"

# =============================
# 2. 输入文件
# =============================
sample_files <- list.files(
  path = "/sibcb1/bioinformatics/yangyue/project/immunotherapy/1.new_data/",
  pattern = "_seurat_Neutrophils.rds$",
  full.names = TRUE
)

usage_data <- list.files(
  path = "/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.5.cNMF_Neutrophils/4.3.Neutrophils_GEP30_usage/",
  pattern = "_cGEP_usage_cell_level.csv$",
  full.names = TRUE
)

# 建立文件索引 (关键步骤)
usage_names <- sub(
  "_cGEP_usage_cell_level$",
  "",
  file_path_sans_ext(basename(usage_data))
)
names(usage_data) <- usage_names

# =============================
# 3. 输出目录
# =============================
outdir <- "/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.5.cNMF_Neutrophils/4.3.2.Neutrophils_usage_immu_pheatmap_GEP30/"
cell_dir <- file.path(outdir, "cell_level")
patient_dir <- file.path(outdir, "patient_level")

if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
dir.create(cell_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(patient_dir, recursive = TRUE, showWarnings = FALSE)

# =============================
# 4. 工具函数
# =============================
get_stat_value <- function(df, group_name, col) {
  v <- df[df$group == group_name, col]
  if (length(v) == 0) NA else as.numeric(v)
}

# =============================
# 5. 主循环
# =============================
stat_all <- list()

for (file in sample_files) {
  
  # 提取样本名
  sample_name <- sub("_seurat_Neutrophils$", "", file_path_sans_ext(basename(file)))
  message("------------------------------------------------")
  message(">>> Processing: ", sample_name)
  
  # 匹配 CSV 文件
  csv_file <- usage_data[sample_name]
  if (is.na(csv_file)) {
    warning("!!! Skipping: No usage file found for sample: ", sample_name)
    next
  }
  
  # 读取数据
  obj <- readRDS(file)
  usage_df <- read.csv(csv_file, row.names = 1, check.names = FALSE)
  
  # 数据转置: cell x cGEP -> cGEP x cell
  Usage_cGEP <- t(as.matrix(usage_df))
  
  # 整理 Metadata
  # 注意：确保 obj@meta.data 中包含 'sample' 和 'Response' 列
  meta <- obj@meta.data %>%
    rownames_to_column("cell_id") %>%
    transmute(
      cell_id,
      sample,
      group = factor(Response, levels = GROUP_LEVELS)
    )
  
  # 整理长表格
  Usage_cGEP_long <- as.data.frame(Usage_cGEP) %>%
    rownames_to_column("cGEP") %>%
    pivot_longer(-cGEP, names_to = "cell_id", values_to = "score") %>%
    left_join(meta, by = "cell_id") %>%
    filter(!is.na(group))
  
  # 计算病人水平数据
  Usage_cGEP_patient <- Usage_cGEP_long %>%
    group_by(sample, cGEP, group) %>%
    summarise(score_mean = mean(score, na.rm = TRUE), .groups = "drop")
  
  # 准备绘图循环
  all_cGEP <- unique(Usage_cGEP_long$cGEP)
  n_pages <- ceiling(length(all_cGEP) / PATHWAYS_PER_PAGE)
  
  for (page in seq_len(n_pages)) {
    
    subset_cGEP <- all_cGEP[
      ((page - 1) * PATHWAYS_PER_PAGE + 1):
        min(page * PATHWAYS_PER_PAGE, length(all_cGEP))
    ]
    
    plot_cell <- list()
    plot_patient <- list()
    
    for (pw in subset_cGEP) {
      
      # -------------------------------------------------
      # A. CELL LEVEL ANALYSIS
      # -------------------------------------------------
      dfc <- Usage_cGEP_long %>%
        filter(cGEP == pw)
      
      stats_c <- dfc %>%
        group_by(group) %>%
        summarise(
          n = n(),
          median = median(score, na.rm = TRUE),
          .groups = "drop"
        )
      
      cell_n_R   <- get_stat_value(stats_c, LABEL_R, "n")
      cell_n_NR  <- get_stat_value(stats_c, LABEL_NR, "n")
      cell_med_R <- get_stat_value(stats_c, LABEL_R, "median")
      cell_med_NR <- get_stat_value(stats_c, LABEL_NR, "median")
      
      cell_p <- NA
      cell_auc <- NA
      roc_c <- NULL
      
      # 至少两组且每组不为空才做检验
      if (n_distinct(dfc$group) == 2) {
        cell_p <- wilcox.test(score ~ group, data = dfc)$p.value
        roc_c <- roc(dfc$group, dfc$score, levels = GROUP_LEVELS, quiet = TRUE)
        cell_auc <- as.numeric(auc(roc_c))
      }
      
      # --- Cell Boxplot ---
      p_cell <- ggplot(dfc, aes(group, score, fill = group)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(width = 0.2, size = 0.3, alpha = 0.3) +
        theme_bw(base_size = 8) +
        theme(legend.position = "none") +
        labs(title = pw, y = "Cell Usage")
      
      if (!is.na(cell_p)) {
        p_cell <- p_cell +
          annotate("text", x = 1.5,
                   y = max(dfc$score, na.rm = TRUE),
                   label = paste0("p=", signif(cell_p, 3)))
      }
      
      # --- Cell ROC ---
      p_cell_roc <- ggplot() + theme_void()
      if (!is.na(cell_auc)) {
        p_cell_roc <- ggplot(
          data.frame(
            fpr = 1 - roc_c$specificities,
            tpr = roc_c$sensitivities
          ),
          aes(fpr, tpr)
        ) +
          geom_line(color = "red", linewidth = 0.4) +
          geom_abline(linetype = 2) +
          annotate("text", x = 0.6, y = 0.15,
                   label = paste0("AUC=", round(cell_auc, 2))) +
          theme_bw(base_size = 8)
      }
      
      plot_cell[[pw]] <- p_cell + p_cell_roc
      
      
      # -------------------------------------------------
      # B. PATIENT LEVEL ANALYSIS
      # -------------------------------------------------
      dfp <- Usage_cGEP_patient %>%
        filter(cGEP == pw)
      
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
      pat_med_NR <- get_stat_value(stats_p, LABEL_NR, "median")
      
      pat_p <- NA
      pat_auc <- NA
      roc_p <- NULL
      
      if (n_distinct(dfp$group) == 2) {
        pat_p <- wilcox.test(score_mean ~ group, data = dfp)$p.value
        roc_p <- roc(dfp$group, dfp$score_mean, levels = GROUP_LEVELS, quiet = TRUE)
        pat_auc <- as.numeric(auc(roc_p))
      }
      
      # --- Patient Boxplot ---
      p_pat <- ggplot(dfp, aes(group, score_mean, fill = group)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(width = 0.15) +
        theme_bw(base_size = 8) +
        theme(legend.position = "none") +
        labs(title = pw, y = "Patient Mean Usage")
      
      if (!is.na(pat_p)) {
        p_pat <- p_pat +
          annotate("text", x = 1.5,
                   y = max(dfp$score_mean, na.rm = TRUE),
                   label = paste0("p=", signif(pat_p, 3)))
      }
      
      # --- Patient ROC ---
      p_pat_roc <- ggplot() + theme_void()
      if (!is.na(pat_auc)) {
        p_pat_roc <- ggplot(
          data.frame(
            fpr = 1 - roc_p$specificities,
            tpr = roc_p$sensitivities
          ),
          aes(fpr, tpr)
        ) +
          geom_line(color = "blue", linewidth = 0.4) +
          geom_abline(linetype = 2) +
          annotate("text", x = 0.6, y = 0.15,
                   label = paste0("AUC=", round(pat_auc, 2))) +
          theme_bw(base_size = 8)
      }
      
      plot_patient[[pw]] <- p_pat + p_pat_roc
      
      # -------------------------------------------------
      # C. 汇总统计数据
      # -------------------------------------------------
      stat_all[[length(stat_all) + 1]] <- data.frame(
        dataset = sample_name,
        cGEP = pw,
        
        cell_N_R = cell_n_R,
        cell_N_NR = cell_n_NR,
        cell_median_R = cell_med_R,
        cell_median_NR = cell_med_NR,
        cell_pvalue = cell_p,
        cell_auc = cell_auc,
        
        patient_N_R = pat_n_R,
        patient_N_NR = pat_n_NR,
        patient_median_R = pat_med_R,
        patient_median_NR = pat_med_NR,
        patient_pvalue = pat_p,
        patient_auc = pat_auc
      )
    } # end pathway loop
    
    # 保存绘图
    ggsave(
      file.path(cell_dir, paste0(sample_name, "_cell_page", page, ".png")),
      wrap_plots(plot_cell, ncol = GEP_NCOL),
      width = FIG_WIDTH,
      height = HEIGHT_PER_GEP * ceiling(length(subset_cGEP) / GEP_NCOL),
      dpi = 300,
      limitsize = FALSE
    )
    
    ggsave(
      file.path(patient_dir, paste0(sample_name, "_patient_page", page, ".png")),
      wrap_plots(plot_patient, ncol = GEP_NCOL),
      width = FIG_WIDTH,
      height = HEIGHT_PER_GEP * ceiling(length(subset_cGEP) / GEP_NCOL),
      dpi = 300,
      limitsize = FALSE
    )
    
  } # end page loop
  
  message("<<< DONE: ", sample_name)
}

# =============================
# 6. 输出统计表
# =============================
stat_df <- bind_rows(stat_all)

write.csv(
  stat_df,
  file.path(outdir, "4.3.2.Neutrophils_cGEP30_usage_cell_patient_stats.csv"),
  row.names = FALSE
)

message("ALL FINISHED")