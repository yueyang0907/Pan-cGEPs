# Pan-cGEPs

# 1. 数据

---

## 1.1训练集(**h5ad格式**)

### CD8T

/sibcb1/bioinformatics/yangyue/project/immunotherapy/1.TISCH_data_CD8T/H5ad_CD8T

### NK

/sibcb1/bioinformatics/yangyue/project/immunotherapy/1.TISCH_data_NK

### Mono_Macrophage:

/sibcb1/bioinformatics/yangyue/project/immunotherapy/1.TISCH_data_Mono_Macro

### Neutrophils:

/sibcb1/bioinformatics/yangyue/project/immunotherapy/1.TISCH_data_Neutrophils

## 1.2免疫治疗数据集

### 癌症的免疫治疗响应数据(rds格式)

/sibcb1/bioinformatics/yangyue/project/immunotherapy/1.new_data

### 癌症的免疫治疗响应数据(h5ad格式)

/sibcb1/bioinformatics/yangyue/project/immunotherapy/1.new_data_h5ad/CD8T

/sibcb1/bioinformatics/yangyue/project/immunotherapy/1.new_data_h5ad/NK

/sibcb1/bioinformatics/yangyue/project/immunotherapy/1.new_data_h5ad/Macro

/sibcb1/bioinformatics/yangyue/project/immunotherapy/1.new_data_h5ad/Neutrophils

汇总：

/sibcb1/bioinformatics/yangyue/project/immunotherapy/1.new_data_h5ad_v2

### 外周血的免疫治疗响应数据(rds格式)

/sibcb1/bioinformatics/yangyue/project/immunotherapy/1.new_data/pbmc

### 外周血的免疫治疗响应数据(h5ad格式)

/sibcb1/bioinformatics/yangyue/project/immunotherapy/1.new_data_h5ad/PBMC

## 1.3 非肿瘤免疫疾病数据

COVID:

/sibcb1/bioinformatics/yangyue/project/241120_TFs/data/COVID_GSE145926_allcell.h5ad

SLE:

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.2.cNMF_CD8T/5.startCAT-SLE/GSE135779_PBMC_SLE_Adult_celltype_harmony_umap.h5ad

organoids  transplant **：**

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.2.cNMF_CD8T/5.3_startCAT-organoids/GSE248788/GSE248789_scRNA_clean_anno.h5ad

Aging:

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.2.cNMF_CD8T/5.5_startCAT_Aging/GSE157007_scRNA_anno.h5ad

# 2. CD8T cGEPs

---

## 2.1 CD8T-run cNMF：

输入:

/sibcb1/bioinformatics/yangyue/project/immunotherapy/1.TISCH_data_CD8T/H5ad_CD8T/

代码：

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.2.cNMF_CD8T/1.cNMF_CD8T.py

log:

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.2.cNMF_CD8T/1.cNMF_CD8T.log

结果：

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.2.cNMF_CD8T/cnmf_CD8T_out

## 2.2 merge CD8T GEPs to cGEPs

输入：

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.2.cNMF_CD8T/cnmf_CD8T_out

代码：

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.2.cNMF_CD8T/2.cluster_cnmf_results_CD8T.py

log：

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.2.cNMF_CD8T/2.cluster_cnmf_results_CD8T.log

输出：

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.2.cNMF_CD8T/Example_refbuilder_95_v2

## 2.3 CD8T-注释

### 2.3.1 CD8T免疫学标志物注释：

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.2.cNMF_CD8T/3.1.cNMF_CD8T_StateAnno.ipynb

### 2.3.2 CD8T-使用了TCAT文章中的52种类型

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.2.cNMF_CD8T/3.2.cNMF_CD8T_FunctionAnno.ipynb

对应的结果：

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.2.cNMF_CD8T/3.cGEP_topgene/test/cGEP1_70_function_CD8_Spectra52_annotation.csv

### 2.3.3 CD8T-三种通路富集分析:

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.2.cNMF_CD8T/3.2.CD8T_GEP_PATHWAY.ipynb

pathway结果：

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.2.cNMF_CD8T/3.2.CD8T_GEP_pathway/3.CD8T_Cluster_Pathway_Master_Annotation.xlsx 

## 2.4 CD8T-最终注释结果

综合上述免疫学标志物与三种富集分析方法，（cGEPs名称，marker, 类别）：

包含非编码基因：

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.2.cNMF_CD8T/3.cGEP_topgene/3.3.CD8T_cGEP_Anno_Complete_With_Genes-allgene.csv

过滤非编码基因:

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.2.cNMF_CD8T/3.cGEP_topgene/3.3.CD8T_cGEP_Anno_Complete_With_Genes.csv

## 2.5 CD8T-cGEPs-QC

### 2.5.1 癌肿在cGEPs中分布情况

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.2.cNMF_CD8T/3.4.2.QC_CD8T_癌肿样本名.ipynb

结果：

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.2.cNMF_CD8T/cnmf_CD8T_out_Sample_counts/3.4.CD8TcGEP_in_cancer.png

### 2.5.2 cGEP的代表性marker绘制热图:

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.2.cNMF_CD8T/3.5.2.CD8T.ExampleMarkGene.ipynb

结果
/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.2.cNMF_CD8T/3.5.CD8T.QC.MarkerGene/3.5.CD8T_Markers_Functional.png

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.2.cNMF_CD8T/3.5.CD8T.QC.MarkerGene/3.5.CD8T_Markers_Lineage.png

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.2.cNMF_CD8T/3.5.CD8T.QC.MarkerGene/3.5.CD8T_Markers_Doublet_and_Artifact.png

# 3.NK_cGEPs

---

## 3.1 NK-run cNMF

输入：

/sibcb1/bioinformatics/yangyue/project/immunotherapy/1.TISCH_data_NK/

代码：

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.3.cNMF_NK/1.cNMF_NK.py

log：

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.3.cNMF_NK/1.cNMF_NK.log

输出：

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.3.cNMF_NK/cnmf_NK_out/

## 3.2 merge NK GEPs to cGEPs

代码：

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.3.cNMF_NK/2.cluster_cnmf_results_NK.py

对应的log：

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.3.cNMF_NK/2.cluster_cnmf_results_NK_v2.log

输出的结果文件夹：

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.3.cNMF_NK/Example_refbuilder_NK_v2

## 3.3 NK 注释

免疫学标志物注释：

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.3.cNMF_NK/3.1.NK_Anno.ipynb

免疫学标志物注释结果：

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.3.cNMF_NK/3.cGEP_Anno_NK/3.NK_cGEP_Annotation_Results.csv

三种富集分析注释:

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.3.cNMF_NK/3.2.NK_GEP_PATHWAY与整合注释.ipynb

三种富集分析结果:

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.3.cNMF_NK/3.2.NK_GEP_PATHWAY/3.2.NK_Cluster_Pathway_Master_Annotation.xlsx 

## **3.4 NK-最终注释结果**

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.3.cNMF_NK/3.3.NK_GEP_最终表格整理.ipynb

包含非编码基因：

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.3.cNMF_NK/3.3.NK_cGEP_ALL_ANNO/3.3.NK_cGEP_Anno_Complete_With_Genes_final-allgene.csv

过滤非编码基因：

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.3.cNMF_NK/3.3.NK_cGEP_ALL_ANNO/3.3.NK_cGEP_Anno_Complete_With_Genes_final.csv

## 3.5 NK-cGEPs-QC

### 3.5.1 癌肿在cGEPs中分布情况

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.3.cNMF_NK/3.4.2.QC_NK_癌肿样本名.ipynb

### 3.5.2 cGEP的代表性marker绘制热图:

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.3.cNMF_NK/3.5.2.NK.ExampleMarkGene.ipynb

# 4. Mono-Macro_cGEPs

---

## 4.1 Mono-Macro_run cNMF

输入:

/sibcb1/bioinformatics/yangyue/project/immunotherapy/1.TISCH_data_Mono_Macro/

代码：

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.4.cNMF_Mono_Macro/1.cNMF_Mono_Macro.py

log:

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.4.cNMF_Mono_Macro/1.cNMF_Mono_Macro.log

输出:

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.4.cNMF_Mono_Macro/cnmf_Mono_Macro_out/

## 4.2 merge Mono-Macro GEPs

代码：

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.4.cNMF_Mono_Macro/2.cluster_cnmf_results_Macro.py

log:

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.4.cNMF_Mono_Macro/2.cluster_cnmf_results_Macro.log

结果：

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.4.cNMF_Mono_Macro/Example_refbuilder_Macro

## 4.3 Mono-Macro cGEPs注释

免疫学标志物：

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.4.cNMF_Mono_Macro/3.1.Macro_Anno.ipynb

三种富集分析：

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.4.cNMF_Mono_Macro/3.2.Macro_GEP_PATHWAY与整合注释.ipynb

三种富集分析结果：

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.4.cNMF_Mono_Macro/3.2.Macro_GEP_pathway/3.3.Macro_Cluster_Pathway_Master_Annotation.xlsx

## 4.4 Mono-Macro最终注释结果

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.4.cNMF_Mono_Macro/3.3.Macro_GEP_最终表格整理.ipynb

包含非编码基因：

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.4.cNMF_Mono_Macro/3.1.Macro_GEP_Anno/3.3.Macro_cGEP_Anno_Complete_With_Genes-allgene.csv

过滤非编码基因:

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.4.cNMF_Mono_Macro/3.1.Macro_GEP_Anno/3.3.Macro_cGEP_Anno_Complete_With_Genes.csv

## 4.5 Mono-Macro cGEPs-QC

### 4.5.1 癌肿在cGEPs中的分布情况

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.4.cNMF_Mono_Macro/3.4.2.QC_Mono_Macro_癌肿样本名.ipynb

### 4.5.2  cGEP的代表性marker绘制热图:

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.4.cNMF_Mono_Macro/3.5.2.Mono_Macro_ExampleMarkGene.ipynb

# 5. Neutrophils_cGEPs

---

## 5.1 Neutrophils run cNMF

输入:

/sibcb1/bioinformatics/yangyue/project/immunotherapy/1.TISCH_data_Neutrophils/

代码：

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.5.cNMF_Neutrophils/1.cNMF_Neutrophils.py

log：

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.5.cNMF_Neutrophils/1.cNMF_Neutrophils.log

结果:

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.5.cNMF_Neutrophils/cnmf_Neutrophils_out/

## 5.2 Neutrophils merge GEPs

输入：

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.5.cNMF_Neutrophils/cnmf_Neutrophils_out

代码：

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.5.cNMF_Neutrophils/2.cluster_cnmf_results_Neutrophils.py

log：

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.5.cNMF_Neutrophils/2.cluster_cnmf_results_Neutrophils.log

输出：

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.5.cNMF_Neutrophils/Example_refbuilder_Neutrophils

## 5.3 Neutrophils cGEPs注释

免疫学标志物注释：

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.5.cNMF_Neutrophils/3.1.Neu_Anno.ipynb

免疫学标志物注释结果：

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.5.cNMF_Neutrophils/3.Neutrophils_GEP_Anno/3.Neutrophil_cGEP30_Annotation_Results.csv

三种富集分析注释：

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.5.cNMF_Neutrophils/3.2.2.Neu_GEP30_pathway.ipynb

三种富集分析结果：

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.5.cNMF_Neutrophils/3.2.Neu_pathway/3.2.2.Neutrophils_GEP30_Cluster_Pathway_Master_Annotation.xlsx 

## 5.4 Neutrophils cGEPs 最终注释结果

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.5.cNMF_Neutrophils/3.3.Neutrophil_GEP_最终表格整理.ipynb

包含非编码基因：

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.5.cNMF_Neutrophils/3.Neutrophils_GEP_Anno/3.3.Neutrophil_cGEP30_Anno_Complete_With_Genes-allgene.csv

过滤非编码基因：

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.5.cNMF_Neutrophils/3.Neutrophils_GEP_Anno/3.3.Neutrophil_cGEP30_Anno_Complete_With_Genes.csv

## 5.5 Neutrophils cGEPs QC

### 5.5.1 癌肿在cGEPs中的分布情况

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.5.cNMF_Neutrophils/3.4.2.QC_Neutrophils_癌肿样本名.ipynb

### 5.5.2  cGEP的代表性marker绘制热图:

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.5.cNMF_Neutrophils/3.5.2.Neutrophils_ExampleMarkGene.ipynb

# 6. 应用一：预测肿瘤免疫治疗响应

---

## 6.1 CD8T

---

代码1（cGEPs在不同肿瘤免疫治疗数据中的活性得分）：

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.2.cNMF_CD8T/4.3.CD8T_GEP_usage.ipynb

代码2（取病人水平）:

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.2.cNMF_CD8T/4.3.2.CD8T_GEP_uasge_immu_pheatmap.r  

代码2对应的log文件:

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.2.cNMF_CD8T/4.3.2.CD8T_GEP_uasge_immu_pheatmap.log 

可视化（绘制热图）：

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.2.cNMF_CD8T/4.3.3.CD8T_GEP70_uasge_pheatmap.ipynb

## 6.2 NK

---

代码1（cGEPs在不同肿瘤免疫治疗数据中的活性得分）：

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.3.cNMF_NK/4.3.NK_GEP_usage.ipynb

代码2（取病人水平）:

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.3.cNMF_NK/4.3.2.NK_GEP_uasge_immu_pheatmap.r

代码2对应的log文件:

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.3.cNMF_NK/4.3.2.NK_GEP_uasge_immu_pheatmap.log

可视化（绘制热图）：

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.3.cNMF_NK/4.3.3.NK_GEP41_uasge_pheatmap-0310.ipynb

## 6.3 Mono-Macro

---

代码1（cGEPs在不同肿瘤免疫治疗数据中的活性得分）：

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.4.cNMF_Mono_Macro/4.3.Macro_GEP_usage.ipynb

代码2（取病人水平）:

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.4.cNMF_Mono_Macro/4.3.2.Macro_GEP_uasge_immu_pheatmap.r

代码2对应的log文件:

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.4.cNMF_Mono_Macro/4.3.2.Macro_GEP_uasge_immu_pheatmap.log

可视化（绘制热图）：

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.4.cNMF_Mono_Macro/4.3.3.Mono_Macro_GEP77_uasge_pheatmap-0310.ipynb

## 6.4 Neutrophils

---

代码1（cGEPs在不同肿瘤免疫治疗数据中的活性得分）：

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.5.cNMF_Neutrophils/4.3.Neutrophils_GEP_usage.ipynb

代码2（取病人水平）:

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.5.cNMF_Neutrophils/4.3.2.Neutrophils_GEP_uasge_immu_pheatmap.r

代码2对应的log文件:

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.5.cNMF_Neutrophils/4.3.2.Neutrophils_GEP_uasge_immu_pheatmap.log

可视化（绘制热图）：

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.5.cNMF_Neutrophils/4.3.3.Neutrophils_GEP30_uasge_pheatmap-0310.ipynb

# 7. 应用二：非肿瘤相关的免疫疾病中cGEPs的表达差异

---

## 7.1 CD8T

---

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.2.cNMF_CD8T/5.disease/5.1.CD8T_COVID_GSE145926.ipynb

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.2.cNMF_CD8T/5.disease/5.2.CD8T_SLE.ipynb

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.2.cNMF_CD8T/5.disease/5.3.CD8T_organoids.ipynb

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.2.cNMF_CD8T/5.disease/5.4.CD8T_Aging.ipynb

对应的结果：

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.2.cNMF_CD8T/5.disease/5.1.COVID

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.2.cNMF_CD8T/5.disease/5.2.SLE

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.2.cNMF_CD8T/5.disease/5.3.organ

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.2.cNMF_CD8T/5.disease/5.4.Aging

## 7.2 NK

---

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.3.cNMF_NK/5.1.NK_COVID_GSE145926.ipynb

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.3.cNMF_NK/5.2.NK_SLE_GSE135779.ipynb

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.3.cNMF_NK/5.3.NK_organoids_GSE248789.ipynb

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.3.cNMF_NK/5.4.NK_Aging.ipynb

对应的结果：

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.3.cNMF_NK/5.disease/5.1.COVID

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.3.cNMF_NK/5.disease/5.2.SLE

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.3.cNMF_NK/5.disease/5.3.organoids

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.3.cNMF_NK/5.disease/5.4.Aging

## 7.3 Mono-Macro

---

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.4.cNMF_Mono_Macro/5.disease/5.1.Macro_COVID.ipynb

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.4.cNMF_Mono_Macro/5.disease/5.2.Macro_SLE.ipynb

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.4.cNMF_Mono_Macro/5.disease/5.3.Macro_oganoids.ipynb

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.4.cNMF_Mono_Macro/5.disease/5.4.Macro_Aging.ipynb

对应结果：

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.4.cNMF_Mono_Macro/5.disease/5.1.COVID

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.4.cNMF_Mono_Macro/5.disease/5.2.SLE

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.4.cNMF_Mono_Macro/5.disease/5.3.organ

/sibcb1/bioinformatics/yangyue/project/immunotherapy/7.4.cNMF_Mono_Macro/5.disease/5.4.Aging