rm(list = ls())
setwd("d://Datas/Transcriptome/GSE279704/")

# Step1 数据预处理-------------------------------------------------------------------
library(data.table)
files <- list.files("GSE279704_RAW/", pattern = "*.txt.gz", full.names = TRUE)

exp_list <- lapply(files, function(file) {
  data <- fread(file)
  data_clean <- data[, 1:2]
  return(data_clean)
})
exp_combined <- Reduce(function(x, y) merge(x, y, by = "Gene_Symbol", all = TRUE), exp_list)

exp_combined


library(clusterProfiler)
library(org.Mm.eg.db)
library(AnnotationDbi)
all_genes <- keys(org.Mm.eg.db, keytype = "ENTREZID")
gene_types <- AnnotationDbi::select(org.Mm.eg.db, 
                                    keys = all_genes, 
                                    columns = c("SYMBOL", "GENETYPE"),  # 不用 GENETYPE
                                    keytype = "ENTREZID")

table(is.na(gene_types$GENETYPE))
table(gene_types$GENETYPE)

gene_types = gene_types %>% 
  filter(GENETYPE == "protein-coding")

exp_mrna = exp_combined %>% 
  filter(Gene_Symbol %in% gene_types$SYMBOL)


# 临床文件下载
library(GEOquery)
options(timeout = 3000)
gse = "GSE279704"
eSet <- getGEO(gse,
               destdir = '.',
               getGPL = F)
pd <- pData(eSet[[1]])
gpl <- eSet[[1]]@annotation



library(tidyverse)
pd1 = pd %>% 
  select(geo_accession, title)

pd1$group= c(rep("Control",3),rep("Levo",3), rep("MPTP",3), rep("BTK_High",3), rep("BTK_Low",3))


rownames(exp_mrna) = exp_mrna$Gene_Symbol
exp_mrna = as.data.frame(exp_mrna)

exp = exp_mrna %>% 
  select(-Gene_Symbol)

colnames(exp) <- str_remove(colnames(exp), "_Read_Count")
colnames(exp) <- pd$geo_accession[match(colnames(exp), pd1$title)]





# 过滤低质基因
nonzero_genes <- rowSums(exp != 0) > 0
exp_count <- exp[nonzero_genes, ]
# 基因至少在30%的样本中有表达，根据实际情况调整
exp_count = exp_count[apply(exp_count, 1, function(x) sum(x > 0) > 0.3*ncol(exp_count)), ]
nrow(exp_count)

# filter_count <- apply(exp, MARGIN = 1, FUN = function(x){
#   sum(x > 1) >= ncol(exp) * 0.3
# }) 
# exp_counts <- exp[filter_count, ]



# rm(list = ls())
# load("step1_output.rdata")
### 获取基因长度
# Count矩阵转换
# <https://mp.weixin.qq.com/s/QlcVet7nyxeZUlPFFw07TA>
library(GenomicFeatures)

## 1.读取GTF文件
txdb <- makeTxDbFromGFF("d://Chromosomal gene annotation/gencode.vM36.annotation.gtf.gz")
#TxDb：R语言中用于储存gtf文件的一种格式

## 2.提取外显子信息
exonic <- exonsBy(txdb, by="gene") 
#可以提取基因外显子部分，计算counts数据比对的为基因外显子序列

## 3.外显子长度求和
exonic.gene.sizes <- sum(width(reduce(exonic))) 
#reduce可以删除外显子重复序列部分

## 4.结果整理
mygeneids <- data.frame(gene_id=names(exonic.gene.sizes),width=exonic.gene.sizes)

library(org.Mm.eg.db) #人org.Hs.eg.db
library(AnnotationDbi)
mygeneids$gene_id <- gsub("\\..*", "",mygeneids$gene_id) #去除版本号
mygeneids$gene_symbol <- mapIds(org.Mm.eg.db,
                                keys=mygeneids$gene_id,
                                column="SYMBOL",
                                keytype="ENSEMBL",
                                multiVals="first") #名称转换由ENSEMBL转换为SYMBOL
rownames(mygeneids) <- c(1:nrow(mygeneids))

library(dplyr)
gene_length <- mygeneids %>% 
  dplyr::distinct(gene_symbol,.keep_all = T) %>% 
  dplyr:: filter(!is.na(gene_symbol)) %>%
  dplyr::select(-gene_id)

gene_length <- cbind(gene_length$gene_symbol,gene_length$width)
colnames(gene_length) <- c("gene_symbol","width")
gene_length <- data.frame(gene_length)



# load("step1_output.rdata")
count1 = exp_count
count1$gene_symbol = rownames(count1) 
gene_length_both = gene_length %>% 
  filter(gene_symbol %in% count1$gene_symbol)

count1 = count1 %>% 
  filter(gene_symbol %in% gene_length$gene_symbol)

count2 = count1 %>% 
  left_join(gene_length_both, count1, by = "gene_symbol")

rownames(count2) = count2$gene_symbol
count2 = count2 %>% 
  dplyr::select(-gene_symbol)
# 
# # Count转CPM
# cpm <- apply(X =subset(count2, select = c(-width)), 
#              MARGIN =2, 
#              FUN =function(x){
#                x/sum(as.numeric(x)) * 10^6
#              })


# # Count转RPKM\FPKM
# geneLengths <- as.vector(subset(count2, select = c(width))) 
# geneLengths <- as.numeric(count2$width)
# # 创建基因长度的向量，为下一步标准化做准备  
# rpkm <- apply(X = subset(count2, select = c(-width)),
#               MARGIN = 2, 
#               FUN = function(x) {
#                 10^9 * x / geneLengths / sum(as.numeric(x))
#               }) # 计算rpkm


# Count转TPM
geneLengths <- as.vector(subset(count2, select = c(width))) 
geneLengths <- as.numeric(count2$width)
rpk <- apply( subset(count2, select = c(-width)), 2, 
              function(x) x/(geneLengths/1000)) #求基因长度归一化值
tpm <- apply(X = rpk, 
             MARGIN = 2, 
             FUN = function(x){
               x / sum(as.numeric(x)) * 10^6
             })#使用 rpk 值根据样本大小进行标准化
# exp_count = exp_filtered
exp_tpm = tpm

save(exp_count,exp_tpm, pd1,file = "step1_output.rdata")


# Step2 差异分析-------------------------------------------------------------------
rm(list = ls())
load("step1_output.rdata")
exp_count[] <- lapply(exp_count, as.numeric)

library(tidyverse)
pd2 <- pd1 %>% filter(group %in% c('Control',"MPTP"))

exp_count_filtered <- exp_count[,colnames(exp_count) %in% pd2$geo_accession]
exp_tpm_filtered <- exp_tpm[,colnames(exp_tpm) %in% pd2$geo_accession]

library(DESeq2)
library(tidyverse)
group_list = pd2$group

# 设置条件变量
condition = factor(group_list,levels = c("Control", "MPTP"))
coldata <- data.frame(row.names = colnames(exp_count_filtered), condition)
# 创建 DESeq2 数据集对象
dds <- DESeqDataSetFromMatrix(countData = exp_count_filtered, colData = coldata, design = ~  condition)
# 进行差异分析
dds <- DESeq(dds)   
result <- as.data.frame(results(dds))
result$genes <- rownames(result)
# 根据 log2FoldChange 和 pvalue 计算基因的变化情况
log2FC_t= 1.5
change = ifelse(result$pvalue>0.05,'Stable',
                ifelse(abs(result$log2FoldChange) < log2FC_t,'Stable',
                       ifelse(result$log2FoldChange >= log2FC_t,'Up','Down') ))
table(change)
result <- result %>% 
  mutate(result, change) %>% 
  arrange(desc(log2FoldChange))


# 统计各类变化的基因数目
table(result$change)
table(is.na(result$change))

result <- na.omit(result)
table(is.na(result$change))

# heatmap
fivenum(result$log2FoldChange)
fivenum(-log10(result$pvalue))

library(ggplot2)
volcano_plot = ggplot(result, aes(x = log2FoldChange, 
                                  y = -log10(pvalue))) + 
  geom_point(size = 2.5, 
             alpha = 0.5, 
             aes(color = change), 
             show.legend = TRUE) + 
  scale_color_manual(values = c('#00468BFF','gray','#ED0000FF'))+
  ylim(0, 30)+
  xlim(-10, 10)+
  labs(x = 'Log2FoldChange',y = '-Log10(P.Value)')+
  geom_hline(yintercept = -log10(0.05),
             linetype = 2,
             color = 'black',lwd = 0.5)+
  geom_vline(xintercept = c(-1.5, 1.5),
             linetype = 2, 
             color = 'black', lwd = 0.5)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
volcano_plot
ggsave(filename = "volcano_plot.pdf",dpi = 300,
       device = "pdf", width = 6,height = 6)
# write_excel_csv(result,file = "result.csv")

# Heatmap
#略

save(result, exp_count_filtered,group_list, exp_tpm_filtered,pd2, file = "step2_output.rdata")


####### PCA
rm(list = ls())
load("step2_output.rdata")
library(FactoMineR)
library(factoextra)
library(tidyverse)
library(ggsci)
cors <- pal_jama()(5)
exp_pca = log2(as.data.frame(exp_tpm_filtered) + 1)

dat=as.data.frame(t(exp_pca))
dat.pca <- PCA(dat, scale.unit = F, graph = FALSE)
pca_plot <- fviz_pca_ind(dat.pca,
                         geom.ind = "point", #使用点来表示样本
                         col.ind = group_list, #根据分组进行着色
                         palette = cors,#####色彩颜色根据分组个数决定
                         addEllipses = TRUE, #添加椭圆，以表示各组样本的分布区域
                         legend.title = "Groups") #图例标题
print(pca_plot)


# Step3 富集分析 ------------------------------------------------------------------
rm(list = ls())
load("step2_output.rdata")

library(clusterProfiler)
library(org.Mm.eg.db)
library(tidyverse)
library(enrichplot)

degs_result = result %>% 
  filter(change %in% c("Up",'Down'))


diff_entrez<-bitr(
  degs_result$genes,
  fromType = 'SYMBOL',
  toType = 'ENTREZID',
  OrgDb = 'org.Mm.eg.db'
)


go_enrich<-clusterProfiler::enrichGO(gene = diff_entrez$SYMBOL,
                                     ont = 'all',#可选'BP','CC','MF' or 'all'
                                     keyType = "SYMBOL",
                                     OrgDb = org.Mm.eg.db,
                                     pAdjustMethod = "BH",#p值矫正方法
                                     pvalueCutoff = 0.05,
                                     qvalueCutoff = 0.05)

go_geo<- clusterProfiler::simplify(go_enrich, 
                                   cutoff=0.7, 
                                   by="pvalue",
                                   select_fun=min)

go_result<-go_geo@result
# write_excel_csv(file = "go_result.csv", go_result)
head(go_result)

# upsetplot(go_geo)

eGoBP <- go_result %>%  
  filter(ONTOLOGY=="BP") %>% 
  arrange(desc(FoldEnrichment)) %>% 
  head(10)

eGoCC <- go_result %>%  
  filter(ONTOLOGY=="CC") %>% 
  arrange(desc(FoldEnrichment)) %>% 
  head(10)

eGoMF <- go_result %>%  
  filter(ONTOLOGY=="MF") %>% 
  arrange(desc(FoldEnrichment)) %>%
  head(10)

eGo10 <- rbind(eGoBP,eGoMF,eGoCC)
eGo10$Log10Pvalue = -log10(eGo10$pvalue)


library(ggplot2)
library(ggthemes)
GO_p <- ggplot(eGo10,aes(FoldEnrichment,
                         fct_reorder(factor(Description),FoldEnrichment)
                         )) +  
  geom_point(aes(size = Count,color = -Log10Pvalue)) +  
  scale_color_gradient(low = "#ED000099",high = "#00468BFF") + 
  scale_y_discrete(position = 'right') +
  labs(color = "-log10(Pvalue)",size="Count", shape="Ontology",       
       x = "Enrichment Factor",
       y = "GO term",
       title="GO enrichment") +  
  theme_test() + 
  theme(legend.position = "left")
GO_p

p1 = GO_p + facet_wrap( ~ ONTOLOGY)
p1
p2 = GO_p + facet_wrap( ~ ONTOLOGY,
                        ncol= 1,scale='free')
p2

ggsave("Go_dot_ggplot.pdf", plot = p2, device = "pdf", 
       width = 9,  height = 9, 
       units = "in",dpi = 300)


# KEGG
KEGG_enrich <- clusterProfiler::enrichKEGG(gene = diff_entrez$ENTREZID,
                                           organism = "mmu", #物种Homo sapiens 
                                           pvalueCutoff = 0.05,#pvalue阈值
                                           qvalueCutoff = 0.05,#qvalue阈值
                                           pAdjustMethod = "BH",
                                           minGSSize = 10,#富集分析中考虑的最小基因集合大小
                                           maxGSSize = 500)#富集中考虑的最大基因集合大小
KEGG_enrich<-setReadable(KEGG_enrich,
                         OrgDb = org.Mm.eg.db,
                         keyType = 'ENTREZID')
KEGG_result<-KEGG_enrich@result
head(KEGG_result)


KEGG_result1 = KEGG_result %>% 
  filter(pvalue < 0.05) %>% 
  arrange(desc(FoldEnrichment)) 

# write_excel_csv(file = "KEGG_result.csv", KEGG_result1)

KEGG_p <- ggplot(KEGG_result1,aes(FoldEnrichment,fct_reorder(factor(Description), FoldEnrichment))) +  
  geom_point(aes(size=Count,color= -log10(pvalue))) +  
  scale_color_gradient(low="#ED000099",high = "#00468BFF") +  
  labs(color="-log10(Pvalue)",size="Count",       
       x="Enrichment Factor",
       y="KEGG term",
       title="KEGG enrichment") +  
  theme_test() 
KEGG_p

ggsave("KEGG_dot_ggplot.pdf", plot = KEGG_p, device = "pdf", 
       width = 9,  height = 9, 
       units = "in",dpi = 300)


# Step4 GSEA分析 ------------------------------------------------------------
rm(list = ls())
load("step2_output.rdata")

library(clusterProfiler)
library(org.Mm.eg.db)
library(tidyverse)
library(enrichplot)

entrez<-bitr(
  result$genes,
  fromType = 'SYMBOL',
  toType = 'ENTREZID',
  OrgDb = 'org.Mm.eg.db'
)

gene_list <- result$log2FoldChange
names(gene_list) <- rownames(result)

gene_list <- gene_list[names(gene_list) %in% entrez[,1]]
names(gene_list) <- entrez[match(names(gene_list), entrez[, 1]), 2]
length(gene_list)
head(gene_list)


KEGG_gse <- gseKEGG(geneList = gene_list, 
                    organism = "mmu", 
                    minGSSize = 10, 
                    maxGSSize = 500, 
                    pvalueCutoff = 0.05, 
                    pAdjustMethod = "BH", 
                    verbose = FALSE, 
                    eps = 0)

KEGG_gse <- setReadable(KEGG_gse, 
                        OrgDb = org.Mm.eg.db, 
                        keyType = "ENTREZID")

KEGG_gse_result <- KEGG_gse@result

KEGG_gse_result <- KEGG_gse_result %>%
  arrange(desc(abs(NES)))

# write_excel_csv(file = "KEGG_gse_result.csv", KEGG_gse_result)

fivenum(KEGG_gse_result$NES)
p_gsea = ggplot(KEGG_gse_result, 
           aes(x = NES, y = reorder(Description, NES))) +
  geom_point(aes(size=setSize,color=-log10(pvalue))) +  
  scale_color_gradient(low="#ED000099",high = "#00468BFF") +  
  xlim(-2,2) + 
  labs(color="-log10(Pvalue)",size="setSize",       
       x="Normalized Enrichment Score (NES)",
       y="GSEA KEGG term",
       title="GSEA KEGG enrichment") +  
  theme_test() 
p_gsea

ggsave("KEGG_gsea_ggplot.pdf", plot = p_gsea, device = "pdf", 
       width = 6,  height = 6, 
       units = "in",dpi = 300)

# ridgeplot(KEGG_gse)
# 
# gseaplot2(KEGG_gse, 
#           geneSetID = 1, 
#           title = KEGG_gse$Description[1])
# 
# gseaplot2(KEGG_gse, 
#           geneSetID = 2, 
#           title = KEGG_gse$Description[2])
# 
# gseaNb(object = KEGG_gse,
#        geneSetID = 'mmu04080')
# 
# gseaNb(object = KEGG_gse,
#        geneSetID = 'mmu03040')


library(GseaVis)

gsea_plot1 = gseaNb(object = KEGG_gse,
       geneSetID = 'mmu04080', 
       addPval = T, 
       pvalX =  0.95,
       pvalY =  0.85,
       subPlot = 2)
gsea_plot1
ggsave(filename = "gsea_plot1.pdf",plot = gsea_plot1, 
       device = 'pdf', dpi = 300, width = 6, height = 6)

gsea_plot2 = gseaNb(object = KEGG_gse,
       geneSetID = 'mmu03040', 
       addPval = T, 
       pvalX =  0.65,
       pvalY =  0.55,
       subPlot = 2)
gsea_plot2
ggsave(filename = "gsea_plot2.pdf",plot = gsea_plot2, 
       device = 'pdf', dpi = 300, width = 6, height = 6)


# Step5 GSVA分析 ------------------------------------------------------------
rm(list = ls())
library(tidyverse)
load("step2_output.rdata")
library(tidyverse)
library(clusterProfiler)
library(msigdbr) 
library(GSVA)
library(GSEABase)

# gene_set <- getGmt("D://Msigdb/mmu/mh.all.v2024.1.Mm.symbols.gmt")
gene_set <- getGmt("D://Msigdb/mmu/m2.all.v2024.1.Mm.symbols.gmt")

exp_gsva = log2(exp_tpm_filtered + 1)

params <- gsvaParam(as.matrix(exp_gsva), gene_set, 
                    minSize = 1, maxSize = Inf, 
                    kcdf = "Gaussian", tau = 1, 
                    maxDiff = TRUE, 
                    absRanking = FALSE)
gsva_result <- gsva(param = params, 
                    verbose = TRUE, 
                    BPPARAM = BiocParallel::SerialParam(progressbar = TRUE))

## GSVA差异
# identical(colnames(gsva_result),colnames(exp_gsva))
# identical(colnames(gsva_result),rownames(pd2))
library(limma)
Group <- factor(pd2$group, levels = c("Control", "MPTP"))
design <- model.matrix(~0 + Group)  

fit <- lmFit(gsva_result, design)
fit <- eBayes(fit)
GSVA_degs <- topTable(fit, coef = 2, number = Inf)

fivenum(GSVA_degs$logFC)
fivenum(-log10(GSVA_degs$P.Value))
fivenum(-log10(GSVA_degs$adj.P.Val))

# write_excel_csv(file = "GSVA_degs_M2.csv", GSVA_degs)
# 一般取0.5
GSVA_degs_filter = GSVA_degs %>% 
  filter(P.Value < 0.05 & abs(logFC)>0.5 )



# 热图
gsva_result_filter = gsva_result[rownames(gsva_result) %in% rownames(GSVA_degs_filter),]



library(pheatmap)
annotation_cols <- data.frame(Group = group_list)
rownames(annotation_cols) <- colnames(gsva_result)
annotation_colors = list(
  Group = c(Control = "#4500AC7F", MPTP = "#FF88887F")
)


pheatmap::pheatmap(gsva_result_filter,
                   cluster_rows = T,
                   cluster_cols = F,
                   show_colnames = F,
                   clustering_method =  "ward.D2",
                   display_numbers = T,
                   main = "Correlation Heatmap",
                   annotation_col = annotation_cols,
                   annotation_colors = annotation_colors,
                   color = colorRampPalette(c('#00468BFF', "white", '#ED000099'))(50),
                   width = 20,height = 12, 
                   filename = "gsva_mmu_M2.pdf")


# Step6 通路 ------------------------------------------------------------
rm(list = ls())
library(data.table)
library(tidyverse)
degs = fread("result.csv")

degs_filter = degs %>% 
  filter(change != "Stable")

# OS通路
OS_genes = clusterProfiler::read.gmt("d://Msigdb/mmu/BIOCARTA_ARENRF2_PATHWAY.v2024.1.Mm.gmt")

os_both = intersect(degs_filter$genes, OS_genes$gene)
# 无
# 
# 
# # Hal os
# Hal_OS_genes = clusterProfiler::read.gmt("d://Msigdb/mmu/HALLMARK_OXIDATIVE_PHOSPHORYLATION.v2024.1.Mm.gmt")
# hal_both = intersect(degs_filter$genes, Hal_OS_genes$gene)
# # 无
# 
# 
# # Cal
# cal_genes = clusterProfiler::read.gmt("d://Msigdb/mmu/REACTOME_CALCITONIN_LIKE_LIGAND_RECEPTORS.v2024.1.Mm.gmt")
# cal_os_both = intersect(degs$genes, cal_genes$gene)
# # 无
# 
# 
# # Ace
# Ace_genes = clusterProfiler::read.gmt("d://Msigdb/mmu/WP_ACE_INHIBITOR_PATHWAY.v2024.1.Mm.gmt")
# ace_both = intersect(degs_filter$genes, Ace_genes$gene)
# # 无
# 
# 
# # Beta ox
# Box_genes = clusterProfiler::read.gmt("d://Msigdb/mmu/REACTOME_BETA_OXIDATION_OF_HEXANOYL_COA_TO_BUTANOYL_COA.v2024.1.Mm.gmt")
# Box_both = intersect(degs_filter$genes, Box_genes$gene)
# # 无

