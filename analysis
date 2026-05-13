library(tximport)
library(dplyr)
library(ggplot2)
library(DESeq2)
library(readxl)
library(ggthemes)
library(rtracklayer)
library(tibble)
library(ggrepel)
library(EnhancedVolcano)
library(pheatmap)
library(RColorBrewer)
library(clusterProfiler)
library(enrichplot)
library(wordcloud)
library(cowplot)
library(ggpubr)
library(ggthemes)
library(org.Hs.eg.db)
library(scales)
library(stringr)
library(rrvgo)

###### variables

work_dir <- "~/path/to/project/"
counts_dir <- "raw_data/"
meta_file <- "metadata/sample_metadata.csv"
gtf_path <- "ref/gencode.annotation.gtf.gz"

outdir <- "results/"
plotdir <- "plots/"

#design params
design_factor <- "Condition" 
trt_level <- "Treated"
ctrl_level <- "Vehicle" 

## thresholds
fdr_cutoff <- 0.05
lfc_cutoff <- 1.0


setwd(work_dir)

dir.create(outdir, showWarnings = FALSE)
dir.create(plotdir, showWarnings = FALSE)

files <- list.files(path = paste0(getwd(), "/", counts_dir), pattern = ".sf", full.names = TRUE, recursive = TRUE)

folder_names <- basename(dirname(files))
folder_names <- gsub("\\_.*","",folder_names)
names(files) <- folder_names


gtf <- rtracklayer::import(gtf_path)
gtf_df <- as.data.frame(gtf)

tx2gene <- gtf_df %>%
  filter(type == "transcript") %>%
  dplyr::select(transcript_id, gene_id)

gene_name_map <- gtf_df %>% 
  filter(type == "gene") %>% 
  dplyr::select(gene_id, gene_name) 

txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene)

##load meta
#meta <- read_xlsx(meta_file)
meta <- read.csv(meta_file, stringsAsFactors = FALSE)
rownames(meta) <- meta$Sample_ID # adjust 'Sample_ID' to whatever your col is called

#ensure sample order matches exactly between counts and metadata
sample_order <- colnames(txi.salmon$counts)
meta <- meta[match(sample_order, rownames(meta)), ]


add_gene_name <- function(res_obj, gene_map = gene_name_map) {
  df <- as.data.frame(res_obj) %>%
    rownames_to_column(var = "gene_id") %>%
    left_join(gene_map, by = "gene_id") %>%
    arrange(padj, abs(log2FoldChange))
  return(df)
}

#deseq2 analysis

meta[[design_factor]] <- as.factor(meta[[design_factor]])
meta[[design_factor]] <- relevel(meta[[design_factor]], ref = ctrl_level)

design_formula <- as.formula(paste("~", design_factor))

dds <- DESeqDataSetFromTximport(txi.salmon, meta, design_formula)

keep <- rowSums(counts(dds)) > 10
dds <- dds[keep, ]

dds <- DESeq(dds)
#dds <- DESeq(dds, parallel=T)

plotDispEsts(dds)

res <- results(dds, contrast = c(design_factor, trt_level, ctrl_level))
res <- lfcShrink(dds, contrast = c(design_factor, trt_level, ctrl_level), res=res, type="ashr")

#quick pvalue check
res %>% as.data.frame() %>%
  arrange(padj) %>%
  ggplot(aes(x = pvalue)) +
  geom_histogram(color = "white", bins = 50) +
  ggtitle(paste("P-value histogram:", trt_level, "vs", ctrl_level)) -> pval_hist
print(pval_hist)


res_df <- add_gene_name(res)
sig_genes_df <- res_df %>% filter(!is.na(padj) & padj < fdr_cutoff)

write.table(res_df, paste0(outdir, "All_DEGs_", trt_level, "_vs_", ctrl_level, ".txt"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(sig_genes_df, paste0(outdir, "Sig_DEGs_", trt_level, "_vs_", ctrl_level, ".txt"), sep = "\t", row.names = FALSE, quote = FALSE)

### PCA

vsdata <- vst(dds, blind = TRUE)

png(filename = paste0(plotdir,"PCA_plot.png"), res = 600, width = 18, height = 18, units = "cm")
plotPCA(vsdata, intgroup = c(design_factor)) + 
  theme_clean() +
  labs(color = "Condition")
dev.off()

#heatmap of sample distances

rld <- rlog(dds, blind = TRUE)
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld[[design_factor]], rownames(meta), sep = "_")
colnames(sampleDistMatrix) <- rownames(sampleDistMatrix)
colours <- colorRampPalette(rev(brewer.pal(9, "Greens")))(255)

png(filename = paste0(plotdir,"Sample_Distance_Heatmap.png"), res = 600, width = 34, height = 18, units = "cm")
pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colours, fontsize_col = 8, fontsize_row = 8)
dev.off()

#### export normalised counts for external plotting
norm_counts <- as.data.frame(counts(dds, normalized = TRUE))
norm_counts$gene_id <- rownames(norm_counts)
norm_counts <- norm_counts %>% left_join(gene_name_map, by = "gene_id") %>% relocate(gene_name, .before = everything())
write.csv(norm_counts, file = paste0(outdir,"normalised_counts.csv"), row.names = FALSE)


#Plots

res_df <- res_df[!is.na(res_df$padj), ]
sig_vars <- sum(res_df$padj < fdr_cutoff & abs(res_df$log2FoldChange) > lfc_cutoff, na.rm = TRUE)
my_subtitle <- paste0("Total ", nrow(res_df), " genes. ", sig_vars, " were significant.")

png(filename = paste0(plotdir,"Volcano_plot.png"), res = 600, width = 28, height = 16, units = "cm")
EnhancedVolcano(res_df,
                lab = res_df$gene_name,          
                x = 'log2FoldChange',
                y = 'padj',                    
                title = paste0(trt_level, ' vs ', ctrl_level),
                pCutoff = fdr_cutoff,                
                FCcutoff = lfc_cutoff,                  
                pointSize = 3, selectLab = NA,
                legendPosition = "right", subtitle = "", caption = my_subtitle)
dev.off()

# top 50 variable genes heatmap
vsd <- vst(dds, blind=FALSE)
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing=TRUE), 50)
mat  <- assay(vsd)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)

annotation_col <- as.data.frame(colData(vsd)[, design_factor, drop = FALSE]) 

png(filename = paste0(plotdir,"Top50VariableGenes.png"), res = 600, width = 28, height = 16, units = "cm")
pheatmap(mat, 
         annotation_col = annotation_col, 
         show_rownames = FALSE,  
         show_colnames = FALSE,
         fontsize_row = 8,
         main = "Top 50 Variable Genes", name="Expression")
dev.off()


###### GO analysis

# Upregulated
up_genes <- res_df$gene_name[res_df$padj < fdr_cutoff & res_df$log2FoldChange > lfc_cutoff]

ego_up <- enrichGO(gene          = up_genes, 
                   universe      = res_df$gene_name,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = 'SYMBOL', 
                   ont           = "BP",              
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05)

ego_up2 <- simplify(ego_up, cutoff=0.5, by="p.adjust", select_fun=min)

plot_df_up <- as.data.frame(ego_up2) %>% slice_head(n = 10)

g_up <- ggplot(plot_df_up, aes(x = -log10(p.adjust), 
                               y = reorder(Description, -log10(p.adjust)),
                               fill = -log10(p.adjust))) +
  geom_col(width = 0.8, color = "black") +
  scale_fill_gradient(low = "#d2e9f7", high = "#08519c") +
  scale_y_discrete(labels = label_wrap(40)) +
  ylab("") + xlab("-log10(p.adjust)") +
  theme_pubr() +
  ggtitle("Upregulated GO Terms") +
  theme(
    text = element_text(size = 14, colour = "black"),
    legend.position = "right"
  )

ggsave(paste0(plotdir,"enrichGO_upregulated.png"), plot = g_up, width = 14, height = 6, dpi = 300)

# Downregulated
down_genes <- res_df$gene_name[res_df$padj < fdr_cutoff & res_df$log2FoldChange < -lfc_cutoff]

ego_down <- enrichGO(gene          = down_genes, 
                     universe      = res_df$gene_name,
                     OrgDb         = org.Hs.eg.db,
                     keyType       = 'SYMBOL', 
                     ont           = "BP",              
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.05)

ego_down2 <- simplify(ego_down, cutoff=0.5, by="p.adjust", select_fun=min)

plot_df_down <- as.data.frame(ego_down2) %>% slice_head(n = 10)

g_down <- ggplot(plot_df_down, aes(x = -log10(p.adjust), 
                                   y = reorder(Description, -log10(p.adjust)),
                                   fill = -log10(p.adjust))) +
  geom_col(width = 0.8, color = "black") +
  scale_fill_gradient(low = "#d2e9f7", high = "#08519c") +
  scale_y_discrete(labels = label_wrap(40)) +
  ylab("") + xlab("-log10(p.adjust)") +
  theme_pubr() +
  ggtitle("Downregulated GO Terms") +
  theme(
    text = element_text(size = 14, colour = "black"),
    legend.position = "right"
  )

ggsave(paste0(plotdir,"enrichGO_downregulated.png"), plot = g_down, width = 14, height = 6, dpi = 300)

save.image("objects/analysis_complete.RData")
