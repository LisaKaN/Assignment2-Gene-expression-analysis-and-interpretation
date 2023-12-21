########################################################
#               1.) Preprocessing                      #
########################################################


# Untar folder
#folder_name = "brca_tcga_pan_can_atlas_2018.tar.gz"

#untar(folder_name)

# Go to new path
new_dir = paste(getwd(),"brca_tcga_pan_can_atlas_2018", sep = "/" )
setwd(new_dir)

# inspecting clinical file
clinical = read.delim("data_clinical_patient.txt")
rnaseq = read.delim("data_mrna_seq_v2_rsem.txt")
cna = read.delim("data_cna.txt")

# Preprocessing rnaseq data
keep = !duplicated(rnaseq[,1])
rnaseq = rnaseq[keep,]
rownames(rnaseq)  = rnaseq[,1]


rna_cna_sub = as.matrix(rnaseq[,-c(1,2)]) # col 1&2 gene names

# extract ERBB2 data 
erbb2_idx = which(cna[,1] == "ERBB2")#

# Histogram of ERBB2 expression
png("../ERBB2.png")
hist(as.numeric(cna[erbb2_idx,-c(1,2)]), main="Histogram of ERBB2 expression",
     xlab = "CNA level for ERBB2")
dev.off()

# Select rna cases which have cna data
rna_cna_id = which(is.element(colnames(rnaseq[,-c(1,2)]), colnames(cna[,-c(1,2)])))
rna_cna_sub = rnaseq[,2+rna_cna_id]


# Build metadata
meta_erbb2 = matrix(0,length(rna_cna_id),1)
colnames(meta_erbb2) = "ERBB2Amp"


# Store for which patients erbb2 is amplified
for (i in 1:length(rna_cna_id)){
  col_i = colnames(rna_cna_sub)[i]
  col_cna = which(colnames(cna)==col_i)
  meta_erbb2[i,] = 1*(cna[erbb2_idx,col_cna]>0)
}


#pos_example = which(meta_erbb2==1)[1]

#col_i = colnames(rna_cna_sub)[pos_example]
#col_cna = which(colnames(cna)==col_i)



########################################################
#    2.) Differential Gene Expression Analysis         #
########################################################

# BiocManager
if(!require("BiocManager", quietly = TRUE))
 install.packages("BiocManager")
 
# Install DeSeq2
#BiocManager::install("DESeq2")
library(DESeq2)
 
# Build DESeq Object
rna_cna_sub[is.na(rna_cna_sub)] = 0 
rna_cna_sub[rna_cna_sub<0] = 0
 
dds = DESeqDataSetFromMatrix(countData = round(rna_cna_sub),
                            colData = meta_erbb2,
                            design = ~ ERBB2Amp)

# Filter
smallestGroupSize = 3
keep = rowSums(counts(dds) >= 10) >= smallestGroupSize
dds = dds[keep,]

# Normalize
dds = DESeq(dds)

# Get Results
res = results(dds)
res_05 = results(dds, alpha = 0.05)
 
# Summary
summary(res)
summary(res_05)
rownames(res) = rnaseq[keep, 1]
 
# Significantly Differentially Expressed
signif = which(res$padj<0.05)
deg = res[signif,]
 
# Separate into downregulated and upregulated
dup = deg[deg[,2]>0.,]
ddown = deg[deg[,2]<0.,]
print('Number of overexpressed genes:')
print(length(dup$log2FoldChange))
print('Number of underexpressed genes:')
print(length(ddown$log2FoldChange))

# Get top ten genes in dup and ddown by log2FoldChange
topten_up = dup[order(-dup$log2FoldChange),][1:10,]
topten_down = ddown[order(ddown$log2FoldChange),][1:10,]
 
# Visualize top differentially expressed genes
## up-regulated
library(ggplot2)
df_up = data.frame(topten_up)
p_up = ggplot(data=df_up, aes(x=reorder(rownames(df_up), -log2FoldChange), y=log2FoldChange)) +
        geom_bar(stat="identity")
p_up

## down-regulated
df_down = data.frame(topten_down)
p_down = ggplot(data=df_down, aes(x=reorder(rownames(df_down), -log2FoldChange), y=log2FoldChange)) +
  geom_bar(stat="identity")
p_down

## up and down
df_res = rbind(df_up, df_down)
p_all = ggplot(data=df_res, aes(x=reorder(rownames(df_res), -log2FoldChange), y=log2FoldChange)) +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
  xlab("Gene") +
  ylab("log2 Fold Change") +
  ggtitle("Top over- and underexpressed genes")
ggsave("../genes.png")
p_all


########################################################
#               3.) Pathway Enrichment                 #
########################################################

# Getting entrez ids of significantly differentially expressed genes
entrez_ids = rnaseq[keep,2]

entrez_all = entrez_ids[signif]
entrez_up = entrez_all[signif[deg[,2]>0.]]
entrez_down = entrez_all[signif[deg[,2]<0.]]

#BiocManager::install("clusterProfiler")
library(clusterProfiler)
 
# KEGG pathway over-representation analysis
all_paths = enrichKEGG(gene = entrez_all, organism = 'hsa', pvalueCutoff = 0.05)
head(all_paths)

# Visualize top pathways
top_pathways = head(all_paths[order(-all_paths$Count),], n=10)
p_pathways = ggplot(data=data.frame(top_pathways), aes(x=reorder(Description, Count), y=Count)) +
  geom_bar(stat="identity") + 
  coord_flip() +
  xlab("Pathway description") +
  ggtitle("Top overrepresented pathways") +
  theme(axis.text = element_text(size = 7))
ggsave("../top_pathways.png")
p_pathways

# Visualization signaling pathway "Pathways of neurodegeneration - multiple diseases"

# browseKEGG(all_paths, 'hsa05022')

BiocManager::install("pathview")
library(pathview)

neurodeg = pathview(gene.data  = entrez_all,
                     pathway.id = 'hsa05022',
                     species    = 'hsa',
                     limit      = list(gene=max(abs(entrez_all)), cpd=1))

########################################################
#               4.) PCA                                #
########################################################

# Transform to visualize
rld = vst(dds, blind = FALSE)

# Do PCA
pc = prcomp(t(assay(rld)))

# PCA plot
png("../PCA.png")
plotPCA(rld, intgroup="ERBB2Amp")
dev.off()





