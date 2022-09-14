library(Seurat)
library(plyr)
library(DESeq2)

seurat.obj <- readRDS("Seurat.epithelial.RDS")

cpm.df <- data.frame()

for (temp.sample in unique(seurat.obj@meta.data$orig.ident)) {
    temp.subset <- subset(seurat.obj, orig.ident == temp.sample)

    if (nrow(temp.subset@meta.data) >= 1)   {
        temp.cpm <- rowMeans(as.matrix(RelativeCounts(data = temp.subset@assays$RNA@counts, scale.factor = 1e6)))
        temp.cpm <- data.frame(Genes = names(temp.cpm), CPM = as.numeric(temp.cpm))
        colnames(temp.cpm)[2] <- temp.sample

        if(ncol(x.mat) == 0)    {
            cpm.df <- temp.cpm
        }else{
            cpm.df <- join(cpm.df, temp.cpm, by = "Genes")
        }
    }
}
cpm.df <- cpm.df[rowSums(cpm.df) > 10,]

group.df <- data.frame(do.call(rbind, strsplit(x = colnames(cpm.df), split = "\\.")))
colnames(group.df) <- c("Groups", "Patients")

coldata <- data.frame(Conditions = group.df$Groups)
coldata$Conditions <- factor(x = coldata$Conditions, levels = c("SE", "PR"))

dds <- DESeqDataSetFromMatrix(countData = cpm.df, colData = coldata, design = ~Conditions)
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds)
res <- results(dds)
