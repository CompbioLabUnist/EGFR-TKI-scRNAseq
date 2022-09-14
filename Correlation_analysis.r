library(Seurat)
library(data.table)
library(plyr)
library(reshape)

ref.exp.file <- "./exprMatrix.tsv"
ref.meta.file <- "./meta.tsv"
seurat.file <- "./Seurat.RDS"

ref.exp <- data.frame(fread(file = ref.exp.file, sep = "\t", data.table = F), row.names = 1, check.names = F)
ref.meta <- read.table(file = ref.meta.file, sep = "\t", header = T, stringsAsFactors = F, row.names = 1)

seurat.obj <- readRDS(seurat.file)
cell.meta <- data.frame(Cells = rownames(seurat.obj@meta.data), orig.ident = seurat.obj@meta.data$orig.ident)

ref.df <- data.frame()

for (temp.celltype in unique(ref.meta$CellType))  {
    temp.cells <- rownames(subset(ref.meta, CellType == temp.celltype))

    temp.mean <- rowMeans(ref.exp[,temp.cells])
    temp.df <- data.frame(Genes = names(temp.mean), Value = as.numeric(temp.mean))
    colnames(temp.df) <- c("Genes", temp.celltype)

    if (nrow(ref.df) == 0)    {
        ref.df <- temp.df
    }else{
        ref.df <- join(ref.df, temp.df, by = "Genes")
    }
}
rownames(ref.df) <- ref.df$Genes
ref.df$Genes <- NULL

exp.df <- as.matrix(seurat.obj@assays$RNA@data)
genes <- rownames(exp.df)[rownames(exp.df) %in% rownames(ref.df)]

exp.df <- exp.df[genes,]
ref.df <- ref.df[genes,]

cor.df <- cor(exp.df, ref.df)

cor.melt <- melt(cor.df)
colnames(cor.melt) <- c("Cells", "Celltypes", "Corr")
cor.melt <- join(cor.melt, meta, by = "Cells")

cor.best.df <- data.frame()

for (temp.cell in as.character(unique(cor.melt$Cells))) {
	temp.subset.df <- subset(cor.melt, Cells == temp.cell)
	cor.best.df <- rbind(cor.best.df, temp.subset.df[temp.subset.df$Corr == max(temp.subset.df$Corr),])
}
