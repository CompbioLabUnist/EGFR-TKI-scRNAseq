library(Seurat)
library(infercnv)
library(plyr)

epithelial.file <- "Seurat.epithelial.RDS"
myeloid.file <- "Seurat.meyloid.RDS"


epithelial.obj <- readRDS(epithelial.file)
myeloid.obj <- readRDS(myeloid.file)

seurat.obj <- merge(epithelial.obj, myeloid.obj)

count.df <- seurat.obj@assays$RNA@counts
celltype.df <- data.frame(id = rownames(seurat.obj@meta.data), celltype = seurat.obj@meta.data$celltypes)

write.table(x = as.matrix(count.df), file = "./InferCNV.input.matrix", quote = F, sep = "\t")
write.table(x = celltype.df, file = "./InferCNV.annotation", quote = F, sep = "\t", row.names = F, col.names = F)


infercnv.obj <- CreateInfercnvObject(raw_counts_matrix = "./InferCNV.input.matrix", gene_order_file = "./genes.gtf", annotations_file = "./InferCNV.annotation",
									 ref_group_names = c("Monocyte"), 
									 delim = "\t", chr_exclude = c("Y", "MT"))

infercnv.obj <- infercnv::run(infercnv_obj = infercnv.obj, cutoff = 0.1, cluster_by_groups = T, denoise = T, HMM = F, out_dir = "./")