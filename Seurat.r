
library(Seurat)

seurat.files <- Sys.glob("./*.qc.scater.Rout")
doublet.file <- "./scrublet.doublet.result"

seurat.obj.list <- list()

for (temp.file in seurat.files) {
  temp.id <- gsub(pattern = ".qc.scater.Rout", replacement = "", x = basename(temp.file))
  temp.id <- gsub(pattern = "\\_", replacement = "\\-", x = temp.id)

  load(temp.file)
  seurat.obj.list[[temp.id]] <- seurat.obj
}

seurat.obj <- merge(x = seurat.obj.list[[1]], y = seurat.obj.list[2:length(seurat.obj.list)], add.cell.ids = names(seurat.obj.list), min.cells = 0, min.features = 0, do.normalize = F)
seurat.obj[["percent.mt"]] <- PercentageFeatureSet(object = seurat.obj, pattern = "^MT-")
seurat.obj <- NormalizeData(object = seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)
seurat.obj <- FindVariableFeatures(object = seurat.obj, selection.method = "vst")
seurat.obj <- CellCycleScoring(object = seurat.obj, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = T)
seurat.obj <- ScaleData(object = seurat.obj, features = rownames(seurat.obj), vars.to.regress = c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score"))
seurat.obj <- RunPCA(object = seurat.obj, features = VariableFeatures(object = seurat.obj), npcs = 50)
seurat.obj <- FindNeighbors(object = seurat.obj, dims = 1:24, force.recalc = T)
seurat.obj <- FindClusters(object = seurat.obj, resolution = 0.6)
seurat.obj <- RunTSNE(object = seurat.obj, dims = 1:24)
x.markers <- FindAllMarkers(object = seurat.obj.filter, min.pct = 0.25, logfc.threshold = 0.25)


cca.obj.list <- list()
cca.obj.list[["PT"]] <- CreateSeuratObject(counts = "LUAD.umi.matrix", meta.data = "LUAD.meta")
cca.obj.list[["PE"]] <- seurat.obj


for (i in 1:length(cca.obj.list)) {
  cca.obj.list[[i]] <- NormalizeData(cca.obj.list[[i]], verbose = F)
  cca.obj.list[[i]] <- FindVariableFeatures(cca.obj.list[[i]], selection.method = "vst", nfeatures = 2000, verbose = F)
}
anchors <- FindIntegrationAnchors(object.list = cca.obj.list, dims = 1:20)
integrated <- IntegrateData(anchorset = anchors, dims = 1:20)

DefaultAssay(object = integrated) <- "integrated"
integrated <- ScaleData(object = integrated, verbose = F)
integrated <- RunPCA(object = integrated, npcs = 20, verbose = F)
integrated <- FindNeighbors(object = integrated, dims = 1:20, force.recalc = T)
integrated <- FindClusters(object = integrated, resolution = 0.6)
integrated <- RunTSNE(object = integrated, reduction = "pca", dims = 1:20)
