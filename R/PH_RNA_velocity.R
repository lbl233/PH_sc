setwd("~/Mydata/PAH/")

PAH_Hy.list <- readRDS("Hypoxia.rds")
PAH_HySu.list <- readRDS("HySu.rds")

PAH_mus.list <- PAH_Hy.list

for (i in 1:length(PAH_HySu.list)) {
  PAH_mus.list[[i+5]] <- PAH_HySu.list[[i]]
}
for (i in 1:11){
  seur <- PAH_mus.list[[i]]
  seur[['percent.mt']] = PercentageFeatureSet(seur, pattern = '^mt-')
  DefaultAssay(seur) <- "RNA"
  seur <- NormalizeData(seur)
  seur <- FindVariableFeatures(seur)
  PAH_mus.list[[i]] <- seur
}
PAH_mus.anchors <- FindIntegrationAnchors(object.list = PAH_mus.list, dims = 1:20)
PAH_mus <- IntegrateData(anchorset = PAH_mus.anchors, dims = 1:20)


VlnPlot(PAH_mus, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0, group.by = "orig.ident")
PAH_mus <- subset(PAH_mus, percent.mt < 20 & nFeature_RNA < 5000)
DefaultAssay(PAH_mus) <- "RNA"
PAH_mus <- FindVariableFeatures(PAH_mus, selection.method = "vst", nfeatures = 2000)
DefaultAssay(PAH_mus) <- "integrated"
all.genes = rownames(PAH_mus)
PAH_mus = ScaleData(PAH_mus, verbose = TRUE)
PAH_mus <- RunPCA(PAH_mus, npcs = 30, verbose = FALSE)
PAH_mus <- RunTSNE(PAH_mus, reduction = "pca", dims = 1:20) ##### option1 #########
PAH_mus <- RunUMAP(PAH_mus, reduction = "pca", dims = 1:20) ##### option2 #########
PAH_mus <- FindNeighbors(PAH_mus, reduction = "pca", dims = 1:20)
PAH_mus <- FindClusters(PAH_mus, resolution = 0.3)
DimPlot(PAH_mus, label = TRUE)
DimPlot(PAH_mus, label = TRUE, split.by = "Group")
EMT <- subset(PAH_mus_RV, idents = c(1,6,7,15,16,20,21,22))
EMT <- subset(EMT, nCount_spliced > 1000)
DefaultAssay(EMT) <- "RNA"
EMT <- FindVariableFeatures(EMT, selection.method = "vst", nfeatures = 2000)
DefaultAssay(EMT) <- "integrated"
all.genes = rownames(EMT)
EMT = ScaleData(EMT, verbose = FALSE)
EMT <- RunPCA(EMT, npcs = 30, verbose = FALSE)
EMT <- RunTSNE(EMT, reduction = "pca", dims = 1:20) ##### option1 #########
EMT <- RunUMAP(EMT, reduction = "pca", dims = 1:20) ##### option2 #########
EMT <- FindNeighbors(EMT, reduction = "pca", dims = 1:20)
EMT <- FindClusters(EMT, resolution = 0.3)

DimPlot(EMT)
#calculate cell distance based on embedding 
EMT_VF <- subset(EMT, features = VariableFeatures(EMT))

emb <- EMT@reductions$umap@cell.embeddings
emb <- EMT_VF@reductions$umap@cell.embeddings ##umap embedding##
emb <- EMT_VF@reductions$pca@cell.embeddings ##pca embedding##
cell.dist <- as.dist(1-armaCor(t(emb)))
fit.quantile <- 0.02


emat <- EMT@assays$spliced@counts[VariableFeatures(EMT),]
nmat <- EMT@assays$unspliced@counts[VariableFeatures(EMT),]
cluster.label <- EMT_VF$seurat_clusters
emat <- filter.genes.by.cluster.expression(emat,cluster.label,min.max.cluster.average = 0.2)
nmat <- filter.genes.by.cluster.expression(nmat,cluster.label,min.max.cluster.average = 0.05)
length(intersect(rownames(emat),rownames(nmat)))
ncol(EMT_VF)
rvel.cd <- gene.relative.velocity.estimates(emat,nmat,deltaT = 1,
                                            kCells=10,
                                            cell.dist=cell.dist,
                                            fit.quantile=fit.quantile,
                                            n.cores=24)
gg <- PCAPlot(EMT_VF)
gg <- UMAPPlot(EMT)
head(ggplot_build(gg)$data)
colors <- as.list(ggplot_build(gg)$data[[1]]$colour)
names(colors) <- rownames(emb)
show.velocity.on.embedding.cor(emb,rvel.cd,n=30,scale='sqrt',
                                     cell.colors=ac(colors,alpha=0.5),
                                     cex=0.8,arrow.scale=2,show.grid.flow=T,
                                     min.grid.cell.mass=1.0,grid.n=50,arrow.lwd=1,
                                     do.par=F,cell.border.alpha = 0.1,
                                     n.cores=24,main="Endo Epi Velocity")

EMT.list <- SplitObject(EMT, split.by = "Group")
emb1 <- EMT.list[[3]]@reductions$umap@cell.embeddings
emat1 <- EMT.list[[3]]@assays$spliced@counts[VariableFeatures(EMT),]
nmat1 <- EMT.list[[3]]@assays$unspliced@counts[VariableFeatures(EMT),]

cell.dist <- as.dist(1-armaCor(t(emb1)))

rvel.cd <- gene.relative.velocity.estimates(emat1,nmat1,deltaT = 1,
                                            kCells=10,
                                            cell.dist=cell.dist,
                                            fit.quantile=fit.quantile,
                                            n.cores=24)
head(ggplot_build(gg)$data)
gg <- UMAPPlot(EMT.list[[3]])
colors <- as.list(ggplot_build(gg)$data[[1]]$colour)
names(colors) <- rownames(emb1)
show.velocity.on.embedding.cor(emb1,rvel.cd,n=30,scale='sqrt',
                               cell.colors=ac(colors,alpha=0.5),
                               cex=0.8,arrow.scale=2,show.grid.flow=T,
                               min.grid.cell.mass=1.0,grid.n=50,arrow.lwd=1,
                               do.par=F,cell.border.alpha = 0.1,
                               n.cores=24,main="Endo Epi Velocity")
