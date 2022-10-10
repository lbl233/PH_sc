library(Seurat)
library(dplyr)
library(ggplot2)
library(NMF)
library(parallel)

load("/data3/zhouyt/R-image/PAH_sc/5_integrate.RData")
S$annot_new <- as.character(S$annot)
S$annot_new[which(S$annot_new == "SMC_vascular")] <- "SMC"
S$annot_new[which(S$annot_new == "SMC_airway")] <- "SMC"
S$annot_new[which(S$annot_new == "Macrophage_alvelor")] <- "Macrophage_alveolar"
S$annot_new[which(S$annot_new == "SMC/Pericyte")] <- "Pericyte"
S$annot_new[which(S$annot_new == "Endothelial progenitor")] <- "Endothelial capillary"
S$annot_new[which(S$annot_new == "Monocyte")] <- "Monocyte_classical"
S$annot_new[which(S$annot_new == "Fibroblast_alvelor")] <- "Fibroblast_alveolar"

#Find highly expressed genes (HEGs)
Idents(S) <- "annot_new"
DefaultAssay(S) <- "RNA_homolog"
S.marker <- FindAllMarkers(S,only.pos = T)
S.marker = mclapply(Cl, mc.cores = 8, function(x){
  markers <- FindMarkers(S, ident.1 = x, only.pos = T)
  markers$cluster <- x
  return(markers)
})
S_marker_sub_1 <- subset(S.marker, subset = (cluster %in% Cl) & (avg_logFC > 1))
gene.list.1 <- unique(S_marker_sub_1$gene)

length(unique(S$annot_clean_for_corr))
S$annot_clean_for_corr[which(S$annot_clean_for_corr == "SMC_vascular")] <- "SMC"
S$annot_clean_for_corr[which(S$annot_clean_for_corr == "SMC_airway")] <- "SMC"
S$annot_clean_for_corr[which(S$annot_clean_for_corr == "Macrophage_alvelor")] <- "Macrophage_alveolar"
S$annot_clean_for_corr[which(S$annot_clean_for_corr == "SMC/Pericyte")] <- "Pericyte"
S$annot_clean_for_corr[which(S$annot_clean_for_corr == "Endothelial progenitor")] <- "Endothelial capillary"
Cl <- levels(factor(as.vector(S$annot_clean_for_corr)))
Species <- levels(factor(as.vector(S$species)))
DefaultAssay(S) <- "RNA_homolog"
S <- subset(S, annot_clean_for_corr %in% Cl)
S <- NormalizeData(S)
Idents(S) <- "annot_clean_for_corr"

#calculate correlation of each cell type cross species
S.list <- SplitObject(S, split.by = "species")
mouse_seur <- S.list[["mouse"]]
human_seur <- S.list[["human"]]
rat_seur <- S.list[["rat"]]

#NMF reduction and correlation 
#get expression matrix
mouse_seur.list <- SplitObject(mouse_seur, split.by = "annot_clean_for_corr")
mouse_mat.list <- sapply(c(1:length(Cl)), function(i) as.matrix(mouse_seur.list[[Cl[i]]]@assays$RNA_homolog@data))
rm(mouse_seur, mouse_seur.list)

human_seur.list <- SplitObject(human_seur, split.by = "annot_clean_for_corr")
human_mat.list <- sapply(c(1:length(Cl)), function(i) human_seur.list[[Cl[i]]]@assays$RNA_homolog@data)
rm(human_seur, human_seur.list)

rat_seur.list <- SplitObject(rat_seur, split.by = "annot_clean_for_corr")
rat_mat.list <- sapply(c(1:length(Cl)), function(i) rat_seur.list[[Cl[i]]]@assays$RNA_homolog@data)
rm(rat_seur, rat_seur.list)
save(human_mat.list, mouse_mat.list, gene.list.1, Cl, file = 'image/mat_for_NMF.RData')

#run NMF
ptm <- proc.time()
w_human = lapply(human_mat.list, function(x){
  x = as.data.frame(x)
  x = x[gene.list.1,]
  x = x[rowSums(x)!=0,]
  x = basis(nmf(x, 5))
})
print(proc.time() - ptm)

ptm <- proc.time()
w_mouse = lapply(mouse_mat.list, function(x){
  x = as.data.frame(x)
  x = x[gene.list.1,]
  x = x[rowSums(x)!=0,]
  x = basis(nmf(x, 5))
})
print(proc.time() - ptm)

ptm <- proc.time()
w_rat = lapply(rat_mat.list, function(x){
  x = as.data.frame(x)
  x = x[gene.list.1,]
  x = x[rowSums(x)!=0,]
  x = basis(nmf(x, 5))
})
print(proc.time() - ptm)


human_mouse_nmf <- matrix(0, nrow = length(w_human), ncol = length(w_mouse))
for (i in 1:length(w_human)) {
  for (j in 1:length(w_mouse)) {
    gene.overlapped <- intersect(rownames(w_human[[i]]), rownames(w_mouse[[j]]))
    human_mouse_nmf[i,j] <- cor(as.vector(w_human[[i]][gene.overlapped,]), as.vector(w_mouse[[j]][gene.overlapped,]))
  }
}

human_rat_nmf <- matrix(0, nrow = length(w_human), ncol = length(w_rat))
for (i in 1:length(w_human)) {
  for (j in 1:length(w_rat)) {
    gene.overlapped <- intersect(rownames(w_human[[i]]), rownames(w_rat[[j]]))
    human_rat_nmf[i,j] <- cor(as.vector(w_human[[i]][gene.overlapped,]), as.vector(w_rat[[j]][gene.overlapped,]))
  }
}

#calculate correlation
human_mouse_cor <- sapply(1:length(w_human), function(i) human_mouse_nmf[i,i])
human_rat_cor <- sapply(1:length(w_human), function(i) human_rat_nmf[i,i])
PH_corr <- data.frame(mouse_corr = human_mouse_cor, 
                      rat_corr = human_rat_cor,
                      celltype = Cl)
pdf("/data4/lbl/PH/Figure/Fig1F.pdf", width = 10, height = 8)
ggplot(PH_corr, aes(x = mouse_corr, y = rat_corr, colour = main)) + geom_point() + 
  geom_abline(slope = 1, intercept = 0) + 
  geom_text_repel(aes(mouse_corr, rat_corr, label = celltype)) + coord_equal()
dev.off()