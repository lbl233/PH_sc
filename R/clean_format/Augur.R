library(Seurat)
library(Augur)

load("/data3/zhouyt/R-image/PAH_sc/5_integrate.RData")
S$annot_new <- as.character(S[['annot']][,1])
S[["annot_new"]][which(S$annot_new == "SMC_vascular"),] <- "SMC"
S[["annot_new"]][which(S$annot_new == "SMC_airway"),] <- "SMC"
S[["annot_new"]][which(S$annot_new == "Macrophage_alvelor"),] <- "Macrophage_alveolar"
S[["annot_new"]][which(S$annot_new == "SMC/Pericyte"),] <- "Pericyte"

S[["cell_type"]] <- as.factor(S$annot_new)
table(S$cell_type)
S[["label"]] <- as.factor(S$id1)
table(S$id1)
S_PH <- subset(S, species != "human")
rm(S)
table(S_PH$batch)
S_PH.list <- SplitObject(S_PH, split.by = "batch")
rm(S_PH)
augur.list <- list()
for(i in 1:3){
  seur <- S_PH.list[[i]]
  stat <- as.data.frame(table(seur$annot_new))
  Cl <- as.character(stat[stat$Freq > 100,1])
  print(length(Cl))
  seur <- subset(seur, annot_new %in% Cl)
  augur.list[[i]] <- calculate_auc(seur)
}
names(augur.list) <- names(S_PH.list)
saveRDS(augur.list, file = "/data4/lbl/PH/augur.rds")
print("finish")

augur.list[["mouse_1"]] <- augur_hy
augur.list[["mouse_2"]] <- augur_HySu
augur.list[["rat"]] <- augur_mct
for (i in 1:3) {
  seur <- S_PH.list[[i]]
  seur[["augur_auc"]]<- 0
  augur <- augur.list[[i]]
  augur_auc <- as.data.frame(augur$AUC)
  stat <- as.data.frame(table(seur$annot_new))
  Cl <- as.character(stat[stat$Freq > 100,1])
  rownames(augur_auc) <- augur_auc[,1]
  for(j in 1:length(Cl)){
    seur$augur_auc[which(seur$annot_new == Cl[j])] = augur_auc[Cl[j],2]
  }
  seur$augur_auc[is.na(seur$augur_auc)] <- 0
  S_PH.list[[i]] <- seur
}

pdf("/data4/lbl/PH/Figure/Fig2B.pdf")
FeaturePlot(S_PH.list[[1]], features = "augur_auc",) + scale_color_gradientn(colors=c(alpha("#5383c8",1), "white","firebrick")) + coord_equal() + ggtitle("Hx")
FeaturePlot(S_PH.list[[2]], features = "augur_auc") + scale_color_gradientn(colors=c(alpha("#5383c8",1), "white","firebrick")) + coord_equal() + ggtitle("SuHx")
FeaturePlot(S_PH.list[[3]], features = "augur_auc") + scale_color_gradientn(colors=c(alpha("#5383c8",1), "white","firebrick")) + coord_equal() + ggtitle("MCT")
dev.off()