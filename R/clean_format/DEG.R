######################
####calculate DEGs####
######################
load("/data3/zhouyt/R-image/PAH_sc/5_integrate.RData")
S$annot_new <- as.character(S[['annot']][,1])
S[["annot_new"]][which(S$annot_new == "SMC_vascular"),] <- "SMC"
S[["annot_new"]][which(S$annot_new == "SMC_airway"),] <- "SMC"
S[["annot_new"]][which(S$annot_new == "Macrophage_alvelor"),] <- "Macrophage_alveolar"
S[["annot_new"]][which(S$annot_new == "SMC/Pericyte"),] <- "Pericyte"
S[["annot_new"]][which(S$annot_new == "Endothelial progenitor"),] <- "Endothelial capillary"
DefaultAssay(S) <- "RNA_homolog"
S <- NormalizeData(S)
table(S$id1)

S_mouse_HySu <- S_PH.list$mouse_2
S_mouse_hy <- S_PH.list$mouse_1
S_rat <- S_PH.list$rat
# table(S$annot_new)
# S_mouse_HySu
table(S_mouse_HySu$id1)
S_mouse_HySu$annot_id1 <- paste0(S_mouse_HySu$annot_new, "_", S_mouse_HySu$id1)
Idents(S_mouse_HySu) <- "annot_id1"
S_mouse_stat <- as.data.frame(table(S_mouse_HySu$annot_new, S_mouse_HySu$id1))
tmp1 <- subset(S_mouse_stat, Var2 == "c")
tmp2 <- subset(S_mouse_stat, Var2 == "s")
Cl1 <- intersect(as.character(tmp1[tmp1$Freq > 20,1]), #remove the rare cell types while comparing the similarity
                 as.character(tmp1[tmp2$Freq > 20,1]))

table(S_mouse_hy$id1)
S_mouse_hy$annot_id1 <- paste0(S_mouse_hy$annot_new, "_", S_mouse_hy$id1)
Idents(S_mouse_hy) <- "annot_id1"
S_mouse_stat <- as.data.frame(table(S_mouse_hy$annot_new, S_mouse_hy$id1))
tmp1 <- subset(S_mouse_stat, Var2 == "c")
tmp2 <- subset(S_mouse_stat, Var2 == "h")
Cl2 <- intersect(as.character(tmp1[tmp1$Freq > 20,1]),
                 as.character(tmp1[tmp2$Freq > 20,1]))

table(S_rat$id1)
S_rat$annot_id1 <- paste0(S_rat$annot_new, "_", S_rat$id1)
Idents(S_rat) <- "annot_id1"
S_rat_stat <- as.data.frame(table(S_rat$annot_new, S_rat$id1))
tmp1 <- subset(S_rat_stat, Var2 == "c")
tmp2 <- subset(S_rat_stat, Var2 == "m")
Cl3 <- intersect(as.character(tmp1[tmp1$Freq > 20,1]),
                 as.character(tmp1[tmp2$Freq > 20,1]))

rm(HySu_DEG.list)
HySu_DEG.list <- lapply(1:length(Cl1), function(i) {
  ident1 <- paste0(Cl1[i], "_", "s")
  ident2 <- paste0(Cl1[i], "_", "c")
  FindMarkers(S_mouse_HySu, ident.1 = ident1, ident.2 = ident2, only.pos = F, assay = "RNA_homolog", slot = "data")
})
names(HySu_DEG.list) <- Cl1

rm(hy_DEG.list)
hy_DEG.list <- lapply(1:length(Cl2), function(i) {
  ident1 <- paste0(Cl2[i], "_", "h")
  ident2 <- paste0(Cl2[i], "_", "c")
  FindMarkers(S_mouse_hy, ident.1 = ident1, ident.2 = ident2, only.pos = F, assay = "RNA_homolog", slot = "data")
})
names(hy_DEG.list) <- Cl2

rm(MCT_DEG.list)
MCT_DEG.list <- lapply(1:length(Cl3), function(i) {
  ident1 <- paste0(Cl3[i], "_", "m")
  ident2 <- paste0(Cl3[i], "_", "c")
  FindMarkers(S_rat, ident.1 = ident1, ident.2 = ident2, only.pos = F, assay = "RNA_homolog", slot = "data")
})
names(MCT_DEG.list) <- Cl3

Cl <- Reduce(union, list(names(hy_DEG.list), names(HySu_DEG.list), names(MCT_DEG.list)))

HySu_DEG_1.list <- lapply(1:32, function(i) subset(HySu_DEG.list[[i]], subset = (p_val_adj < 0.05)))
names(HySu_DEG_1.list) <- names(HySu_DEG.list)
hy_DEG_1.list <- lapply(1:37, function(i) subset(hy_DEG.list[[i]], subset = (p_val_adj < 0.05)))
names(hy_DEG_1.list) <- names(hy_DEG.list)
MCT_DEG_1.list <- lapply(1:32, function(i) subset(MCT_DEG.list[[i]], subset = (p_val_adj < 0.05)))
names(MCT_DEG_1.list) <- names(MCT_DEG.list)
uDEG.list <- list()
uDEG.list <- sapply(1:42, function(i) length(unique(c(Hypoxia=rownames(hy_DEG_1.list[[Cl[i]]]), HySu=rownames(HySu_DEG_1.list[[Cl[i]]]), MCT=rownames(MCT_DEG_1.list[[Cl[i]]])))))
main <- read.table("/data4/lbl/PH/main.txt", header = F)
DEG_total_genes <- lapply(1:42, function(i) unique(c(Hypoxia=rownames(hy_DEG_1.list[[Cl[i]]]), HySu=rownames(HySu_DEG_1.list[[Cl[i]]]), MCT=rownames(MCT_DEG_1.list[[Cl[i]]]))))
uEDG_cEC <- DEG_total_genes[[10]]
uEDG_AM <- DEG_total_genes[[23]]
GO_uEDG_AM <- enrichGO(gene = uEDG_AM, OrgDb= org.Mm.eg.db, ont = "BP",pAdjustMethod = "BH", keyType = "SYMBOL",pvalueCutoff = 0.05,qvalueCutoff = 0.05)
GO_uEDG_AM_rst <- GO_uEDG_AM@result
GO_uEDG_AM_rst$celltype <- "Macrophage_alveolar"

GO_uEDG_cEC <- enrichGO(gene = uEDG_cEC, OrgDb= org.Mm.eg.db, ont = "BP",pAdjustMethod = "BH", keyType = "SYMBOL",pvalueCutoff = 0.05,qvalueCutoff = 0.05)
GO_uEDG_cEC_rst <- GO_uEDG_cEC@result
GO_uEDG_cEC_rst$celltype <- "Endothelial capillary"

tmp <- rbind(GO_uEDG_cEC_rst[1:10,],GO_uEDG_AM_rst[1:10,])
tmp$Description <- factor(tmp$Description, levels = tmp$Description)

DEG_total_num <- data.frame(celltype = Cl,
                            total = uDEG.list)
DEG_total_num$main <- main$V1
DEG_stat <- DEG_total_num[order(DEG_total_num$main),]
DEG_stat$celltype <- as.character(DEG_stat$celltype)
DEG_stat$celltype[2] <- "Endothelial capillary"
DEG_stat$celltype[39] <- "Fibroblast_alveolar"
DEG_stat$celltype <- factor(DEG_stat$celltype, levels = DEG_stat$celltype)

S_PH$annot_new[which(S_PH$annot_new == "SMC_vascular")] <- "SMC"
S_PH$annot_new[which(S_PH$annot_new == "SMC_airway")] <- "SMC"
S_PH$annot_new[which(S_PH$annot_new == "Macrophage_alvelor")] <- "Macrophage_alveolar"
S_PH$annot_new[which(S_PH$annot_new == "SMC/Pericyte")] <- "Pericyte"
S_PH$annot_new[which(S_PH$annot_new == "Endothelial progenitor")] <- "Endothelial capillary"
S_PH$annot_new[which(S_PH$annot_new == "Monocyte")] <- "Monocyte_classical"
S_PH$annot_new[which(S_PH$annot_new == "Fibroblast_alvelor")] <- "Fibroblast_alveolar"
S_PH$annot_new <- as.factor(S_PH$annot_new)
Idents(S_PH) <-"annot_new"
pdf("/data4/lbl/PH/Figure/Fig2A.pdf", width = 15)
DimPlot(S_PH, split.by = "id1", ncol = 2)
dev.off()


pdf("/data4/lbl/PH/Figure/Fig2C1.pdf", width = 5)
ggplot(DEG_stat, aes(celltype, total)) +
  geom_bar(aes(fill = main), stat="identity", position="dodge", width=.5) + theme(axis.text.x=element_text(size=10, angle = 90, vjust=0.3, hjust = 1)) +
  theme(panel.background = element_rect(fill = "white"), panel.grid.major= element_line(),
        axis.line.x = element_line(color="black"),axis.line.x.top = element_line(color="black"), axis.line.y.right = element_line(color="black"),
        axis.line.y = element_line(color="black"))  + coord_flip()
dev.off()
pdf("/data4/lbl/PH/Figure/Fig2C2.pdf", width = 10)
ggplot(tmp, aes(x=-log10(p.adjust),y=rev(Description), fill = celltype)) +
  geom_bar(stat="identity", position="dodge", width=.5) + theme(axis.text.x=element_text(size=10, angle = 90, vjust=0.3, hjust = 1)) +
  theme(panel.background = element_rect(fill = "white"), panel.grid.major= element_line(),
        axis.line.x = element_line(color="black"),axis.line.x.top = element_line(color="black"), axis.line.y.right = element_line(color="black"),
        axis.line.y = element_line(color="black")) 
dev.off()