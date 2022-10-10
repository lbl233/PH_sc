library(Seurat)
library(clustree)

#load integrated singel cell data
load("/data3/zhouyt/R-image/PAH_sc/5_integrate.RData")
S$annot_new <- as.character(S$annot)
S$annot_new[which(S$annot_new == "SMC_vascular")] <- "SMC"
S$annot_new[which(S$annot_new == "SMC_airway")] <- "SMC"
S$annot_new[which(S$annot_new == "Macrophage_alvelor")] <- "Macrophage_alveolar"
S$annot_new[which(S$annot_new == "SMC/Pericyte")] <- "Pericyte"
S$annot_new[which(S$annot_new == "Endothelial progenitor")] <- "Endothelial capillary"
S$annot_new[which(S$annot_new == "Monocyte")] <- "Monocyte_classical"
S$annot_new[which(S$annot_new == "Fibroblast_alvelor")] <- "Fibroblast_alveolar"
DefaultAssay(S) <- "RNA_homolog"
S <- NormalizeData(S)
S$annot_new[colnames(PH_macro)] <- paste0("Macrophage", "_", PH_macro$res_used, PH_macro$fun_annot)
S_PH <- subset(S, species != "human")

#get macrophage barcodes
macrophages <- rownames(S_PH[["annot_new"]])[grep("Macrophage", S_PH$annot_new)]
#subset macrophage
PH_macro <- subset(S_PH, cells = macrophages)

#re-run PCA and re-clustring
DefaultAssay(PH_macro) <- "RNA_homolog"
PH_macro = FindVariableFeatures(PH_macro)
DefaultAssay(PH_macro) <- "integrated"
DefaultAssay(PH_macro)
all.genes = rownames(PH_macro)
PH_macro = ScaleData(PH_macro, features = all.genes)
PH_macro = RunPCA(PH_macro,features = VariableFeatures(PH_macro))
PH_macro = FindNeighbors(PH_macro, dims = 1:50)
# S = FindClusters(S, resolution = .1)
PH_macro = RunUMAP(PH_macro, dims = 1:15)
PH_macro = RunTSNE(PH_macro, dims = 1:20)
DimPlot(PH_macro, reduction = "tsne")
DimPlot(PH_macro, reduction = "tsne", split.by = "id1", ncol = 2)
DimPlot(PH_macro, reduction = "tsne", group.by = "id1")
DimPlot(PH_macro, reduction = "umap")
#Fig3A
pdf("/data4/lbl/PH/Figure/Fig3A.pdf", width = 6, height = 5)
DimPlot(PH_macro, reduction = "tsne")
DimPlot(PH_macro, reduction = "tsne",group.by = "id1")
dev.off()
#rename old seurat clusters
res_tree = c(0,.03,.06,.12,.2,.3,.4,.5,.8,1,1.5,2,3,5,8)
for (i in 1:length(res_tree)){
  PH_macro[[paste('initial_cluster_res.', res_tree, sep = '')]] <- PH_macro[[paste('integrated_snn_res.', res_tree, sep = '')]]
}
pdf("/data4/lbl/PH/Figure/Fig3A.pdf", width = 6, height = 5)
DimPlot(PH_macro, reduction = "tsne", label = T)
DimPlot(PH_macro, reduction = "tsne",group.by = "id1")
dev.off()

# PH_macro = FindClusters_re(PH_macro, resolution = .1)
for (i in res_tree) {
  PH_macro = FindClusters(PH_macro, resolution = i)
}
res_tree_data = PH_macro[[paste('integrated_snn_res.', res_tree, sep = '')]]
clustree(res_tree_data, prefix = 'integrated_snn_res.')

Idents(PH_macro) <- "integrated_snn_res.0.3"

#find HEGs of macrophage subtypes
Macro.marker = FindAllMarkers(PH_macro, assay = 'RNA_homolog', only.pos = T)

#functional identification of macrophage subtypes

library(clusterProfiler)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(stringr)
library(enrichplot)

#GSEA
Macro_Go.list <- lapply(c(1:9), function(i){
  df <- subset(Macro.marker, cluster == (i))
  gene <-str_trim(df$gene,"both")
  gene =bitr(gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Mm.eg.db")
  gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)
  gene_df <- data.frame(logFC= df$avg_logFC,
                        SYMBOL = df$gene) 
  gene_df <- merge(gene_df,gene,by="SYMBOL")
  geneList<-gene_df$logFC
  names(geneList)=gene_df$SYMBOL
  geneList=sort(geneList,decreasing = TRUE)
  GO <- gseGO(geneList, ont = "BP", OrgDb = 'org.Mm.eg.db', pvalueCutoff = 1, keyType = "SYMBOL")
})

Macro_Go_0.0.6.list <- lapply(c(1:6), function(i){
  df <- subset(Macro.marker.0.06, cluster == (i-1))
  gene <-str_trim(df$gene,"both")
  gene =bitr(gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Mm.eg.db")
  gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)
  gene_df <- data.frame(logFC= df$avg_logFC,
                        SYMBOL = df$gene) 
  gene_df <- merge(gene_df,gene,by="SYMBOL")
  geneList<-gene_df$logFC
  names(geneList)=gene_df$SYMBOL
  geneList=sort(geneList,decreasing = TRUE)
  GO <- gseGO(geneList, ont = "BP", OrgDb = 'org.Mm.eg.db', pvalueCutoff = 1, keyType = "SYMBOL")
})
GO_fun <- c("GO:0031589", "GO:0042060", "GO:0001525", "GO:0001568", "GO:0033993", 
            "GO:1902105", "GO:0030595", "GO:0032640", "GO:0032635", "GO:0032652",
            "GO:0001816", "GO:0008285", "GO:0010556", "GO:0001944", "GO:0002694", 
            "GO:0000281", "GO:0000278", "GO:0000280", "GO:0051301", "GO:1902850"
            )
Macro_GO_fun_annot.list <- lapply(c(1:4), function(i){Macro_Go_0.0.6.list[[i]]@result[GO_fun,]})
Macro_GO_fun_annot <- rbind(Macro_Go_0.0.6.list[[1]]@result[GO_fun[1:5],],
                            Macro_Go_0.0.6.list[[2]]@result[GO_fun[6:10],],
                            Macro_Go_0.0.6.list[[3]]@result[GO_fun[11:15],],
                            Macro_Go_0.0.6.list[[4]]@result[GO_fun[16:20],])
Macro_GO_fun_annot$Description <- factor(Macro_GO_fun_annot$Description, levels = rev(Macro_GO_fun_annot$Description))
Macro_GO_fun_annot$Cluster <- factor(Macro_GO_fun_annot$Cluster, levels = c("Tissue remodeling", "Pro-inflammation", "Cytokine producing", "Proliferating"))
# Macro_GO_fun_annot <- Reduce(rbind, Macro_GO_fun_annot.list)
Macro_GO_fun_annot$Cluster <- rep(c("Tissue remodeling", "Pro-inflammation", "Cytokine producing", "Proliferating"), each = 5)
pdf("/data4/lbl/PH/Figure/Fig3A2.pdf")
ggplot(Macro_GO_fun_annot, aes(x=NES,y=Description, fill = Cluster)) +
  geom_bar(stat="identity", position="dodge", width=.5) + #theme(axis.text.x=element_text(size=10, angle = 90, vjust=0.3, hjust = 1)) +
  theme(panel.background = element_rect(fill = "white"), panel.grid.major= element_line(),
        axis.line.x = element_line(color="black"),axis.line.x.top = element_line(color="black"), axis.line.y.right = element_line(color="black"),
        axis.line.y = element_line(color="black")) 
dev.off()
#DEGs
S_PH.list <- SplitObject(PH_macro, split.by = "batch")
S_mouse_HySu <- S_PH.list$mouse_2
S_mouse_hy <- S_PH.list$mouse_1
S_rat <- S_PH.list$rat
rm(S_PH.list)
table(PH_macro$res_used, PH_macro$id1)
PH_macro$res_used_id1 <- paste0(PH_macro$res_used, "_", PH_macro$id1)
Idents(PH_macro) <- "res_used_id1"
HySu_DEG.list <- lapply(1:3, function(i) {
  ident1 <- paste0(i, "_", "s")
  ident2 <- paste0(i, "_", "c")
  FindMarkers(S_mouse_HySu, ident.1 = ident1, ident.2 = ident2, only.pos = F, assay = "RNA_homolog", slot = "data")
})

rm(hy_DEG.list)
hy_DEG.list <- lapply(1:3, function(i) {
  ident1 <- paste0(i, "_", "h")
  ident2 <- paste0(i, "_", "c")
  FindMarkers(S_mouse_hy, ident.1 = ident1, ident.2 = ident2, only.pos = F, assay = "RNA_homolog", slot = "data")
})

rm(MCT_DEG.list)
MCT_DEG.list <- lapply(1:3, function(i) {
  ident1 <- paste0(i, "_", "m")
  ident2 <- paste0(i, "_", "c")
  FindMarkers(S_rat, ident.1 = ident1, ident.2 = ident2, only.pos = F, assay = "RNA_homolog", slot = "data")
})

#Signature Genes
PH_macro_sub <- subset(PH_macro, res_used %in% c(1,2,3,4,5,6))
PH_macro_sub$res_used_reorder <- factor(PH_macro_sub$res_used, levels = rev(c(1,3,2,4,6,5)))
Macro_SG <- c("Fabp1", "Cd9", "Ctsd", "Spp1", "Dhcr24", "Sqle", "Il1b", "Ccl2", "Cxcl16", "Trem2",  "Prg4", "Igf1", "Actb", "Mki67", "Top2a")
p <- DotPlot(PH_macro_sub, features = Macro_SG, assay = "RNA_homolog", group.by = "res_used_reorder") + RotatedAxis()
#Fig3C
pdf("/data4/lbl/PH/Figure/Fig3C.pdf", width = 10,height = 5)
p + scale_color_distiller(palette="RdYlBu",type = "div")+
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major= element_line(color="grey", linetype = "dashed"))+
  theme(axis.title.x=element_blank(),  # X axis title
        axis.title.y=element_blank(),  # Y axis title
        # axis.text.x=element_text(size=10, 
        #                          angle = 45,
        #                          vjust=.6),  # X axis text
        axis.text.y=element_text(size=10),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank())+
  scale_size_continuous(range = c(1,10))
dev.off()

# compare cell proportion
PH_macro$id3 <- paste0(PH_macro$batch, "_", PH_macro$id2)
df_stat <- as.data.frame(table(PH_macro$res_used, PH_macro$id3))
group <- c("Control", "Control","Hx", "Hx","Hx","Control","Control","Control","SuHx","SuHx","SuHx","Control","Control","Control","MCT","MCT","MCT")
df_stat$group <- rep(group, each = K)
K <- length(unique(PH_macro$res_used))
N <- length(unique(PH_macro$id3))
end_idx <- 0 
for(i in 1:N){
  start_idx <- end_idx+1
  end_idx <- i*K
  df_stat$Freq[start_idx:end_idx] <- (df_stat$Freq[start_idx:end_idx]/sum(df_stat$Freq[start_idx:end_idx]))
}
group <- c("Control", "Control","Hx", "Hx","Hx","Control","Control","Control","SuHx","SuHx","SuHx","Control","Control","Control","MCT","MCT","MCT")
df_stat$group <- rep(group, each = K)
library(ggpubr)
library(farver)
compare_means(Freq ~ group, data = df_stat_sub)
my_comparisons <- list( c('Control','Hx'), c('Control','SuHx') , c('Control','MCT'))

df_stat_sub <- subset(df_stat, Var1 %in% c(1,3))
df_stat_sub$Var1 <- factor(df_stat_sub$Var1, levels = c(1,3))
p <- ggboxplot(df_stat_sub, x = "group", y = "Freq",
               color = "group", palette = "jco",
               add = "jitter") + facet_wrap(~Var1, nrow =1)
p + stat_compare_means(comparisons = my_comparisons)
#Figure3D
pdf("/data4/lbl/PH/Figure/Figure3D.pdf", width = 8)
p + stat_compare_means(comparisons = my_comparisons)
dev.off()

pdf("/data4/lbl/PH/Figure/Supplement/SupFig4A.pdf")
DimPlot(Macro, reduction = "tsne", label = T) + coord_equal()
dev.off()
Macro$species <- factor(Macro$species, levels = c("human", "rat", "mouse"))
pdf("/data4/lbl/PH/Figure/Supplement/SupFig4B.pdf", width = 10,height = 8)
DimPlot(Macro, reduction = "tsne", split.by = "species", label = T, group.by = "annot_rp3", ncol = 2) + NoLegend()
dev.off()

setwd("/data4/lbl/PH/")
Macro$combine <- paste0("combined", Macro$integrated_snn_res.0.3)
tmp1 <- subset(Macro, species == "human")
tmp2 <- subset(Macro, species != "human")
df1 <- as.data.frame(table(tmp1$combine, tmp1$annot_rp3))
K <- length(unique(tmp1$combine))
N <- length(unique(tmp1$annot_rp3))
end_idx <- 0 
for(i in 1:N){
  start_idx <- end_idx+1
  end_idx <- i*K
  df1$Freq[start_idx:end_idx] <- (df1$Freq[start_idx:end_idx]/sum(df1$Freq[start_idx:end_idx]))
}
df1 <- df1[,c(2,1,3)]
names(df1) <- c("Source", "Target", "Weight")
df2 <- as.data.frame(table(tmp2$combine, tmp2$annot_rp3))
K <- length(unique(tmp2$combine))
N <- length(unique(tmp2$annot_rp3))
end_idx <- 0 
for(i in 1:N){
  start_idx <- end_idx+1
  end_idx <- i*K
  df2$Freq[start_idx:end_idx] <- (df2$Freq[start_idx:end_idx]/sum(df2$Freq[start_idx:end_idx]))
}
names(df2) <- c("Source", "Target", "Weight")
df_rp <- rbind(df1, df2)
# df_rp <- as.data.frame(table(Macro$annot_rp3, Macro$integrated_snn_res.0.3))
df_rp <- subset(df_rp, Weight != 0)
# df_rp$Var2 <- paste0("combined", df_rp$Var2)
names(df_rp) <- c("Source", "Target", "Weight")
df_rp <- df_rp[-grep("undefined", df_rp$Target),]
Sankeylinks <- df_rp
Sankeynodes<-data.frame(name=unique(c(as.character(Sankeylinks$Source),as.character(Sankeylinks$Target))),stringsAsFactors=FALSE)
Sankeynodes$index<-0:(nrow(Sankeynodes) - 1)
Sankeylinks<-merge(Sankeylinks,Sankeynodes,by.x="Source",by.y="name")
Sankeylinks<-merge(Sankeylinks,Sankeynodes,by.x="Target",by.y="name")
Sankeydata<-Sankeylinks[,c(4,5,3)]
names(Sankeydata)<-c("Source","Target","Value")
Sankeyname<-Sankeynodes[,1,drop=FALSE]
#Fig3F
sankeyNetwork(Links=Sankeydata,Nodes=Sankeyname, Source ="Source",
              Target = "Target", Value = "Value", NodeID = "name",
              #colourScale = JS(ColourScale),
             # colourScale = JS('d3.scaleOrdinal().domain(["Fristående","Nästan_öppet","Halvöppet","Slutet"]).range(["#EDF8E9","#BAE4B3","#74C476","#238B45"])'),
              units = "Quads", fontSize = 12, nodeWidth = 30)




library(monocle)
#SupFig1A
pdf("/data4/lbl/PH/Figure/Supplement/SupFig1A.pdf")
DimPlot(PH_macro, label = T) + coord_equal()
dev.off()
DimPlot(PH_macro, label = T, reduction = "tsne")
table(PH_EC$integrated_snn_res.0.3)
table(PH_macro$res_used)
table(EndoMT_pEC$res_used)
EndoMT_pEC <- subset(PH_macro, res_used %in% c(1,3))
EndoMT_pEC <- subset(EndoMT_pEC, species == "mouse")
EndoMT_pEC_matrix<-as(as.matrix(EndoMT_pEC@assays$RNA_homolog@counts), 'sparseMatrix')
feature_ann<-data.frame(gene_id=rownames(EndoMT_pEC_matrix),gene_short_name=rownames(EndoMT_pEC_matrix))
rownames(feature_ann)<-rownames(EndoMT_pEC_matrix)
EndoMT_pEC_fd<-new("AnnotatedDataFrame", data = feature_ann)
sample_ann<-EndoMT_pEC@meta.data
rownames(sample_ann)<-colnames(EndoMT_pEC_matrix)

EndoMT_pEC_pd<-new("AnnotatedDataFrame", data =sample_ann)

EndoMT_pEC.cds<-newCellDataSet(EndoMT_pEC_matrix,phenoData =EndoMT_pEC_pd,featureData =EndoMT_pEC_fd,expressionFamily=negbinomial.size())
EndoMT_pEC.cds <- estimateSizeFactors(EndoMT_pEC.cds)
EndoMT_pEC.cds <- estimateDispersions(EndoMT_pEC.cds)
EndoMT_pEC.cds <- detectGenes(EndoMT_pEC.cds, min_expr = 0.1)

expressed_genes <- row.names(subset(fData(EndoMT_pEC.cds), num_cells_expressed >= 10))
diff_test_res <- differentialGeneTest(EndoMT_pEC.cds,fullModelFormulaStr = "~res_used") 

ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
ordering_genes <- ordering_genes[-grep("mt-", ordering_genes)]
EndoMT_pEC.cds <- setOrderingFilter(EndoMT_pEC.cds,ordering_genes)

EndoMT_pEC.cds <- reduceDimension(EndoMT_pEC.cds, max_components = 2, reduction_method = "DDRTree")
EndoMT_pEC.cds <- orderCells(EndoMT_pEC.cds)
EndoMT_pEC.cds$Cluster <- paste0("C", EndoMT_pEC.cds$res_used)
#SupFig1B
pdf("/data4/lbl/PH/Figure/Supplement/SupFig1B.pdf")
plot_cell_trajectory(EndoMT_pEC.cds, color_by = "Cluster")
plot_cell_trajectory(EndoMT_pEC.cds, color_by = "Pseudotime")
dev.off()
#Branch expression analysis
BEAM_res <- BEAM(EndoMT_pEC.cds, branch_point = 1, cores = 4)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res_sub <- subset(BEAM_res, pval < 0.05)
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
gl_lasso <- rownames(BEAM_res_sub)
# gl_lasso <- rownames(EndoMT_pEC.cds)
barcode <- colnames(EndoMT_pEC.cds)
# tmp <- colnames(subset(EndoMT_pEC, res_used %in% c(1,3)))
# EndoMT_pEC.cds$Pseudotime_reverse <- max(EndoMT_pEC.cds$Pseudotime) - EndoMT_pEC.cds$Pseudotime
lasso_y <- as.matrix(EndoMT_pEC.cds$Pseudotime)
rownames(lasso_y) <- colnames(EndoMT_pEC.cds)
# barcode <- intersect(tmp, barcode)
lasso_y <- lasso_y[barcode,]
max(lasso_y)
min(lasso_y)
# lasso_y <- max(lasso_y) - lasso_y
lasso_x <- expm1(EndoMT_pEC@assays$RNA_homolog@data[gl_lasso,barcode])
lasso_x <- as.matrix(t(lasso_x))
lasso_rst <- cv.glmnet(x=lasso_x,y=lasso_y,family="gaussian",alpha=1,nfolds = 10)
pdf("/data4/lbl/PH/Figure/Supplement/SupFig1C.pdf")
plot(lasso_rst)
dev.off()

lasso_rst1 <- glmnet(x=lasso_x,y=lasso_y,family="gaussian",alpha=1,lambda = lasso_rst$lambda.min)
fit<-glmnet(x=lasso_x,y=lasso_y,alpha=1,family="gaussian")
plot(fit)

tmp <- coef(lasso_rst1)
tmp_sub <- tmp[which(tmp[,1] != 0),]
plot(lasso_rst1)
tmp_sub <- tmp_sub[order(tmp_sub,decreasing = TRUE)]

lasso_rst3 <- glmnet(x=lasso_x,y=lasso_y,family="gaussian",alpha=1,lambda = lasso_rst$lambda.1se)
tmp3 <- coef(lasso_rst3)
Active.coefficients<-tmp3[which(tmp3[,1] != 0),]
new_tmp <- Active.coefficients
new_tmp <- new_tmp[order(new_tmp,decreasing = TRUE)]

df <- new_tmp
gene <-str_trim(names(df),"both")
gene =bitr(gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Mm.eg.db")
gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)
gene_df <- data.frame(coefficient= df,
                      coefficient_abs = abs(df),
                      SYMBOL = names(df)) 
gene_df <- merge(gene_df,gene,by="SYMBOL")
geneList<-gene_df$coefficient
names(geneList)=gene_df$SYMBOL
geneList=sort(geneList,decreasing = TRUE)
GO <- gseGO(geneList, ont = "BP", OrgDb = 'org.Mm.eg.db', pvalueCutoff = 1, keyType = "SYMBOL")
Macro3_barcode <- colnames(subset(PH_macro, res_used %in% c(1,3)))
valid_x <- expm1(PH_macro@assays$RNA_homolog@data[names(new_tmp)[-1],Macro3_barcode])
valid_x <- as.matrix(valid_x)
Pseudotime_predict <- new_tmp[-1]%*%valid_x + 10.3201
Macro3_pseudo <-Pseudotime_predict 
Macro3_pseudo_df <- as.data.frame(t(Macro3_pseudo))
Macro3_pseudo_df$id1 <- as.character(PH_macro$id1[Macro3_barcode])
Macro3_pseudo_df$res_used_id1 <- as.character(PH_macro$res_used_id1[Macro3_barcode])
Macro3_pseudo_df$res_used <- as.character(PH_macro$res_used[Macro3_barcode])
#SupFig1E
pdf("/data4/lbl/PH/Figure/Supplement/SupFig1E.pdf")
ggplot(Macro3_pseudo_df, aes(x = Macro3_pseudo, y = res_used_id1)) + 
  geom_density_ridges_gradient(aes(fill = id1), scale = 3, size = 0.3) +theme(panel.grid.major =element_blank(), 
                                                                              panel.grid.minor = element_blank(),
                                                                              panel.background = element_blank(),#background
                                                                              panel.border = element_blank())#border
dev.off()

GO_1se_rst <- GO@result
View(GO_1se_rst)
GO_1se_rst <- GO_1se_rst[order(GO_1se_rst$NES),]
GO_1se_rst$Description1 <- factor(GO_1se_rst$Description, levels = GO_1se_rst$Description)
GO_labeled <- GO_1se_rst[c("GO:0032635", "GO:0071356", "GO:0030258", "GO:0071345", "GO:0007179"),]
#SupFig1D2
pdf("/data4/lbl/PH/Figure/Supplement/SupFig1D2.pdf",width = 10)
ggplot(GO_1se_rst, aes(x=Description,y=NES)) +
  geom_bar(stat="identity", position="dodge", width=.5) + #theme(axis.text.x=element_text(size=10, angle = 90, vjust=0.3, hjust = 1)) +
  theme(panel.background = element_rect(fill = "white"), panel.grid.major= element_line(), axis.ticks.x = element_blank(),axis.line.x.bottom = element_line(color="black"),
        axis.line.y.left = element_line(color="black"), axis.line.x.top= element_line(color="black"),
        axis.line.y.right = element_line(color="black"), axis.text.x = element_blank()) + geom_label_repel(data = GO_labeled, aes(label = Description,x=Description,y=NES),box.padding = unit(2, "lines"),fill = "yellow",color = "red",direction = "y")
dev.off()

#SupFig1D1
pdf("/data4/lbl/PH/Figure/Supplement/SupFig1D1.pdf")
plot_genes_branched_pseudotime(EndoMT_pEC.cds[intersect(gl_mouse,names(new_tmp[-1])),],
                               branch_point = 2,
                               color_by = "Cluster",
                               ncol = 1)
plot_genes_branched_pseudotime(EndoMT_pEC.cds[names(new_tmp[2:6]),],
                               branch_point = 1,
                               color_by = "Cluster",
                               ncol = 1)
plot_genes_branched_pseudotime(EndoMT_pEC.cds[names(new_tmp[419:423]),],
                               branch_point = 1,
                               color_by = "Cluster",
                               ncol = 1)
dev.off()
#SupFig1D3
pdf("/data4/lbl/PH/Figure/Supplement/SupFig1D3.pdf",width = 8)
plot_genes_branched_heatmap(EndoMT_pEC.cds[names(new_tmp[c(2:31,394:423)]),], num_clusters = 2,show_rownames = TRUE)
dev.off()