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

#get EC barcodes
EC <- rownames(S_PH[["annot"]])[grep("Endothelial", S_PH$annot)]
EC <- rownames(S[["annot"]])[grep("Endothelial", S$annot)]
#subset 
load("/ssd/zhouyt/R/PAH_sc/image/4_integrate.RData")
PH_EC <- subset(S, cells = EC)

#re-run PCA and re-clustring
DefaultAssay(PH_EC) <- "RNA_homolog"
PH_EC = FindVariableFeatures(PH_EC)
DefaultAssay(PH_EC) <- "integrated"
DefaultAssay(PH_EC)
all.genes = rownames(PH_EC)
PH_EC = ScaleData(PH_EC, features = all.genes)
PH_EC = RunPCA(PH_EC,features = VariableFeatures(PH_EC))
PH_EC = FindNeighbors(PH_EC, dims = 1:50)
PH_EC = RunUMAP(PH_EC, dims = 1:50)
PH_EC = RunTSNE(PH_EC, dims = 1:50)
DimPlot(PH_EC, reduction = "tsne")
DimPlot(PH_EC, reduction = "tsne", split.by = "id1", ncol = 2, label = T)
DimPlot(PH_EC, reduction = "tsne", group.by = "id1")
DimPlot(PH_EC, reduction = "umap")
DimPlot(PH_EC, reduction = "umap", split.by = "id1", ncol = 2)
#rename old seurat clusters
res_tree = c(0,.03,.06,.12,.2,.3,.4,.5,.8,1,1.5,2,3,5,8)
for (i in 1:length(res_tree)){
  S_EC[[paste('initial_cluster_res.', res_tree, sep = '')]] <- S_EC[[paste('integrated_snn_res.', res_tree, sep = '')]]
}
S_EC[["res_used"]] <- "human"
S_EC[["res_used"]][Cells(PH_EC),] <- PH_EC[["res_used"]][Cells(PH_EC),]
DimPlot(S_EC, group.by = "res_used", reduction = "tsne")

for (i in res_tree) {
  S_EC = FindClusters(S_EC, resolution = i)
}
res_tree_data = S_EC[[paste('integrated_snn_res.', res_tree, sep = '')]]
clustree(res_tree_data, prefix = 'integrated_snn_res.')

Idents(S_EC) <- "integrated_snn_res.0.5"
#table(PH_EC[["integrated_snn_res.0.3"]])
DimPlot(PH_EC, label = T, reduction = "tsne")
DimPlot(PH_EC, label = T, reduction = "tsne", split.by = "id1", ncol = 2)
DimPlot(PH_EC, label = T, reduction = "umap", split.by = "id1", ncol = 2)

pdf("/data4/lbl/PH/Figure/Fig5A.pdf", width = 6, height = 5)
DimPlot(PH_EC, reduction = "tsne", group.by = "res_used", label = T)
DimPlot(PH_EC, reduction = "tsne",group.by = "id1")
dev.off()


EC_cellcount <- as.data.frame(table(PH_EC$integrated_snn_res.1, PH_EC$id1))
EC_stat <- data.frame(control = EC_cellcount[1:17, 3], hypoxia = EC_cellcount[18:34, 3], mct = EC_cellcount[35:51, 3], HySu = EC_cellcount[52:68, 3])
EC_cellcount <- apply(EC_stat, 2, function(x)round((x/sum(x))*100, digits = 2))

#
EC_cellcount <- apply(EC_stat, 2, function(x)round((x/sum(x))*100, digits = 2))
PH_EC[["res_used"]] <- as.numeric(PH_EC$integrated_snn_res.1)
PH_EC[["res_used"]] <- as.factor(PH_EC$res_used)
Idents(PH_EC) <- "res_used"
table(Idents(PH_EC))
table(PH_EC$res_used)

PH_EC$annot_for_fig <- as.character(PH_EC$res_used)
PH_EC$annot_for_fig[which(PH_EC$res_used %in% c(7,8,9,15,17))] <- PH_EC$fun_annot[which(PH_EC$res_used %in% c(7,8,9,15,17))]

library(future)
plan("multisession", workers = 8)
EC.marker = FindAllMarkers(PH_EC, assay = 'RNA_homolog', only.pos = T)
plan("sequential")
save(list = c("PH_EC", "EC.marker"), file = '/data4/lbl/PH/EC.RData')
DimPlot(PH_EC, label = T, reduction = "tsne")
Idents(PH_EC) <- "integrated_snn_res.0.3"
EC.marker.3 = FindAllMarkers(PH_EC, assay = 'RNA_homolog', only.pos = T)



#functional identification of EC subtypes

library(clusterProfiler)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(stringr)
library(enrichplot)

#GSEA
EC_Go.list <- lapply(c(1:17), function(i){
  df <- subset(EC.marker, cluster == i)
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

#functional annotation
PH_EC$fun_annot <- "undefined"
PH_EC$fun_annot[which(PH_EC$res_used == 11)] <- "angiogenesis cap EC"
PH_EC$fun_annot[which(PH_EC$res_used == 2|PH_EC$res_used == 4|PH_EC$res_used == 13)] <- "apoptotic EC"
PH_EC$fun_annot[which(PH_EC$res_used == 1|PH_EC$res_used == 3|PH_EC$res_used == 14|PH_EC$res_used == 10)] <- "Normal EC"
PH_EC$fun_annot[which(PH_EC$res_used == 17)] <- "Nox2+ EC"
PH_EC$fun_annot[which(PH_EC$res_used == 15)] <- "Lymphatic EC"
PH_EC$fun_annot[which(PH_EC$res_used == 5|PH_EC$res_used == 6)] <- "Aerocyte EC"
PH_EC$fun_annot[which(PH_EC$res_used == 8)] <- "Pericyte-derived EC"
PH_EC$fun_annot[which(PH_EC$res_used == 9)] <- "Artery EC"
PH_EC$fun_annot[which(PH_EC$res_used == 7)] <- "Vein EC: placenta development"
PH_EC$fun_annot[which(PH_EC$res_used == 16)] <- "lipid metabolic EC"
PH_EC$fun_annot[which(PH_EC$res_used == 12)] <- "circulatory EC"
DimPlot(PH_EC, group.by = "fun_annot", label = T, reduction = "tsne")
table(PH_EC$fun_annot)

#GSEA
Idents(PH_EC) <- "fun_annot"
EC.marker.new <- FindAllMarkers(PH_EC, only.pos = T)
Cl <- levels(EC.marker.new$cluster)
EC_Go_new.list <- lapply(c(1:11), function(i){
  df <- subset(EC.marker.new, cluster == Cl[i])
  gene <-str_trim(df$gene,"both")
  gene =bitr(gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Mm.eg.db")
  gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)
  gene_df <- data.frame(logFC= df$avg_logFC,
                        SYMBOL = df$gene) 
  gene_df <- merge(gene_df,gene,by="SYMBOL")
  geneList<-gene_df$logFC
  names(geneList)=gene_df$SYMBOL
  geneList=sort(geneList,decreasing = T)
  GO <- gseGO(geneList, ont = "BP", OrgDb = 'org.Mm.eg.db', pvalueCutoff = 1, keyType = "SYMBOL")
})

PH_EC$annot_new[which(PH_EC$annot_new == "Endothelial progenitor")] <- "Endothelial capillary"
PH_EC$annot_new <- as.character(PH_EC$annot)
PH_EC$fun_annot[which(PH_EC$fun_annot == "circulatory EC")] <- "Circ vEC"
PH_EC$fun_annot[which(PH_EC$fun_annot == "Aerocyte EC")] <- "Circ aEC"
PH_EC$fun_annot[which(PH_EC$fun_annot == "lipid metabolic EC")] <- "Apoptotic EC"
PH_EC$fun_annot[which(PH_EC$fun_annot == "apoptotic EC")] <- "Apoptotic EC"
PH_EC$fun_annot[which(PH_EC$fun_annot == "angiogenesis cap EC")] <- "Angiogenesis cEC"
Cl <- levels(as.factor(PH_EC$fun_annot))
Cl_cEC <- Cl[c(1,2,4,5,7,8)]
Idents(PH_EC) <- "fun_annot"
PH_EC_sub <- subset(PH_EC, fun_annot %in% Cl_cEC)
cEC_HEG.list <- list()
for(i in 1:length(Cl_cEC)){
  ident1 = Cl_cEC[i]
  tmp <- FindMarkers(PH_EC, ident.1 = ident1, only.pos = T)
  tmp$gene <- rownames(tmp)
  cEC_HEG.list[[i]] <- tmp
}
names(cEC_HEG.list) <- Cl_cEC
EC_fun_annot.list <- lapply(c(1:6), function(i){
  df <- cEC_HEG.list[[i]]
  gene <-str_trim(df$gene,"both")
  gene =bitr(gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Mm.eg.db")
  gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)
  gene_df <- data.frame(logFC= df$avg_logFC,
                        SYMBOL = df$gene) 
  gene_df <- merge(gene_df,gene,by="SYMBOL")
  geneList<-gene_df$logFC
  names(geneList)=gene_df$SYMBOL
  geneList=sort(geneList,decreasing = T)
  GO <- gseGO(geneList, ont = "BP", OrgDb = 'org.Mm.eg.db',pvalueCutoff = 1, keyType = "SYMBOL")
})
names(EC_fun_annot.list) <- Cl_cEC
for(i in 1:6){
  pdf(paste0("/data4/lbl/PH/cEC_", names(cEC_HEG.list)[i], "_HEG_GO.pdf"), width = 10)
  p <- plotGOgraph(EC_fun_annot.list[[i]])
  print(p)
  dev.off()
}

GO_fun <- c("GO:0042127", "GO:0043066", "GO:0043069", "GO:0001525", "GO:0072593", # angiogenesis EC 
            "GO:2001233", "GO:0060326", "GO:0010647", "GO:0030336", "GO:0009967", # apoptotic
            "GO:0010811", "GO:0008015", "GO:0043534", "GO:0090132", "GO:0001570", # Circ aEC
            "GO:0042692", "GO:0003013", "GO:0008015", "GO:0045785", "GO:0043066", # Circ vEC
            "GO:0030316", "GO:0072593", "GO:0070555", "GO:1902105", "GO:0007154"  # Nox2+ EC
            
)

EC_GO_fun_annot <- rbind(EC_fun_annot.list[[1]]@result[GO_fun[1:5],],
                            EC_fun_annot.list[[2]]@result[GO_fun[6:10],],
                            EC_fun_annot.list[[3]]@result[GO_fun[11:15],],
                            EC_fun_annot.list[[4]]@result[GO_fun[16:20],],
                            EC_fun_annot.list[[6]]@result[GO_fun[21:25],])
EC_GO_fun_annot$Description <- factor(EC_GO_fun_annot$Description, levels = rev(EC_GO_fun_annot$Description))
EC_GO_fun_annot$Cluster <- rep(names(EC_fun_annot.list)[-5], each = 5)
EC_GO_fun_annot$Cluster <- factor(EC_GO_fun_annot$Cluster, levels = names(EC_fun_annot.list)[-5])
#Fig5A2
pdf("/data4/lbl/PH/Figure/Fig5A2.pdf")
ggplot(EC_GO_fun_annot, aes(x=NES,y=Description, fill = Cluster)) +
  geom_bar(stat="identity", position="dodge", width=.5) + #theme(axis.text.x=element_text(size=10, angle = 90, vjust=0.3, hjust = 1)) +
  theme(panel.background = element_rect(fill = "white"), panel.grid.major= element_line(),
        axis.line.x = element_line(color="black"),axis.line.x.top = element_line(color="black"), axis.line.y.right = element_line(color="black"),
        axis.line.y = element_line(color="black")) 
dev.off()

EC.genelist <- c("Tek", "Cldn5", "Vcam1", "Vwf", "Plvap", "Ntrk2", "Cxcl12", "Ednrb", "Kdr", "Hopx", "Car4", "Npr3","Ackr3", "Cybb", "Il1b")

PH_EC_sub$fun_annot_reorder <- factor(PH_EC_sub$fun_annot, levels = rev(c("Normal EC","Apoptotic EC", "Circ aEC", "Circ vEC", "Angiogenesis cEC", "Nox2+ EC")))
p <- DotPlot(PH_EC_sub, features = EC.genelist, assay = "RNA_homolog", group.by = "fun_annot_reorder") + RotatedAxis()
#Fig5C
pdf("/data4/lbl/PH/Figure/Fig5C.pdf", width = 10,height = 5)
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

#calculate cell proportion
PH_EC$id3 <- paste0(PH_EC$batch, "_", PH_EC$id2)
df_stat <- as.data.frame(table(PH_EC$fun_annot, PH_EC$id3))
group <- c("Control", "Control","Hx", "Hx","Hx","Control","Control","Control","SuHx","SuHx","SuHx","Control","Control","Control","MCT","MCT","MCT")

K <- length(unique(PH_EC$fun_annot))
N <- length(unique(PH_EC$id3))
df_stat$group <- rep(group, each = K)
end_idx <- 0 
for(i in 1:N){
  start_idx <- end_idx+1
  end_idx <- i*K
  df_stat$Freq[start_idx:end_idx] <- (df_stat$Freq[start_idx:end_idx]/sum(df_stat$Freq[start_idx:end_idx]))
}
library(ggpubr)
library(farver)
compare_means(Freq ~ group, data = df_stat)
my_comparisons <- list( c('Control','Hx'), c('Control','SuHx') , c('Control','MCT'))

df_stat_sub <- subset(df_stat, Var1 %in% c("Normal EC","Apoptotic EC", "Circ aEC", "Circ vEC", "Angiogenesis cEC", "Nox2+ EC"))
df_stat_sub$Var1 <- factor(df_stat_sub$Var1, levels = c("Normal EC","Apoptotic EC", "Circ aEC", "Circ vEC", "Angiogenesis cEC", "Nox2+ EC"))
p <- ggboxplot(df_stat_sub, x = "group", y = "Freq",
               color = "group", palette = "jco",
               add = "jitter") + facet_wrap(~Var1, nrow =2)
p + stat_compare_means(comparisons = my_comparisons)
pdf("/data4/lbl/PH/Figure/Figure5D.pdf", width = 10, height = 10)
p + stat_compare_means(comparisons = my_comparisons)
dev.off()

pdf("/data4/lbl/PH/Figure/Supplement/SupFig2.pdf", width = 8, height = 10)
FeaturePlot(PH_EC, features = EC.genelist, ncol = 3, cols=c("darkblue", "goldenrod","firebrick"), reduction = "tsne") 
dev.off()
PH_EC$fun_annot[which(PH_EC$fun_annot == "Vein EC: placenta development")] <- "Vein EC"

#Fig5F
pdf("/data4/lbl/PH/Figure/Fig5F.pdf", width = 8, height = 6)
p <- DotPlot2(PH_EC, features = "Cortellis_EC1", group.1 = "id1", group.2 = "fun_annot")
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

#####################################################
##### river plot connect human and animal model #####
#####################################################
library(riverplot)
library(networkD3)
load("/data3/zhouyt/R-image/PAH_sc/5_integrate.RData")
#get macrophage barcodes
S$annot_new <- as.character(S$annot)
S$annot_new[which(S$annot_new == "SMC_vascular")] <- "SMC"
S$annot_new[which(S$annot_new == "SMC_airway")] <- "SMC"
S$annot_new[which(S$annot_new == "Macrophage_alvelor")] <- "Macrophage_alveolar"
S$annot_new[which(S$annot_new == "SMC/Pericyte")] <- "Pericyte"
S$annot_new[which(S$annot_new == "Endothelial progenitor")] <- "Endothelial capillary"
S$annot_new[which(S$annot_new == "Monocyte")] <- "Monocyte_classical"
S$annot_new[which(S$annot_new == "Fibroblast_alvelor")] <- "Fibroblast_alveolar"
Endothelial <- rownames(S[["annot_new"]])[grep("Endothelial", S$annot_new)]
#subset macrophage
EC<- subset(S, cells = Endothelial)
rm(S)
DefaultAssay(EC) <- "RNA_homolog"
EC = FindVariableFeatures(EC)
DefaultAssay(EC) <- "integrated"
DefaultAssay(EC)
all.genes = rownames(EC)
EC = ScaleData(EC, features = all.genes)
EC = RunPCA(EC,features = VariableFeatures(PH_EC))
EC = FindNeighbors(EC, dims = 1:50)
EC = RunTSNE(EC, dims = 1:30)
EC$annot_rp <- EC$annot_new
EC$annot_rp2 <- EC$annot_new
EC$annot_rp[colnames(PH_EC)] <- PH_EC$fun_annot
EC$annot_rp[colnames(PH_macro)] <- PH_EC$res_used
EC$annot_rp3 <- EC$annot_new
EC$annot_rp3[colnames(PH_EC)] <- paste0(PH_EC$fun_annot, PH_EC$res_used)
EC = FindClusters(EC, resolution = 0.3)
pdf("/data4/lbl/PH/Figure/Supplement/SupFig3A.pdf")
DimPlot(EC, reduction = "tsne", label = T) + coord_equal()
dev.off()
EC$species <- factor(EC$species, levels = c("human", "rat", "mouse"))
pdf("/data4/lbl/PH/Figure/Supplement/SupFig3B.pdf", width = 10,height = 8)
DimPlot(EC, reduction = "tsne", split.by = "species", label = T, group.by = "annot_rp", ncol = 2) + NoLegend()
dev.off()
DimPlot(EC, reduction = "tsne", group.by = "annot_rp")

setwd("/data4/lbl/PH/")
EC$combine <- paste0("combined", EC$integrated_snn_res.0.3)
tmp1 <- subset(EC, species == "human")
tmp2 <- subset(EC, species != "human")
df1 <- as.data.frame(table(tmp1$combine, tmp1$annot_rp))
K <- length(unique(tmp1$combine))
N <- length(unique(tmp1$annot_rp))
end_idx <- 0 
for(i in 1:N){
  start_idx <- end_idx+1
  end_idx <- i*K
  df1$Freq[start_idx:end_idx] <- (df1$Freq[start_idx:end_idx]/sum(df1$Freq[start_idx:end_idx]))
}
df1 <- df1[,c(2,1,3)]
names(df1) <- c("Source", "Target", "Weight")
df2 <- as.data.frame(table(tmp2$combine, tmp2$annot_rp))
K <- length(unique(tmp2$combine))
N <- length(unique(tmp2$annot_rp))
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
#Fig5E
sankeyNetwork(Links=Sankeydata,Nodes=Sankeyname, Source ="Source",
              Target = "Target", Value = "Value", NodeID = "name",
              #colourScale = JS(ColourScale),
              # colourScale = JS('d3.scaleOrdinal().domain(["Fristående","Nästan_öppet","Halvöppet","Slutet"]).range(["#EDF8E9","#BAE4B3","#74C476","#238B45"])'),
              units = "Quads", fontSize = 12, nodeWidth = 30)


#Figure in discussion
angio_GO <- c("GO:0043066", "GO:0042127", "GO:0045596", "GO:0043069", "GO:0001525")
circv_GO <- c("GO:0042692", "GO:0003013")
GO <- EC_fun_annot.list$`Angiogenesis cEC`
GO@result <- GO@result[angio_GO,]
pdf("/data4/lbl/PH/Figure/Supplement/SupFig7A.pdf", width = 10)
cnetplot(GO,circular = TRUE, colorEdge = TRUE, colorDot = "pvalue", foldChange = geneList)+coord_equal()+ggtitle("Angiogenesis EC")
dev.off()