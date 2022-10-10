setwd("~/Mydata/PAH/")
library(Seurat)
library(dplyr)
library(ggplot2)

#load preprocessed data
load("/home/zhouyt/R/PAH_sc/image/2.RData")
load("/home/zhouyt/R/marker.v2.RData")
load("/home/zhouyt/R/PAH4/annot_metadata.RData")

old.barcode <- intersect(rownames(temp), colnames(S))
S$zyt_annot <- "unannot"

S$zyt_annot[old.barcode] <- as.character(temp[old.barcode,1])
S$zyt_annot2 <- "unannot"
S$zyt_annot2[old.barcode] <- as.character(temp[old.barcode,2])

save.image("20201213.RData")
DimPlot(S, label = T)
DimPlot(S, split.by = "id1", ncol = 4)
table(S$batch)
DimPlot(S, group.by = "zyt_annot", split.by = "id1")
DimPlot(S, group.by = "zyt_annot", label = T)
FeaturePlot(S, features = c("Ptprc", "Pecam1", "Col1a2", #identify immune cell and non-immune cell
                            "Cd68", "Adgre1", "Cd14", "Itgax", "Itgam", "Fcgr2b", "Siglecf", #macrophage and myeloid markers
                            "Cd3e", "Cd4", "Cd8a", "Klrd1", "Foxp3", #T cell markers
                            "Cd19", "Mzb1", "Cd79a", "Ighm", # B cell markers
                            "Vwf", "Kit", "Vcam1", "Plvap", 'Itga2b', "Pf4", "Itgb3", #EC and MKC marker
                            "Cd34", "Acta2", "Pdgfra", "Myh11", "Cd248", #SMC and FB marker
                            "Cpa3", "Fcer1a", #MC markers
                            "Sftpc", "Hopx", "Foxj1"), cols=c("darkblue", "goldenrod","firebrick"), ncol = 6)
pdf("featureplot_marker_gene.pdf", width = 20, height = 20)
dev.off()
##### myeloid cell identification #####
#Ly6g and Cd11c are the general markers for Neutrophil, and Cd15 and Cd16 are used to distinguish Neutro from monocyte and eos
#F4/80 and MerTK are used to identify macrophage
FeaturePlot(S, features = c("Ptprc",  #CD45
                            "Cd68", "Il1b",
                            "Adgre1", #F4/80 cDC dont express
                            "Cd14", "Fcgr3", # CD14, CD16
                            "Itgax", "Itgam", # CD11c, CD11b
                            "Siglecf", #specific to eos and Alvelor macrophage
                            "Itgae", "Cd8a", #CD103
                            "Ly6g", "Mertk", "Ly6c2", #monocyte markers
                            "Csf2rb2", "Bst2", #IL3R Bst-2 pDC marker
                            "Itga2b", "Cd200r3", "Cebpa", #basophil markers, lack expression of c-kit
                            "Fcer1a", "Kit", "Cpa3" # inflammatory DC and MC marker
                            ), cols=c("darkblue", "goldenrod","firebrick"), ncol = 5)
