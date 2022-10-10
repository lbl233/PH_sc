library(Seurat)
library(velocyto.R)

file_dir <- "/data3/zhouyt/PAH-cellranger/"
doublet_dir <- "/data3/zhouyt/data/PAH/"
sample_name <- c("E4", "E5", "E6", "E7", "E8")
group <- c(rep("control",2),rep("Hypoxia", 3))
species <- c(rep("mouse", 5))
Batch <- c(rep("E",5))
dataset_list <- list()

for (i in 1:length(sample_name)) {
  filepath = paste(file_dir, paste(sample_name[i],sample_name[i], sep = "/"), sep = "")
  doubletpath = paste(doublet_dir , sample_name[i], sep = "")
  # Seurat
  seur.data = Read10X(data.dir = doubletpath)
  seur = CreateSeuratObject(counts = seur.data, project = sample_name[i])
  print(ncol(seur))
  head(colnames(seur))
  ldat <- read.loom.matrices(file = paste(filepath, "/velocyto/", sample_name[i], ".loom", sep = ""))
  emat <- ldat$spliced
  nmat <- ldat$unspliced
  amat <- ldat$ambiguous
  print(ncol(emat))
  colnames(emat) <- paste(substring(colnames(emat),4,19),"-1",sep="")
  colnames(nmat) <- paste(substring(colnames(nmat),4,19),"-1",sep="")
  colnames(amat) <- paste(substring(colnames(amat),4,19),"-1",sep="")
  head(colnames(emat))
  emat<-emat[!duplicated(rownames(emat)),]
  nmat<-nmat[!duplicated(rownames(nmat)),]
  amat<-amat[!duplicated(rownames(amat)),]
#  emat <- filter.genes.by.cluster.expression(emat,cluster.label,min.max.cluster.average = 0.5)
#  nmat <- filter.genes.by.cluster.expression(nmat,cluster.label,min.max.cluster.average = 0.05)
  length(intersect(rownames(emat),rownames(emat)))
  seur[["spliced"]] <- CreateAssayObject(counts = emat)
  seur[["unspliced"]] <- CreateAssayObject(counts = nmat)
  seur[["ambiguous"]] <- CreateAssayObject(counts = amat)
  doublet = read.table(paste(doubletpath, '/predicted_doublets.txt', sep = ""))
  doublet_score = read.table(paste(doubletpath, '/doublet_scores.txt', sep = ""))
  seur[['doublet_score']] = doublet_score$V1
  seur[['doublet']] = doublet$V1
  seur[['Group']] = group[i]
  seur[['Species']] = species[i]
  seur[['Batch']] = Batch[i]
  seur <- subset(seur, doublet == 0 & nFeature_RNA < 5000)
  dataset_list[[i]] <- seur
}

saveRDS(dataset_list, "~/Mydata/PAH/Hypoxia.rds")
