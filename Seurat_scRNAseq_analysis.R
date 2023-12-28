library(Seurat)
library(dplyr)
library(harmony)

####1-CTRL####
CTRL=Read10X(data.dir = "E:/WWT_TYY_CD34/CTRL/filtered_feature_bc_matrix/")
CTRL=CreateSeuratObject(counts = CTRL, project = "HSC_CTRL", min.cells = 20, min.features = 200)

CTRL[["percent.mt"]] <- PercentageFeatureSet(CTRL, pattern = "^MT-")
ggplot(CTRL@meta.data,aes(x=nCount_RNA,y=percent.mt))+geom_point(size=1,color='navy',shape=20)+scale_x_continuous(breaks=seq(0,170000,5000))+scale_y_continuous(breaks=seq(0,25,1))+theme(axis.text.x= element_text(angle=90))
ggplot(CTRL@meta.data,aes(x=nCount_RNA,y=nFeature_RNA))+geom_point(size=1,color='navy',shape=20)+scale_x_continuous(breaks=seq(0,170000,5000))+scale_y_continuous(breaks=seq(0,9000,250))+theme(axis.text.x= element_text(angle=90))
CTRL <- subset(CTRL, subset = nFeature_RNA > 1750 & nFeature_RNA < 5500 & percent.mt < 3.6& nCount_RNA<30000& nCount_RNA>0)

  ##Normalize
CTRL <- NormalizeData(CTRL, normalization.method = "LogNormalize", scale.factor = 10000)
##Identification of highly variable features (feature selection)
CTRL <- FindVariableFeatures(CTRL, selection.method = "vst", nfeatures = 1500)

CTRL$group<- "CTRL"
CTRL=RenameCells(CTRL,add.cell.id = "CTRL")

####2-TGFB####
TGFB=Read10X(data.dir = "E:/WWT_TYY_CD34/TGFb/filtered_feature_bc_matrix/")
TGFB=CreateSeuratObject(counts = TGFB, project = "HSC_TGFB", min.cells = 20, min.features = 200)

TGFB[["percent.mt"]] <- PercentageFeatureSet(TGFB, pattern = "^MT-")
ggplot(TGFB@meta.data,aes(x=nCount_RNA,y=percent.mt))+geom_point(size=1,color='navy',shape=20)+scale_x_continuous(breaks=seq(0,170000,5000))+scale_y_continuous(breaks=seq(0,25,1))+theme(axis.text.x= element_text(angle=90))
ggplot(TGFB@meta.data,aes(x=nCount_RNA,y=nFeature_RNA))+geom_point(size=1,color='navy',shape=20)+scale_x_continuous(breaks=seq(0,170000,5000))+scale_y_continuous(breaks=seq(0,9000,250))+theme(axis.text.x= element_text(angle=90))
TGFB <- subset(TGFB, subset = nFeature_RNA > 1750 & nFeature_RNA < 5750 & percent.mt < 4 & nCount_RNA<30000& nCount_RNA>0)

  ##Normalize
TGFB <- NormalizeData(TGFB, normalization.method = "LogNormalize", scale.factor = 10000)
##Identification of highly variable features (feature selection)
TGFB <- FindVariableFeatures(TGFB, selection.method = "vst", nfeatures = 1500)

TGFB$group<- "TGFB"
TGFB=RenameCells(TGFB,add.cell.id = "TGFB")

####3-Harmony integrated, batch effect removal####  
Allsample=merge(CTRL,TGFB)
Allsample <- FindVariableFeatures(Allsample, selection.method = "vst", nfeatures = 1500)
     
  ##regress out cell cycle
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
Allsample <- CellCycleScoring(Allsample, s.features = s.genes, g2m.features = g2m.genes)
Allsample <- ScaleData(Allsample, vars.to.regress =c( "percent.mt",'nCount_RNA',"S.Score", "G2M.Score"))

Allsample <- RunPCA(Allsample, npcs = 20, verbose = FALSE)
Allsample=RunHarmony(Allsample,group.by.vars = 'group', plot_convergence = FALSE,dims.use = 1:20,assay.use = 'RNA')

Allsample=RunUMAP(Allsample,reduction = "harmony", dims = 1:20)
Allsample=FindNeighbors(Allsample,reduction = "harmony", dims = 1:20)
Allsample=FindClusters(Allsample,resolution = 0.6)
DimPlot(Allsample, reduction = "umap", label = TRUE, pt.size = 1.3)

marker=FindAllMarkers(Allsample,logfc.threshold = 0.25 ,only.pos = T)


