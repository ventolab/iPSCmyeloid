
#R 3.6.0

library(Seurat)
library(viridis)


'%!in%' <- function(x,y)!('%in%'(x,y))
sufa<-function(x){summary(as.factor(x))}

############################################
############### Import data ################
############################################

### Import metadata from the processed object dowloaded from  www.HiPImmuneatlas.org
AGGprocd<-readRDS("path/to/dowloaded/RDSobject.rds")
MD<-AGGprocd@meta.data
rm(AGGprocd);gc()

### Import cellranger aggr output matrix
path<-"path/to/aggr/folder"
AGG.data <- Read10X(data.dir = paste0(path,"/outs/filtered_feature_bc_matrix"))
AGG <- CreateSeuratObject(counts = AGG.data, project = "AGG")

### Transfer metadata variables into aggregated object
AGG@meta.data$VARIABLEname<-MD[match(rownames(AGG@meta.data),rownames(MD)),"VariableNameInMD"] #CellLineSOUPorCELL will be needed


###################################
############### QC ################
###################################


### Calculate percent of mitochondrial RNA
AGG[["percent.mt"]] <- PercentageFeatureSet(AGG, pattern = "^MT-")

### Check variables for QC
plot1 <- FeatureScatter(AGG, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(AGG, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

VlnPlot(AGG, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

### Subset low QC cells
AGG <- subset(AGG, subset =  nFeature_RNA > 200  & percent.mt < 8.5 & MergedDoublets != "DBL") #shown for completeness, just remove cells NA in CellLineSOUPorCELL variable


### Annotate cell cycle phase
cc.genes[cc.genes=="MLF1IP"]<-"CENPU" #there used to be a warning of 3 genes not found in object, current synonyms are there
cc.genes[cc.genes=="FAM64A"]<-"TENT5A"
cc.genes[cc.genes=="HN1"]<-"JPT1"
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]
AGG <- CellCycleScoring(object = AGG, s.features = s.genes, g2m.features = g2m.genes)
AGG$CC.Difference <- AGG$S.Score - AGG$G2M.Score


RidgePlot(AGG, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2, group.by = "Phase")

#################################################
############### Batch Correction ################
#################################################

Idents(AGG)<-"CellLineSOUPorCELL"

AGG.list <- SplitObject(AGG, split.by = "ident")

for (i in 1:length(AGG.list)) {
  AGG.list[[i]] <- SCTransform(AGG.list[[i]], vars.to.regress = c("nCount_RNA","percent.mt","CC.Difference"),verbose = T)
}

AGG.features <- SelectIntegrationFeatures(object.list = AGG.list, nfeatures = 5000)
AGG.list <- PrepSCTIntegration(object.list = AGG.list, anchor.features = AGG.features, verbose = T)
AGG.anchors <- FindIntegrationAnchors(object.list = AGG.list, normalization.method = "SCT", anchor.features = AGG.features, verbose = T)
AGG <- IntegrateData(anchorset = AGG.anchors, normalization.method = "SCT", verbose = T)

#########################################################
############### Dimensionality reduction ################
#########################################################

AGG <- RunPCA(AGG, verbose = T)
str(AGG[["pca"]])
eigs<-AGG[["pca"]]@stdev^2
AGG[["pca"]]@misc$proportion<-eigs/sum(eigs)
AGG[["pca"]]@misc$cumulative<-cumsum(eigs)/sum(eigs)

SelPC<-min(which(AGG[["pca"]]@misc$cumulative>0.9)) #PCs that explain 90% of variability

AGG <- RunUMAP(AGG, dims = 1:SelPC, verbose = T,n.components = 3L,reduction.name="umap3",reduction.key="UMAP3L_") #2D UMAP plots
AGG <- RunUMAP(AGG, dims = 1:SelPC, verbose = T,n.components = 2L,reduction.name="umap",reduction.key="UMAP2L_") #3D UMAP plots

AGG <- FindNeighbors(AGG, dims = 1:SelPC, verbose = T)
RES<-c(0.5,0.8,1,3,6,10)
AGG <- FindClusters(AGG,resolution=RES, verbose = T)


#######################################################
############### Store processed object ################
#######################################################

saveRDS(AGG, file="path/to/folder/SeuratObjectName.rds")


