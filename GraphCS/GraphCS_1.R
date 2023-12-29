library(Seurat)
library(data.table)
library(stringr)
library(dplyr)
dir='../input'

Immu_Rdata_path = commandArgs(trailingOnly = TRUE)[1]
RNA_Rdata_path = commandArgs(trailingOnly = TRUE)[2]
species = commandArgs(trailingOnly = TRUE)[3]
ref = commandArgs(trailingOnly = TRUE)[4]
print(Immu_Rdata_path)
print(RNA_Rdata_path)
print(species)

#samples_folder=readLines('samples_path.txt')
samples=list.dirs("../input/",full.names=F,recursive=F)
print(samples)

gsmList = lapply(samples,function(pro){
  folder=file.path(dir ,pro)
  #print(pro)
  print(folder)
  print(pro)
  gsm.data = Read10X(folder)
  if (typeof(gsm.data) == 'list'){
  gsm=CreateSeuratObject(counts = gsm.data$'Gene Expression',
                         project =  pro ,
                         min.cells = 5,
                         min.features = 300) 
  }else if (typeof(gsm.data) == 'S4'){
  gsm=CreateSeuratObject(counts = gsm.data,
                         project =  pro ,
                         min.cells = 5,
                         min.features = 300) 
  }
  return(gsm)
})

names(gsmList)=samples
names(gsmList)
######################################################
if (length(gsmList) != 1){
for (i in 1:length(gsmList)) {
  gsmList[[i]] <- NormalizeData(gsmList[[i]], verbose = TRUE)
  gsmList[[i]] <- FindVariableFeatures(gsmList[[i]], selection.method = "vst", 
                                             nfeatures = 4000, verbose = TRUE)
}
gsmList

#########################质量控制（过滤低质量细胞-被线粒体污染的细胞）###########################
library(ggplot2)
for (i in 1:length(gsmList)) {
  mt<-PercentageFeatureSet(gsmList[[i]],pattern = "^MT-")
  gsmList[[i]] <- subset(gsmList[[i]], subset = nFeature_RNA > 200 & nFeature_RNA < 4500 & mt < 10)
  gsmList[[i]] <- subset(gsmList[[i]], features = head(VariableFeatures(gsmList[[i]]),4000))
}
rm(mt)
gc()
library(harmony)
#SCT标准化数据
gsmList <- lapply(X = gsmList, FUN = SCTransform, method = "glmGamPoi")
features <- SelectIntegrationFeatures(object.list = gsmList, nfeatures = 3000)
gsmList <- lapply(X = gsmList, FUN = RunPCA, features = features)
#数据缩放（去除基因效应）
gsmList <- lapply(X = gsmList,FUN = ScaleData)
gsm <- merge(gsmList[[1]], gsmList[2:length(gsmList)],add.cell.ids = samples)
#harmony前置需要pca
gsm <- RunPCA(gsm, npcs=50,features = features, verbose = TRUE)
#使用harmony整合数据
gsm <- RunHarmony(gsm,group.by.vars="orig.ident",assat.use="SCT",max.iter.harmony = 10)
#########降维聚类分群##########
library(cowplot)
#非线性降维聚类
gsm <- RunTSNE(gsm, reduction = "harmony", dims = 1:15) 
gsm <- RunUMAP(gsm, reduction = "harmony", dims = 1:15)
###找出类群###
gsm <- FindNeighbors(gsm, reduction = "harmony", dims = 1:15)
}else{
  gsm <- gsmList[[1]]
  gsm <- NormalizeData(gsm, verbose = TRUE)
  gsm <- FindVariableFeatures(gsm, selection.method = "vst", 
                                   nfeatures = 10000, verbose = TRUE)
gsm
library(ggplot2)
  mt<-PercentageFeatureSet(gsm,pattern = "^MT-")
  gsm <- subset(gsm, subset = nFeature_RNA > 200 & nFeature_RNA < 4500 & mt < 10)
  gsm <- subset(gsm, features = head(VariableFeatures(gsm),10000))
rm(mt)
gc()
library(harmony)
#SCT标准化数据
gsm <- SCTransform(gsm, method = "glmGamPoi")
features <- VariableFeatures(object = gsm)
gsm <- RunPCA(gsm,npcs=50,features = features, verbose = TRUE)
#数据缩放（去除基因效应）
gsm <- ScaleData(gsm)
rm(gsmList)
gc()
#harmony前置需要pca
gsm <- RunPCA(gsm, npcs=50,features = features, verbose = TRUE)
#########降维聚类分群##########
library(cowplot)
#非线性降维聚类
gsm <- RunTSNE(gsm, reduction = "pca", dims = 1:15) 
gsm <- RunUMAP(gsm, reduction = "pca", dims = 1:15)
###找出类群###
gsm <- FindNeighbors(gsm, reduction = "pca", dims = 1:15)
}#########单个List##########
rm(gsmList)
gc()
gsm <- FindClusters(gsm, resolution = 0.8)
table(gsm@meta.data$SCT_snn_res.0.8) 

#########添加标注信息用作训练集###########
library(SingleR)
# library(celldex)

gsm_for_SingleR <- GetAssayData(gsm,slot="data")

######################################

clusters=gsm@meta.data$seurat_clusters

if (species == "Human")
{
load(Immu_Rdata_path)
pred.humanImmu <- SingleR(test = gsm_for_SingleR, ref = humanImmu, labels = humanImmu$label.main,
                          method = "cluster", clusters = clusters, 
                          assay.type.test = "logcounts", assay.type.ref = "logcounts")


load(RNA_Rdata_path)

pred.humanRNA <- SingleR(test = gsm_for_SingleR, ref = humanRNA, labels = humanRNA$label.fine,
                         method = "cluster", clusters = clusters, 
                         assay.type.test = "logcounts", assay.type.ref = "logcounts")
 
CellType=data.frame(ClusterID=levels(gsm),
                    humanImmu=pred.humanImmu$labels,
                    humanRNA=pred.humanRNA$labels )
head(CellType)

gsm@meta.data$CellType=CellType[match(clusters,CellType$ClusterID),'humanRNA']
}else if (species == "Mouse"){
load(Immu_Rdata_path)
pred.mouseImmu <- SingleR(test = gsm_for_SingleR, ref = mouseImmu, labels = mouseImmu$label.main,
                          method = "cluster", clusters = clusters, 
                          assay.type.test = "logcounts", assay.type.ref = "logcounts")


load(RNA_Rdata_path)

pred.mouseRNA <- SingleR(test = gsm_for_SingleR, ref = mouseRNA, labels = mouseRNA$label.fine,
                         method = "cluster", clusters = clusters, 
                         assay.type.test = "logcounts", assay.type.ref = "logcounts")
 
CellType=data.frame(ClusterID=levels(gsm),
                    mouseImmu=pred.mouseImmu$labels,
                    mouseRNA=pred.mouseRNA$labels )
head(CellType)

gsm@meta.data$CellType=CellType[match(clusters,CellType$ClusterID),'mouseRNA']
}
save(gsm,Immu_Rdata_path,RNA_Rdata_path,species,file = 'output.Rdata')
pro='first_anno'

gsm@meta.data$batchlb<-gsm@meta.data$orig.ident

gsm2 <- gsm

gsmList = lapply(c(gsm2,gsm),function(pro){
  list=pro
  return(list)
})
save(gsmList,file = 'outputList.RData')


######################useless###################
#if (species == "Human")
#{
#load(humanref)
#gsm1 <- FindVariableFeatures(gsm@assays$RNA, selection.method = "vst", nfeatures = 5000)
#gsm2 <- subset(x = gsm, features = head(VariableFeatures(gsm1),5000))

#gsm2

#gsmList = lapply(c(schuman,gsm),function(pro){
#  list=pro
#  return(list)
#})
#save(gsmList,file = 'outputList.RData')
#}else if (species == "Mouse"){

#load(mouseref)

#gsm1 <- FindVariableFeatures(gsm@assays$RNA, selection.method = "vst", nfeatures = 2000)
#gsm2 <- subset(x = gsm, features = head(VariableFeatures(gsm1),2000))
#gsm2
#gsmList = lapply(c(scmouse,gsm),function(pro){
#  list=pro
#  return(list)
#})
#save(gsmList,file = 'outputList.RData')
#}