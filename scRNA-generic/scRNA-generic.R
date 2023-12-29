library(Seurat)
library(data.table)
library(stringr)
library(dplyr)
dir='../input'

Immu_Rdata_path = commandArgs(trailingOnly = TRUE)[1]
RNA_Rdata_path = commandArgs(trailingOnly = TRUE)[2]
species = commandArgs(trailingOnly = TRUE)[3]

pvalueCutoff = as.numeric(commandArgs(trailingOnly = TRUE)[4])
qvalueCutoff = as.numeric(commandArgs(trailingOnly = TRUE)[5])
heatmap_width = as.numeric(commandArgs(trailingOnly = TRUE)[6])
heatmap_height = as.numeric(commandArgs(trailingOnly = TRUE)[7])
bubble_width = as.numeric(commandArgs(trailingOnly = TRUE)[8])
bubble_height = as.numeric(commandArgs(trailingOnly = TRUE)[9])

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
                                             nfeatures = 2000, verbose = TRUE)
}
gsmList

#########################质量控制（过滤低质量细胞-被线粒体污染的细胞）###########################
library(ggplot2)
for (i in 1:length(gsmList)) {
  mt<-PercentageFeatureSet(gsmList[[i]],pattern = "^MT-")
  gsmList[[i]] <- subset(gsmList[[i]], subset = nFeature_RNA > 200 & nFeature_RNA < 4500 & mt < 10)
  gsmList[[i]] <- subset(gsmList[[i]], features = head(VariableFeatures(gsmList[[i]]),2000))
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
gsm
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
#最后对聚类的不同细胞亚群进行注释
#前面呢是标准的聚类分群，每个细胞亚群仅仅是一个编号，实际上还需要给予它们一定的生物学意义，我们这里采用SingleR的标准代码：
options(stringsAsFactors = F)
library(Seurat)
library(stringr)
library(dplyr)
library(ggplot2)

# marker gene #
###该步骤应该作为一个表格输出，预计需要输出更多marker基因？比如20个
gsm <- PrepSCTFindMarkers(object=gsm ,assay="SCT" ,verbose=TRUE)
gsm.markers <- FindAllMarkers(gsm, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- gsm.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top10

#############################################
#基因富集（GO和KEGG）

library(clusterProfiler)
options(clusterProfiler.download.method = 'curl')
options(download.file.extra = '-k')
library(topGO)
library(Rgraphviz)
library(GOSemSim)
library(DOSE)
library(GOplot)
library(stringr)
library(enrichplot)


if (species == "Human")
{
library(org.Hs.eg.db)
genelist <- as.character(top10$gene)

eg <- bitr(genelist,fromType="SYMBOL", 
           toType=c("ENTREZID","GENENAME"), 
           OrgDb="org.Hs.eg.db"); 
head(eg)

geneList <- eg$ENTREZID

#查看基因在聚类里的分布
#输入能被map的marker基因
gene <- eg$SYMBOL

go <- enrichGO(gene = geneList,
               OrgDb = org.Hs.eg.db,
               ont='BP',            # replace with BP, CC and MF
               pAdjustMethod = 'BH',
               pvalueCutoff = pvalueCutoff,
               qvalueCutoff = qvalueCutoff,
               keyType = 'ENTREZID')
head(go)
write.csv(go@result, file = "go_result.csv", row.names = F)
go_bar <- barplot(go, showCategory=10,drop=T) 
go_dot <- dotplot(go, showCategory=10) 

ggsave(go_bar,filename = "go_bar.pdf", dpi = 300 , height = 10)
ggsave(go_dot,filename = "go_dot.pdf", dpi = 300 , height = 10)

enrichKG <- enrichKEGG(geneList, 
            organism = 'hsa',
            keyType = 'kegg', 
            pvalueCutoff = pvalueCutoff, 
            pAdjustMethod = 'BH', 
            minGSSize = 3, 
            maxGSSize = 500, 
            qvalueCutoff = qvalueCutoff, 
            use_internal_data = F)
			
head(enrichKG)
write.csv(enrichKG@result, file = "kg_result.csv", row.names = F)

kg_bar <- barplot(enrichKG,showCategory=15)
kg_dot <- dotplot(enrichKG)

ggsave(kg_bar,filename = "kg_bar.pdf", dpi = 300 )
ggsave(kg_dot,filename = "kg_dot.pdf", dpi = 300 )
}else if (species == "Mouse"){
library(org.Mm.eg.db)
genelist <- as.character(top10$gene)

eg <- bitr(genelist,fromType="SYMBOL", 
           toType=c("ENTREZID","GENENAME"), 
           OrgDb="org.Mm.eg.db"); 
head(eg)

geneList <- eg$ENTREZID

#查看基因在聚类里的分布
#输入能被map的marker基因
gene <- eg$SYMBOL

go <- enrichGO(gene = geneList,
               OrgDb = org.Mm.eg.db,
               ont='BP',            # replace with BP, CC and MF
               pAdjustMethod = 'BH',
               pvalueCutoff = pvalueCutoff,
               qvalueCutoff = qvalueCutoff,
               keyType = 'ENTREZID')
head(go)
write.csv(go@result, file = "go_result.csv", row.names = F)
go_bar <- barplot(go, showCategory=10,drop=T) 
go_dot <- dotplot(go, showCategory=10) 

ggsave(go_bar,filename = "go_bar.pdf", dpi = 300 , height = 10)
ggsave(go_dot,filename = "go_dot.pdf", dpi = 300 , height = 10)

enrichKG <- enrichKEGG(geneList, 
            organism = 'mmu',
            keyType = 'kegg', 
            pvalueCutoff = pvalueCutoff, 
            pAdjustMethod = 'BH', 
            minGSSize = 3, 
            maxGSSize = 500, 
            qvalueCutoff = qvalueCutoff, 
            use_internal_data = F)
			
head(enrichKG)
write.csv(enrichKG@result, file = "kg_result.csv", row.names = F)

kg_bar <- barplot(enrichKG,showCategory=15)
kg_dot <- dotplot(enrichKG)

ggsave(kg_bar,filename = "kg_bar.pdf", dpi = 300 )
ggsave(kg_dot,filename = "kg_dot.pdf", dpi = 300 )
}

#############################################
#细胞注释(人和小鼠都可以)
#这里展示的是人的分析，小鼠对应导入的RNA数据名字是另一个，mouse***的

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
 
CellType=data.frame(ClusterID=levels(gsm@meta.data$seurat_clusters),
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
 
CellType=data.frame(ClusterID=levels(gsm@meta.data$seurat_clusters),
                    mouseImmu=pred.mouseImmu$labels,
                    mouseRNA=pred.mouseRNA$labels )
head(CellType)

gsm@meta.data$CellType=CellType[match(clusters,CellType$ClusterID),'mouseRNA']
}

pro='first_anno'

gsm_fin_umap <- DimPlot(gsm,reduction = "umap",label=T, group.by = 'CellType')
ggsave(gsm_fin_umap,filename = "umap.pdf",width = 10, height = 7, units = 'in', dpi = 300 )
gsm_fin_tsne <- DimPlot(gsm,reduction = "tsne",label=T, group.by = 'CellType')
ggsave(gsm_fin_tsne,filename = "tsne.pdf",width = 10, height = 7, units = 'in', dpi = 300 )

######提取非线性低维聚类信息######

tSNE<- as.data.frame(gsm@reductions$tsne@cell.embeddings)
UMAP<- as.data.frame(gsm@reductions$umap@cell.embeddings)
Seurat::Idents(gsm) <- gsm$CellType
celltype<- Idents(gsm)
tSNE<- cbind(tSNE,celltype)
UMAP<- cbind(UMAP,celltype)
write.csv(tSNE,file="tsne.csv")
write.csv(UMAP,file="umap.csv")
rm(tSNE,UMAP)
gc()

#用AverageHeatmap画关于marker基因的图
library(dplyr)
library(scRNAtoolVis)
library(Cairo)

Seurat::Idents(gsm) <- gsm$CellType

CairoPDF(file="heatmap.pdf",width = heatmap_width,height = heatmap_height)
AverageHeatmap(object = gsm, markerGene = head(gene,60),
               cluster_rows =T, row_names_side='right',
               border =T,row_km=length(levels(gsm)),column_split =1:length(levels(gsm)),
               column_title = NULL,row_title=NULL,assays='SCT'
               )
dev.off()


#差异分析，提取出每个细胞的差异基因，整合后画热图#

for(i in 1:length(levels(gsm))) {
data <- FindMarkers(gsm, ident.1 = gsm@meta.data$CellType[i], ident.2 = NULL, only.pos = TRUE, 
min.pct = 0.25, logfc.threshold = 0.25) 
datatop10 <- data %>% top_n(n = 10, wt = avg_log2FC)
if (i == 1){
  markers <- row.names(datatop10)
}else{
  markers <- rbind(markers,row.names(datatop10))
}
} 
markers <- as.character(markers)
CairoPDF(file="variation_heatmap.pdf",width = heatmap_width,height = heatmap_height)
AverageHeatmap(object = gsm, markerGene = markers,
               cluster_rows = T, row_names_side='right',
               border =T,row_km=length(levels(gsm)),column_split =1:length(levels(gsm)),
               column_title = NULL,row_title=NULL,assays='SCT'
               )
dev.off()


##############################################
#细胞通讯(使用cellchat)

library(CellChat)
library(patchwork)
library(Seurat)
#library(ggplot2)
options(stringsAsFactors = FALSE)

#做没做SCT的都能用这句话，因为会匹配默认的assay
data.input = GetAssayData(gsm,slot = "data")

meta.data =  gsm@meta.data
cellchat <- createCellChat(object = data.input, 
                           meta = meta.data, 
                           group.by = "CellType")
rm(data.input)
rm(meta.data)
rm(gsm_for_SingleR)
rm(gsm)
gc()
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
groupSize

#加载CellChatDB数据库

if (species == "Human")
{
CellChatDB <- CellChatDB.human
}else if (species == "Mouse"){
CellChatDB <- CellChatDB.mouse
}


dplyr::glimpse(CellChatDB$interaction)

CellChatDB.use <- CellChatDB # simply use the default CellChatDB
# set the used database in the object
cellchat@DB <- CellChatDB.use

cellchat <- subsetData(cellchat)  # This step is necessary even if using the whole database
library(future)
options(future.globals.maxSize= 40000000000) #大于3407872000即可
future::plan("multicore", workers = 10)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

#计算通讯概率，推断细胞通讯网络
cellchat <- computeCommunProb(cellchat,population.size = F)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

levels(cellchat@idents)
#在信号通路水平推断细胞通讯
cellchat <- computeCommunProbPathway(cellchat)
head(cellchat@net)
head(cellchat@netP)

#计算加和的cell-cell通讯网络
cellchat <- aggregateNet(cellchat)

#细胞互作数量与强度统计分析

groupSize <- as.numeric(table(cellchat@idents))
library(Cairo)
CairoPDF(file="TTL_net_number_strength.pdf",width = 12,height = 7)
par(mfrow = c(1,2), xpd=TRUE)
p1 <- netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale=T,
                  label.edge=F, arrow.size = 1, title.name = "Number of Interactions")
p2 <- netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale=T,
                  label.edge=F, arrow.size = 1,title.name = "Interactions weights/strength")
dev.off()
# save as TTL/net_number_strength.pdf

# show all the significant interactions (L-R pairs)
#需要指定受体细胞和配体细胞
p = netVisual_bubble(cellchat, sources.use = seq(1,length(levels(cellchat@idents)),2),
                     targets.use = seq(2,length(levels(cellchat@idents)),2),,remove.isolate = FALSE,angle.x = 45)
ggsave("bubble.pdf",p ,width = bubble_width,height = bubble_height)
# save as TIL/Mye Lymph bubble.pdf

#计算网络中心性权重
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name="netP")

pathways.show <- cellchat@netP$pathways

##heatmap signal##
CairoPDF(file="cell_heat_signal.pdf",width = 10,height = 8)
par(mfrow = c(1,2), xpd=TRUE)
h1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing",height = 10)
h2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming",height = 10)
h1+h2
dev.off()

CairoPDF(file="SignalingRole.pdf",width = 10,height = 6)
p <- netAnalysis_signalingRole_network(cellchat, signaling= pathways.show, width = 10,
                                  height = 6, font.size = 10)
dev.off()