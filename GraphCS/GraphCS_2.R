options(stringsAsFactors = F)
library(Seurat)
library(stringr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(data.table)

Rdata = commandArgs(trailingOnly = TRUE)[1]
csv = commandArgs(trailingOnly = TRUE)[2]

pvalueCutoff = as.numeric(commandArgs(trailingOnly = TRUE)[3])
qvalueCutoff = as.numeric(commandArgs(trailingOnly = TRUE)[4])
heatmap_width = as.numeric(commandArgs(trailingOnly = TRUE)[5])
heatmap_height = as.numeric(commandArgs(trailingOnly = TRUE)[6])
bubble_width = as.numeric(commandArgs(trailingOnly = TRUE)[7])
bubble_height = as.numeric(commandArgs(trailingOnly = TRUE)[8])

load(Rdata)
gsm
RNA <- read.csv(csv ,header = T,row.names=1)

# marker gene #
###该步骤应该作为一个表格输出，预计需要输出更多marker基因？比如20个
gsm <- PrepSCTFindMarkers(object=gsm ,assay="SCT" ,verbose=TRUE)
gsm.markers <- FindAllMarkers(gsm, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- gsm.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top10


#############################################
#细胞注释(人和小鼠都可以)
#这里目前是使用聚类里出现频次最高的细胞名来注释聚类，后续需要优化的
gsm@meta.data$CellType <- RNA[,1]

CellType=data.frame(ClusterID=gsm@meta.data$seurat_clusters,
                    RNA=gsm@meta.data$CellType)

celltype = lapply(levels(gsm@meta.data$seurat_clusters),function(i){
  string_counts <- table(CellType[CellType$ClusterID==i,2])
  celltype <- names(string_counts)[which.max(string_counts)]
  return(celltype)
})
CellType <- paste(celltype)

CellType=data.frame(ClusterID=levels(gsm),RNA=CellType)

clusters=gsm@meta.data$seurat_clusters
gsm@meta.data$CellType=CellType[match(clusters,CellType$ClusterID),'RNA']

gsm_fin_umap <- DimPlot(gsm,reduction = "umap",label=T, group.by = 'CellType')
ggsave(gsm_fin_umap,filename = "outputs_umap.pdf",width = 10, height = 7, units = 'in', dpi = 300 )
gsm_fin_tsne <- DimPlot(gsm,reduction = "tsne",label=T, group.by = 'CellType')
ggsave(gsm_fin_tsne,filename = "outputs_tsne.pdf",width = 10, height = 7, units = 'in', dpi = 300 )

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


#用AverageHeatmap画关于marker基因的图
library(dplyr)
library(scRNAtoolVis)
library(Cairo)

Seurat::Idents(gsm) <- gsm$CellType

CairoPDF(file="outputs_heatmap.pdf",width = heatmap_width,height = heatmap_height)
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

##heatmap signal##
CairoPDF(file="cell_heat_signal.pdf",width = 10,height = 8)
par(mfrow = c(1,2), xpd=TRUE)
h1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing",height = 10)
h2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming",height = 10)
h1+h2
dev.off()

pathways.show <- cellchat@netP$pathways
CairoPDF(file="SignalingRole.pdf",width = 10,height = 6)
p <- netAnalysis_signalingRole_network(cellchat, signaling= pathways.show, width = 10,
                                  height = 6, font.size = 10)
dev.off()