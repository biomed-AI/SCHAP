#scATAC-seq细胞类型标注
library(Signac)
library(Seurat)
library(GenomeInfoDb)
#library(EnsDb.Hsapiens.v75)
#library(EnsDb.Hsapiens.v86)
library(patchwork)
library(dplyr)
library(ChIPpeakAnno)
library(ChIPseeker)
library(clusterProfiler)
library(GenomicFeatures) 
#library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#library(org.Hs.eg.db)
library(ggplot2)
library(SingleCellExperiment)
set.seed(1234)

species = commandArgs(trailingOnly = TRUE)[1]
anno_data = commandArgs(trailingOnly = TRUE)[2]
FragmentHistogram_region = commandArgs(trailingOnly = TRUE)[3]

print(species)

if (species == "Human")
{
  library(EnsDb.Hsapiens.v86)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  genomename = 'hg38'
}else{
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(EnsDb.Mmusculus.v79)
  library(org.Mm.eg.db)
  genomename = 'mm10'
}




h5=list.files('../input/sample/','.h5',full.names=T,recursive=F)
csv=list.files('../input/sample/','.csv',full.names=T,recursive=F)
tsv=list.files('../input/sample/','.tsv.gz$',full.names=T,recursive=F)

print(h5)
print(csv)
print(tsv)

# file1-读取Peak/Cell矩
#Peak/Cell matrix. 该文件类似于用于分析单细胞RNA-seq的基因表达计数矩阵。然而，矩阵的每一行不是基因，而是代表一个基因组区域（“peak region”，一个峰），
counts <- Read10X_h5(filename = h5)

# 查看peaks和细胞的数目
#dim(counts)
#counts[1:5,1:5]

# file2-读取细胞注释信息
# 该文件储存了所有单个细胞中唯一片段（unique fragments）的完整列表。该文件非常大，处理速度较慢，并且存储在磁盘上（而不是内存中）。
# 但是，保留此文件的优点是它包含与每个细胞相关的所有片段，而不是仅映射到峰的片段。有关片段文件的更多信息
metadata <- read.csv(
  file = csv,
  header = TRUE,
  row.names = 1
)

# 查看metadata信息
#names(metadata)
#head(metadata)[,1:5]


# ATAC-seq数据使用自定义的assay— ChromatinAssay 进行存储。这使得一些专门的功能可以用于分析单细胞的基因组assay，
# 创建ChromatinAssay
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = genomename,
  fragments = tsv,
  min.cells = 10,
  min.features = 200
)

# 构建seurat对象
pbmc <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)

pbmc
#pbmc[['peaks']]
#granges(pbmc)


# 可以为pbmc对象的人类基因组信息添加基因注释。这使得下游函数直接从对象中提取基因注释信息。
#annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
#annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)
load(anno_data)
annotations
# change to UCSC style since the data was mapped to hg19
#seqlevelsStyle(annotations) <- 'UCSC'
# add the gene information to the object
Annotation(pbmc) <- annotations

pbmc <- NucleosomeSignal(object = pbmc)

# compute TSS enrichment score per cell
pbmc <- TSSEnrichment(object = pbmc, fast = FALSE)

# add blacklist ratio and fraction of reads in peaks
pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments / pbmc$passed_filters * 100

pbmc$blacklist_ratio <- pbmc$blacklist_region_fragments / pbmc$peak_region_fragments

pbmc$high.tss <- ifelse(pbmc$TSS.enrichment > 2, 'High', 'Low')

pdf("TSS_enrichment.pdf")
TSSPlot(pbmc, group.by = 'high.tss') + NoLegend()
dev.off()


pdf("violinplot.pdf",width=15, height=5)
VlnPlot(
  object = pbmc,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)
dev.off()

# 根据不同QC指标进行数据过滤
pbmc <- subset(
  x = pbmc,
  subset = peak_region_fragments > 3000 &
    peak_region_fragments < 20000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)
pbmc

# 使用RunTFIDF函数进行数据归一化
pbmc <- RunTFIDF(pbmc)
# 使用FindTopFeatures函数选择可变的features
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
# 使用RunSVD函数进行线性降维
pbmc <- RunSVD(pbmc)

pdf("correlation_depth_components.pdf")
DepthCor(pbmc)
dev.off()


#非线性降维和聚类
pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindNeighbors(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindClusters(object = pbmc, verbose = FALSE, algorithm = 3)
pdf("UMAP.pdf")
DimPlot(object = pbmc, label = TRUE) + NoLegend()
dev.off()


pbmc$nucleosome_group <- ifelse(pbmc$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
pdf("Fragment_length.pdf")
FragmentHistogram(object = pbmc, group.by = 'nucleosome_group',region = FragmentHistogram_region)
dev.off()



## 创建基因表达活性矩阵（gene activity matrix）
#gene.activities <- GeneActivity(pbmc)
# add the gene activity matrix to the Seurat object as a new assay and normalize it 
#pbmc[['RNA']] <- CreateAssayObject(counts = gene.activities) 
#pbmc <- NormalizeData( 
#  object = pbmc, 
#  assay = 'RNA', 
#  normalization.method = 'LogNormalize', 
#  scale.factor = median(pbmc$nCount_RNA) 
#)

#DefaultAssay(pbmc) <- 'RNA' 

#pdf("FeaturePlot.pdf",width=12, height=6)
#FeaturePlot( 
#  object = pbmc, 
#  features = c('MS4A1', 'CD3D', 'LEF1', 'NKG7', 'TREM1', 'LYZ'), 
#  pt.size = 0.1, 
#  max.cutoff = 'q95', 
#  ncol = 3 
#)
#dev.off()