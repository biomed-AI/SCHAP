library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(patchwork)
library(dplyr)
library(ChIPpeakAnno)
library(ChIPseeker)
library(clusterProfiler)
library(GenomicFeatures) 
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
rds=list.files('../input/sample/','.rds',full.names=T,recursive=F)

print(h5)
print(csv)
print(tsv)
print(rds)

## Peak matrix
counts <- Read10X_h5(filename = h5)

## metadata
metadata <- read.csv(
  file = csv,
  header = TRUE,
  row.names = 1
)

## ChromatinAssay
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = genomename,
  fragments = tsv,
  min.cells = 10,
  min.features = 200
)

gsm_atac <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)

gsm_atac

## add gene annotations
load(anno_data)
annotations
# annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
# change to UCSC style since the data was mapped to hg19
# seqlevelsStyle(annotations) <- 'UCSC'
# add the gene information to the object
Annotation(gsm_atac) <- annotations

## compute nucleosome signal score per cell
gsm_atac <- NucleosomeSignal(object = gsm_atac)
gsm_atac$nucleosome_group <- ifelse(gsm_atac$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
pdf("Fragment_length.pdf")
FragmentHistogram(object = gsm_atac, group.by = 'nucleosome_group',region = FragmentHistogram_region)
# FragmentHistogram(object = gsm_atac, group.by = 'nucleosome_group')
dev.off()

## compute Transcriptional start site (TSS) enrichment score
gsm_atac <- TSSEnrichment(object = gsm_atac, fast = FALSE)
gsm_atac$high.tss <- ifelse(gsm_atac$TSS.enrichment > 2, 'High', 'Low')
pdf("TSS_enrichment.pdf")
TSSPlot(gsm_atac, group.by = 'high.tss') + NoLegend()
dev.off()

## add blacklist ratio and fraction of reads in peaks
gsm_atac$pct_reads_in_peaks <- gsm_atac$peak_region_fragments / gsm_atac$passed_filters * 100  
gsm_atac$blacklist_ratio <- gsm_atac$blacklist_region_fragments / gsm_atac$peak_region_fragments  
pdf("violinplot.pdf",width=15, height=5)
VlnPlot(
  object = gsm_atac,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)
dev.off()

## remove cells that are outliers
gsm_atac <- subset(
  x = gsm_atac,
  subset = peak_region_fragments > 3000 &
    peak_region_fragments < 20000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)
gsm_atac

## Normalization and linear dimensional reduction
gsm_atac <- RunTFIDF(gsm_atac) # （Normalization）
gsm_atac <- FindTopFeatures(gsm_atac, min.cutoff = 'q0 ') # (Feature selection)
gsm_atac <- RunSVD(gsm_atac) # (Dimension reduction)
pdf("correlation_depth_components.pdf")
DepthCor(gsm_atac) 
dev.off()

##  Non-linear dimension reduction and clustering
gsm_atac <- RunUMAP(object = gsm_atac, reduction = 'lsi', dims = 2:30)
gsm_atac <- FindNeighbors(object = gsm_atac, reduction = 'lsi', dims = 2:30)
gsm_atac <- FindClusters(object = gsm_atac, verbose = FALSE, algorithm = 3)
pdf("umap_only_clustering.pdf")
DimPlot(object = gsm_atac, label = TRUE) + NoLegend()
dev.off()

## Load the pre-processed scRNA-seq data
if(!is.null(rds)){
  gsm_rna <- readRDS(rds)

  ## Create a gene activity matrix
  gene.activities <- GeneActivity(gsm_atac)
  gsm_atac[['RNA']] <- CreateAssayObject(counts = gene.activities)
  gsm_atac <- NormalizeData(object = gsm_atac,assay = 'RNA',normalization.method = 'LogNormalize',scale.factor = median(gsm_atac$nCount_RNA))
  DefaultAssay(gsm_atac) <- 'RNA'

  # Classical association algorithm，CCA 
  transfer.anchors <- FindTransferAnchors(  
    reference = gsm_rna,
    query = gsm_atac,
    reduction = 'cca'
  ) 

  predicted.labels <- TransferData(
    anchorset = transfer.anchors,
    refdata = gsm_rna$celltype,
    weight.reduction = gsm_atac[['lsi']],
    dims = 2:30
  )

  gsm_atac <- AddMetaData(object = gsm_atac, metadata = predicted.labels)

  # scRNA umap plot
  pdf("scRNA_umap_org.pdf")
  plot1 <- DimPlot(object = gsm_rna,group.by = 'celltype',label = TRUE,repel = TRUE) + NoLegend() +ggtitle('scRNA-seq')
  print(plot1)
  dev.off()

  # scATAC umap plot
  pdf("scATAC_umap_with_integration.pdf")  
  plot2 <- DimPlot(object = gsm_atac,group.by = 'predicted.id',label = TRUE,repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')
  print(plot2)
  dev.off()
}
