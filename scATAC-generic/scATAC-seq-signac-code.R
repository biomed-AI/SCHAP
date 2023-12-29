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

print(h5)
print(csv)
print(tsv)

counts <- Read10X_h5(filename = h5)

metadata <- read.csv(
  file = csv,
  header = TRUE,
  row.names = 1
)

chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = genomename,
  fragments = tsv,
  min.cells = 10,
  min.features = 200
)

pbmc <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)

pbmc

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


pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
pbmc <- RunSVD(pbmc)

pdf("correlation_depth_components.pdf")
DepthCor(pbmc)
dev.off()


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
