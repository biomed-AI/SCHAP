knitr::opts_chunk$set(warning = FALSE, message = FALSE) 
require(magrittr)
require(readr)
require(Matrix)
require(tidyr)
require(dplyr)
library(hdf5r)
library(cicero)
library(Signac)
library(Seurat)
library(EnsDb.Mmusculus.v79)
library(ggplot2)
library(patchwork)
library(monocle3)
library(stringr)

species = commandArgs(trailingOnly = TRUE)[1]
anno_data = commandArgs(trailingOnly = TRUE)[2]
## region format can be changed to 'gene_name'（Cd68） or a range like 'chr1-10-1000'
coveragePlotregion = commandArgs(trailingOnly = TRUE)[3]
## genome can be changed to a'chr1' or 'chr(2:19)'
connections_genome = commandArgs(trailingOnly = TRUE)[4]
#pos1 <- 75863738 ## parameter
#pos2 <- 75884283 ## parameter
pos1 = commandArgs(trailingOnly = TRUE)[5]
pos2 = commandArgs(trailingOnly = TRUE)[6]
pos1 <- as.numeric(pos1)
pos2 <- as.numeric(pos2)

print(species)
print(coveragePlotregion)
print(connections_genome)
print(pos1)
print(pos2)


h5=list.files('../input/sample/','.h5',full.names=T,recursive=F)
csv=list.files('../input/sample/','.csv',full.names=T,recursive=F)
fragment=list.files('../input/sample/','.tsv.gz$',full.names=T,recursive=F)
celllabels=list.files('../input/','.csv',full.names=T,recursive=F)

counts <- Read10X_h5(filename = h5)

barcodes <- read.csv(csv,header=T)


brain_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":","-"),
  genome = "mm10",
  fragments = fragment,
  min.cells = 1
)


brain <- CreateSeuratObject(
  counts = brain_assay,
  assay = 'peaks',
  project = 'ATAC'
)



## the path might need to be changed, to find the macs2
peaks <- CallPeaks(brain,macs2.path='/home/zhanghaokun/anaconda3/envs/mamba_install/envs/monocle/bin/macs2') 


frags <- CreateFragmentObject(fragment,
                              cells = barcodes$barcode)


indata <- FeatureMatrix(fragments = frags,
                     features = peaks,cells = barcodes$barcode)


# format peak info
temp_peakinfo <- data.frame(peaks)
temp_peakinfo <- temp_peakinfo[,c(1:3)]
names(temp_peakinfo) <- c("chr", "bp1", "bp2")

temp_peakinfo$site_name <- paste(temp_peakinfo$chr, temp_peakinfo$bp1, 
                                 temp_peakinfo$bp2, sep="-")
row.names(temp_peakinfo) <-temp_peakinfo$site_name

row.names(indata) <- temp_peakinfo$site_name
colnames(indata) <- barcodes$barcode

brain_assay <- CreateChromatinAssay(
  counts = indata,
  sep = c(":", "-"),
  genome = "mm10",
  fragments = fragment,
  min.cells = 1
)
brain <- CreateSeuratObject(
  counts = brain_assay,
  assay = 'peaks',
  project = 'ATAC')

## Attach predicted cell labels to the data
cell_labels <- read.csv(celllabels)
#cell_labels$Pred <- gsub("^Ex.*", 'Excitatory neurons', cell_labels$Pred)
cell_labels$Pred <- gsub("^Ex.*", 'Excitatory neurons', cell_labels$Pred)
Idents(brain) <- cell_labels$Pred

# extract gene annotations from EnsDb
#annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)

# change to UCSC style since the data was mapped to hg19
#seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
#genome(annotations) <- "mm10"
load(anno_data)

annotations

# add the gene information to the object
Annotation(brain) <- annotations

cell_types = unique(cell_labels$Pred)
levels(brain) <- sort(cell_types)


## region format can be changed to 'gene_name'（Cd68） or a range like 'chr1-10-1000'
pdf('coveragePlot.pdf',width=6, height=6)
Signac::CoveragePlot(
  object = brain,
  region = c(coveragePlotregion),peaks = T,
  ncol = 1
)
dev.off()


cellinfo <- read.csv(celllabels)

cell_list <- unique(cell_labels$Pred)

# format peak info
temp_peakinfo <- data.frame(peaks)
temp_peakinfo <- temp_peakinfo[,c(1:3)]
names(temp_peakinfo) <- c("chr", "bp1", "bp2")
temp_peakinfo$site_name <- paste(temp_peakinfo$chr, temp_peakinfo$bp1, 
                                 temp_peakinfo$bp2, sep="-")

row.names(temp_peakinfo) <-temp_peakinfo$site_name
row.names(indata) <- temp_peakinfo$site_name


cellinfo$Pred <- gsub("^Ex.*", 'Excitatory neurons', cellinfo$Pred)
colnames(indata) <- cellinfo$Pred

colnames(indata)
## only need to change this path to your own path and the cell type that you focus on
celltype_specific_path <- 'EX_connection/' 
if (!dir.exists(celltype_specific_path)) {
    dir.create(celltype_specific_path,recursive=T)
}

## loop through every cell type and compute the connection score

for(cell in cell_list){
  sub_cellinfo <- cellinfo[cellinfo$Pred== cell,]
  sub_cellinfo <- data.frame(sub_cellinfo[, -c(1:3)])
  names(sub_cellinfo) <- 'cells'
  rownames(sub_cellinfo) <- make.names(sub_cellinfo$cells, unique=TRUE)
  sub_cellinfo$cells <- rownames(sub_cellinfo)
  
  
  test = indata[, colnames(indata) == cell ]
  colnames(test) <- sub_cellinfo$cells
  
  # make CDS monocle3
  input_cds <-  suppressWarnings(new_cell_data_set(test,
  cell_metadata = sub_cellinfo,
  gene_metadata = temp_peakinfo))
  
  #Ensure there are no peaks included with zero reads
  #input_cds <- input_cds[Matrix::rowSums(exprs(input_cds)) != 0,] 
  
  set.seed(2023)
  input_cds <- monocle3::detect_genes(input_cds)

  input_cds <- suppressMessages(estimate_size_factors(input_cds))

  input_cds <- suppressMessages(preprocess_cds(input_cds, method = "LSI"))

  input_cds <- suppressMessages(reduce_dimension(input_cds, reduction_method = 'UMAP', 
                                preprocess_method = "LSI"))
  
  umap_coords <- suppressMessages(reducedDims(input_cds)$UMAP)
  
  cicero_cds <- suppressMessages(make_cicero_cds(input_cds, reduced_coordinates = umap_coords))
  
  mm10_chrom_sizes <- read.table('/app/common/bio-platform/SCHAP/SANGO/mm10.chrom.sizes')
  
  mm10_chrom_sizes <- mm10_chrom_sizes[which(mm10_chrom_sizes$V1 %in% str_c("chr",c(1:19,"X","Y"))),]
  
  sample_genome <- subset(mm10_chrom_sizes, V1 == c(connections_genome))
  
  conns <- suppressMessages(run_cicero(cicero_cds, genomic_coords = sample_genome))
  
  file_name <- paste(celltype_specific_path,cell,'_connections.csv',sep='')
  print(file_name)
  write.csv(conns,file_name)
}


## this file needs to be downloaded from online
gene_anno <- rtracklayer::readGFF('/app/common/bio-platform/SCHAP/SANGO/Mus_musculus.GRCm39.110.gtf.gz')
gene_anno$chromosome <- paste0("chr", gene_anno$seqid)
gene_anno$gene <- gene_anno$gene_id
gene_anno$transcript <- gene_anno$transcript_id
gene_anno$symbol <- gene_anno$gene_name

gene_anno <- gene_anno[,c('chromosome','start','end','strand','gene_biotype','gene_id','transcript_id','symbol')]
names(gene_anno)[names(gene_anno) == 'transcript_id'] <- 'transcript'

gene_anno
## need to change the cell type you look at, and regions.
## also change the output path of the plot
specific_cell <- 'EX' ## parameter
chrom <- connections_genome ## parameter
#pos1 <- 75863738 ## parameter
#pos2 <- 75884283 ## parameter

##############draw connection plot######################

for(cell in cell_list){
  fc <- FoldChange(brain, ident.1 = cell,
                 ident.2 = cell_types[!cell_types %in% cell])

  fc <- fc[order(fc$avg_log2FC, decreasing = TRUE), ]
  
  ## Cell-type specific peaks are the ones above a threshold
  top_peaks <- fc[fc$avg_log2FC>=3,]
  
  file_name <- paste(specific_cell,'_connection/',cell,'_connections.csv',sep='')
  conns <- read.csv(file_name)
  conns
  temp <- c()
  for(i in rownames(top_peaks)){
    chr = strsplit(i,'-')[[1]][1]
    int1 = strtoi(strsplit(i,'-')[[1]][2])
    int2 = strtoi(strsplit(i,'-')[[1]][3])
    if((chr == chrom) & (int1 >= pos1) & (int2 <= pos2)){
      temp <- append(temp,i)
    }
  }
  
  conns$peak1_color <- FALSE
  for(con in temp){
    conns$peak1_color[conns$Peak1 == con] <- TRUE
  }
  cell <- gsub(" ", "_", cell)
  if(sum(conns$peak1_color)==0){
    ## the plot output that needs to be changed
    dirname <-paste('Cicero_connectionPlots/',chrom,'_',specific_cell,'_connections/',sep='')
    if (!dir.exists(dirname)) {
    dir.create(dirname,recursive=T)
    }
    plot_name <- paste('Cicero_connectionPlots/',chrom,'_',specific_cell,'_connections/',cell,'_connection.pdf',sep='')
    pdf(plot_name, width=6, height=6)
  
    plot_connections(conns, chrom, pos1, pos2,
                     peak_color = '#ff2606',
                   connection_width = .5, connection_ymax = 1,
                   collapseTranscripts = "longest")
    
    dev.off()
  }else{
    dirname <-paste('Cicero_connectionPlots/',chrom,'_',specific_cell,'_connections/',sep='')
    if (!dir.exists(dirname)) {
    dir.create(dirname,recursive=T)
    }
    plot_name <- paste('Cicero_connectionPlots/',chrom,'_',specific_cell,'_connections/',cell,'_connection.pdf',sep='')
    pdf(plot_name, width=6, height=6)
  
    plot_connections(conns, chrom, pos1, pos2,
                     peak_color = 'peak1_color',
                   connection_width = .5, connection_ymax = 1,
                   collapseTranscripts = "longest")
    
    dev.off()
  }
}

## this function is to write the peaks (cell-type specific and background ones)
## to bed files;
to_txt <- function(peaks,filename,bg=False){
  #split_strings <- strsplit(rownames(peaks), "-")
  split_strings <- strsplit(rownames(peaks), "_")
  peak_matrix <- do.call(rbind, split_strings)
  
  # Create a data frame with the proper column names for a BED file
  bed_df <- data.frame(
    Chromosome = peak_matrix[, 1],
    Start = as.integer(peak_matrix[, 2]),
    End = as.integer(peak_matrix[, 3])
  )

  # Define the directory based on the bg flag
  dir <- ifelse(bg, 'bg_peaks_txt/', 'cell_peaks_txt/')

  ## replace the first input to your own directory
  if (!dir.exists(dir)) {
    dir.create(dir,recursive=T)
  }
  #output <- paste('./',dir,filename,sep='')
  output <- paste0('./',dir,filename)
  print(output)
  write.table(bed_df,output,sep='\t',row.names = F,col.names = F,quote = F)
}
## Code above is to plot coverage plot, code below starts to find cell-type specific/background peaks
## this is to prepare snpsea input
## Find differentially accessible peaks between clusters (cells)

DefaultAssay(brain) <- 'peaks'

## I ran the function by hand one at a time...
## Namely, I ran the code block above with different cell_groups to redefine top_peaks, background_peaks, and then the to_txt function
## This is kinda annoying, so maybe just write a for loop to do this...but I am too lazy to even write a for loop...

for(cell in cell_list){
  ## For example find cell-type specific peaks for Excitatory neurons (cell_group1)
  fc <- FoldChange(brain, ident.1 = cell,
                   ident.2 = cell_types[!cell_types %in% cell])
  
  fc <- fc[order(fc$avg_log2FC, decreasing = TRUE), ]
  ## Cell-type specific peaks are the ones above a threshold
  top_peaks <- fc[fc$avg_log2FC>=3,]
  rownames(top_peaks) <- gsub(x = rownames(top_peaks), pattern = "-", replacement = "_")
  ## Background peaks are the ones below a threshold
  background_peaks <-fc[fc$avg_log2FC<=-3,]
  rownames(background_peaks) <- gsub(x = rownames(background_peaks), pattern = "-", replacement = "_")
  cell <- gsub(" ", "_", cell)
  
  cell_specific_name <- paste('cell_peaks_',cell,'.txt',sep='')
  bg_name <- paste('bg_peaks_',cell,'.txt',sep='')
  to_txt(top_peaks,filename = cell_specific_name,bg=F)
  to_txt(background_peaks,filename=bg_name,bg=T)
}
