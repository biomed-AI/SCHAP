library(Signac)
library(Seurat)
library(JASPAR2020)
library(BSgenome.Mmusculus.UCSC.mm10)
library(BSgenome.Hsapiens.UCSC.hg38)
library(patchwork)
library(tidyr)
library(dplyr)
library(ggplot2)
library(TFBSTools)


species = commandArgs(trailingOnly = TRUE)[1]
da_peaks_ident = commandArgs(trailingOnly = TRUE)[2]

print(species)
print(da_peaks_ident)

h5=list.files('../input/sample/','.h5',full.names=T,recursive=F)
csv=list.files('../input/sample/','.csv',full.names=T,recursive=F)
fragment=list.files('../input/sample/','.tsv.gz$',full.names=T,recursive=F)
celllabels=list.files('../input/','.csv',full.names=T,recursive=F)

func2 <- function(regions, sep = c(":", "-"), ...) {
  ranges.df <- data.frame(ranges = regions)
  ranges.df <- separate(
    data = ranges.df,
    col = "ranges",
    sep = paste0(sep[[1]], "|", sep[[2]]),
    into = c("chr", "start", "end")
  )
  print(ranges.df[1:5, ])
  granges <- makeGRangesFromDataFrame(df = ranges.df, ...)
  return(granges)
}

func1 <- function(
  counts,
  data,
  min.cells = 0,
  min.features = 0,
  max.cells = NULL,
  ranges = NULL,
  motifs = NULL,
  fragments = NULL,
  genome = NULL,
  annotation = NULL,
  bias = NULL,
  positionEnrichment = NULL,
  sep = c("-", "-"),
  validate.fragments = TRUE,
  verbose = TRUE,
  ...
) {
  if (missing(x = counts) && missing(x = data)) {
    stop("Must provide either 'counts' or 'data'")
  } else if (!missing(x = counts) && !missing(x = data)) {
    stop("Either 'counts' or 'data' must be missing; both cannot be provided")
  } else if (!missing(x = counts)) {
    data.use <- counts
  } else {
    data.use <- data
  }
  if (!is.null(x = ranges)) {
    if (length(x = ranges) != nrow(x = data.use)) {
      stop("Length of supplied genomic ranges does not match number
           of rows in matrix")
    }
  } else {
    ranges <- func2(regions = rownames(x = data.use), sep = sep)
  }
  if (!isDisjoint(x = ranges)) {
    warning("Overlapping ranges supplied. Ranges should be non-overlapping.")
  }
  if (!is.null(x = annotation) & !inherits(x = annotation, what = "GRanges")) {
    stop("Annotation must be a GRanges object.")
  }
  # remove low-count cells
  ncount.cell <- colSums(x = data.use > 0)
  data.use <- data.use[, ncount.cell >= min.features]

  if (ncol(x = data.use) == 0) {
    stop("No cells retained due to minimum feature cutoff supplied")
  }

  ncell.feature <- rowSums(x = data.use > 0)
  if (!is.null(x = max.cells)) {
    if (is(object = max.cells, class2 = "character")) {
      percent.cutoff <- as.numeric(
        x = gsub(pattern = "q", replacement = "", x = max.cells)
      )
      max.cells <- (percent.cutoff / 100) * ncol(x = data.use)
    }
  } else {
    max.cells <- ncol(x = data.use)
  }
  features.keep <- (ncell.feature >= min.cells) & (ncell.feature <= max.cells)
  if (sum(features.keep) == 0) {
    stop("No features retained due to minimum cell cutoff supplied")
  }
  data.use <- data.use[features.keep, ]
  ranges <- ranges[features.keep, ]
  # re-assign row names of matrix so that it's a known granges transformation
  new.rownames <- GRangesToString(grange = ranges, sep = c("-", "-"))
  rownames(x = data.use) <- new.rownames
  if (!missing(x = counts)) {
    seurat.assay <- CreateAssayObject(
      counts = data.use,
      data = data,
      min.cells = -1,
      min.features = -1 # min cell/feature filtering already done
    )
  } else {
    seurat.assay <- CreateAssayObject(
      counts = counts,
      data = data.use,
      min.cells = min.cells,
      min.features = min.features
    )
  }
  if (inherits(x = fragments, what = "list")) {
    # check each object in the list is a fragment object
    # fragment list usually supplied when doing object merge,
    # so don't validate cells here, we can assume that was done in
    # individual object creation
    obj.class <- sapply(
      X = fragments, FUN = function(x) inherits(x = x, what = "Fragment")
    )
    if (!all(obj.class)) {
      stop("All objects in fragments list must be Fragment-class objects")
    }
    frags <- lapply(
      X = fragments,
      FUN = AssignFragCellnames,
      cellnames = colnames(x = seurat.assay)
    )
    # subset to cells in the assay
    frags <- lapply(
      X = fragments,
      FUN = subset,
      cells = colnames(x = seurat.assay)
    )
   } else if (inherits(x = fragments, what = "Fragment")) {
    # single Fragment object supplied
    frags <- AssignFragCellnames(
      fragments = fragments, cellnames = colnames(x = seurat.assay)
    )
    # subset to cells in the assay
    frags <- subset(x = frags, cells = colnames(x = seurat.assay))
  } else {
    # path to fragment file supplied, create fragment object
    frags <- list()
    if (!is.null(x = fragments)) {
      if (nchar(x = fragments) > 0) {
        cells <- colnames(x = seurat.assay)
        names(x = cells) <- cells
        frags[[1]] <- CreateFragmentObject(
          path = fragments,
          cells = cells,
          validate.fragments = validate.fragments,
          verbose = verbose,
        )
      }
    }
  }

  chrom.assay <- as.ChromatinAssay(
    x = seurat.assay,
    ranges = ranges,
    seqinfo = genome,
    motifs = motifs,
    fragments = frags,
    annotation = annotation,
    bias = bias,
    positionEnrichment = positionEnrichment
  )
  return(chrom.assay)
}

expression_matrix <- Read10X_h5(filename = h5)
print(ncol(expression_matrix))

index <- rownames(expression_matrix)
index <- sub(c(":"), c("-"), index)
row.names(expression_matrix) <- index

chromatin_assay <- func1(
  counts=expression_matrix, 
  seq=c("-", "-"),
  genome="mm10",
  fragment= fragment,
  # min.cells=10,
  # min.features=200,
)

brain <- CreateSeuratObject(
  counts=chromatin_assay, 
  assay='peaks'
)

pred <- read.csv(celllabels)$Pred
for (i in 1:length(pred)){
  if (grepl('Ex', pred[i])){
    pred[i] <- 'Ex neurons'
  }
}
print(length(pred))


brain <- AddMetaData(
  object = brain,
  metadata = pred,
  col.name = "cell_type"
)

brain <- FindVariableFeatures(brain, nfeatures = 10000)
brain <- subset(x = brain, features = VariableFeatures(object = brain))
print(brain)

pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)


brain <- AddMotifs(
  object = brain,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = pfm
)

cell_type_list <- unique(pred)

Idents(brain) <- pred


#################### can be change by users#################
da_peaks <- FindMarkers(
  object = brain,
  ident.1 = da_peaks_ident, ###parament
  # ident.2 = 'Sst',
  # group.by = 'cell_type',
  only.pos = TRUE,
  test.use = 'LR',
  min.pct = 0.05,
  latent.vars = 'nCount_peaks'
)

top.da.peak <- rownames(da_peaks[da_peaks$p_val < 0.005, ])

enriched.motifs <- FindMotifs(
  object = brain,
  features = top.da.peak
)
head(enriched.motifs)

p1 <- MotifPlot(
  object = brain,
  motifs = head(rownames(enriched.motifs))
)
ggsave("Cell_MotifPlot.pdf", p1, width=10, height=6, dpi = 600)