#
# SCHAP: Single Cell High-performance Analysis Platform

<div align=center>
<img src="https://github.com/biomed-AI/SCHAP/assets/110893478/e17874e6-6311-4bd7-b688-7687c797d18f"/>
</div>


## Overview
SCHAP is a Single-Cell High-performance Analysis Platform primarily designed for analyzing scRNA-seq and scATAC-seq data. It provides users with the option of one-click or step-by-step data analysis control, as well as user-friendly visualization of the results on the website. The platform is supported by the Starlight cloud platform on Tianhe-2, supporting CPU and GPU heterogeneous computing scheduling. It can allocate different types of computing resources for computing tasks based on CPU-intensive and memory-intensive computing requirements.


## Accessing SCHAP
The website link is available at https://bio-web1.nscc-gz.cn/database/SCHAP.
The application is free and open to all users with no login requirement. The analyzed results are easily downloaded via URL. 


## Requirements
You'll need to install the following packages in order to run the codes in your local environment.
#
- Centos 7.6.1810
#
- cellranger 7.1.0
- cellranger-atac 2.1.0
- Seurat v4
- Harmony
- clusterprofile
- CellChat
- velocyto
- Signac
- Monocle 3
- Cireco
- SNPsea
- liftover
- [ADClust](https://github.com/biomed-AI/ADClust)
- SingleR
- [GraphCS](https://github.com/biomed-AI/GraphCS)
- [SANGO](https://github.com/biomed-AI/SANGO)


## Installation
### Algorithm Installation
- To install the ADClust, please follow the guidelines in https://github.com/biomed-AI/ADClust
- To install the GraphCS, please follow the guidelines in https://github.com/biomed-AI/GraphCS
- To install the SANGO, please follow the guidelines in https://github.com/biomed-AI/SANGO

### Environment Installation
- To run scRNA-generic / preprocess & visualization for ADClust and GraphCS

```bash
conda create -n singler --file spec-scrna-list.txt
```
- To run the visualization for SANGO

```bash
conda create -n monocle --file spec-monocle-list.txt
conda create -n motif --file spec-motif-list.txt
conda create -n snpsea --file spec-snpsea-list.txt
```
- To run scATAC-generic

```bash
conda create -n atacenv --file spec-scatac-list.txt
```
- To run velocyto trajectory analysis

```bash
conda create -n R4velocyto --file spec-R4velocyto-list.txt
```

## Usage
Due to the upload size limits for website, you can also utilize our code in your local environment.

- cellranger count

The code to process the fastq files of scRNA-seq data.
```bash
bash entrypoint_cellranger.sh
```

- cellranger-atac count

The code to process the fastq files of scATAC-seq data.
```bash
bash entrypoint_cellranger_atac.sh
```

- scRNA-generic

The workflow includes Seurat, Harmony, SingleR, and CellChat. We used Seurat and Harmony were utilized for data integration and clustering, SingleR for cell annotation, and CellChat for calculating cell communication.
```bash
bash entrypoint_scRNA_generic.sh
```

- scATAC-generic

The workflow includes Seurat and Signac. Seurat and Signac were utilized for data integration and clustering, with Signac used for annotation and computing QC metrics.
```bash
bash entrypoint_scatac.sh
```

- ADClust

The workflow primarily involves the use of Seurat, Harmony, ADClust, SingleR, and CellChat. We used Seurat and Harmony were utilized for data integration, ADClust for clustering, and SingleR for cell annotation. Finally, CellChat was employed to calculate cell communication.

Preprocess the raw data to satisfy the input of ADClust
```bash
bash entrypoint_ADClust_1.sh
```
Visualization for the result generated by ADClust
```bash
bash entrypoint_ADClust_2.sh
```

- GraphCS

The workflow includes Seurat, Harmony, GraphCS, and CellChat. We used Seurat and Harmony were utilized for data integration and clustering, GraphCS for cell annotation, and CellChat for calculating cell communication.

Preprocess the raw data to satisfy the input of GraphCS
```bash
bash entrypoint_GraphCS_1.sh
```
Visualization for the result generated by GraphCS
```bash
bash entrypoint_GraphCS_2.sh
```

- SANGO

This is a workflow for cell annotation of scATAC-seq data using the SANGO tool, and we also perform downstream analysis of the annotated data, including enrichment analysis, chromatin proximity analysis, and chromatin interaction analysis.

The datasets used for the visualization of SANGO are hosted by the following service: [Google Drive](https://drive.google.com/drive/folders/1N478o7W4rvBSk6H48QfH-fKAPxiO6Iyo?usp=sharing)

Visualization for the result generated by SANGO
```bash
bash entrypoint_SANGO.sh
```
Draw the motifplot for the result generated by SANGO
```bash
bash entrypoint_motifplot.sh
```



## Tutorial
For the website tutorial, please refer to: https://bio-web1.nscc-gz.cn/database/SCHAP/documentation#1
