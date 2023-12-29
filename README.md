
# SCHAP: Single Cell High-performance Analysis Platform



## Overview
SCHAP is a Single-Cell High-performance Analysis Platform primarily designed for analyzing scRNA-seq and scATAC-seq data. It provides users with the option of one-click or step-by-step data analysis control, as well as user-friendly visualization of the results on the website. The platform is supported by the Starlight cloud platform on Tianhe-2, supporting CPU and GPU heterogeneous computing scheduling. It can allocate different types of computing resources for computing tasks based on CPU-intensive and memory-intensive computing requirements.


## Accessing SCHAP
The website link is available at https://bio-web1.nscc-gz.cn/database/SCHAP.
The application is free and open to all users with no login requirement. The analyzedresults are easily downloaded via URL. 


## Requirements
You'll need to install the following packages in order to run the codes in your local environment.
- cellranger7.1.0
- cellranger-atac2.1.0
- SeuratV4
- Harmony
- clusterprofile
- CellChat
- velocyto
- Signac
- Cireco
- ADClust
- SingleR
- GraphCS


##Installation
conda ....
conda ....


## Usage
Due to the upload size limits, you can utilize our code in your local environment.

- scRNA-generic
The workflow includes Seurat, Harmony, SingleR, and CellChat. We Seurat and Harmony were utilized for data integration and clustering, SingleR for cell annotation, and CellChat for calculating cell communication.
```  
sh entrypoint SCRNA generic.sh
```

- scATAC-generic
The workflow includes Seurat and Signac. Seurat and Signac were utilized for data integration and clustering, with Signac used for annotation and computing QC metrics.
```  
sh entrypoint SCRNA generic.sh
```

- ADClust
The workflow primarily involves the use of Seurat, Harmony, SingleR, and CellChat. Seurat and Harmony were utilized for data integration, Seurat for clustering, and SingleR for cell annotation. Finally, CellChat was employed to calculate cell communication.
```
sh entrypoint ADClust_1.sh
sh entrypoint ADClust_2.sh
```

- GraphCS

The workflow includes Seurat, Harmony, SingleR, and CellChat. We Seurat and Harmony were utilized for data integration and clustering, GraphCS for cell annotation, and CellChat for calculating cell communication.

```
sh entrypoint_GraphCS 1.sh
sh entrypoint_GraphCS 1.sh
```


## Tutorial
For the step-by-step tutorial, please refer to: https://bio-web1-pre.nscc-gz.cn/database/SCHAP/documentation#1
