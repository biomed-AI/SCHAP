#!/bin/bash

#set -v on
set -x # debug 
set -e # exit script with any err
# set -u # exit script with any use of empty env

mkdir cellranger_atac_input
chmod 755 cellranger_atac_input
#export samplestar=$(ls /GPUFS/sysu_ydyang_5/bio-platform/cellranger-samples/*.tar)

tar -xvf $samplestar -C cellranger_atac_input

export species=$species
echo "$species"
# export submitr_path=/home/zhanghaokun/cellranger_count.sh
export submit_path=/app/common/bio-platform/SCHAP/cellranger_atac_count.sh
export PATH=/app/common/bio-platform/SCHAP/cellranger-atac-2.1.0:$PATH

if [[ $species = Human ]]; then
export transcriptome=/app/common/bio-platform/SCHAP/reference/refdata-cellranger-arc-GRCh38-2020-A-2.0.0
elif [[ $species = Mouse ]]; then
export transcriptome=/app/common/bio-platform/SCHAP/reference/refdata-cellranger-arc-mm10-2020-A-2.0.0
fi

# Run the code
bash $submit_path $transcriptome 2>&1 | tee output.log
