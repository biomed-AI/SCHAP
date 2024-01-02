#!/bin/bash

#set -v on
set -x # debug 
set -e # exit script with any err
# set -u # exit script with any use of empty env

env > env.log

mkdir -p input

tar -xvf $samplestar -C input


cd input && mv * sample && cd ../

cp $celltype_pred input/

export submitr_path=/app/common/bio-platform/SCHAP/SANGO/MouseBrain_DownStreamAnalysis.R

export submitr_path2=/app/common/bio-platform/SCHAP/SANGO/liftOver-snpsea.sh

export species=$species
echo "$species"

if [[ $species = Mouse ]]; then
export anno_data=/app/common/bio-platform/SCHAP/SANGO/mmu_annotation.Rdata
fi

## region format can be changed to 'gene_name'（Cd68） or a range like 'chr1-10-1000'
export coveragePlotregion=$coveragePlot_region
## genome can be changed to a'chr1' or 'chr(2:19)'
export connections_genome=$connections_genome
#pos1 <- 75863738 ## parameter
#pos2 <- 75884283 ## parameter
export pos1=$connections_pos1
export pos2=$connections_pos2


export http_proxy="http://10.20.18.21:3128" && export https_proxy="http://10.20.18.21:3128"

mkdir -p /home/zhanghaokun/anaconda3/envs/mamba_install/envs/
ln -s /app/common/bio-platform/SCHAP/monocle  /home/zhanghaokun/anaconda3/envs/mamba_install/envs/monocle
ln -s /app/common/bio-platform/SCHAP/snpsea  /home/zhanghaokun/anaconda3/envs/mamba_install/envs/snpsea
export PATH=/home/zhanghaokun/anaconda3/envs/mamba_install/envs/snpsea/bin:$PATH

mkdir -p result && cd result && export TZ=Asia/Shanghai && /home/zhanghaokun/anaconda3/envs/mamba_install/envs/monocle/bin/Rscript $submitr_path $species $anno_data $coveragePlotregion $connections_genome $pos1 $pos2 2>&1 | tee ../output.log && $submitr_path2 2>&1 | tee ../output2.log
