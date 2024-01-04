#!/bin/bash

#set -v on
set -x # debug 
set -e # exit script with any err
# set -u # exit script with any use of empty env

mkdir -p input

tar -xvf $samplestar -C input

cd input && mv * sample && cd ../

cp $celltype_pred input/

export submitr_path=/app/common/bio-platform/SCHAP/SANGO/MotifPlot.R

export species=$species
echo "$species"


## (Default 'Microglia')
export da_peaks_ident=$da_peaks_ident

mkdir -p /home/zhanghaokun/anaconda3/envs/mamba_install/envs/
ln -s /app/common/bio-platform/SCHAP/motif  /home/zhanghaokun/anaconda3/envs/mamba_install/envs/motif

export PATH=/home/zhanghaokun/anaconda3/envs/mamba_install/envs/motif/bin:$PATH

mkdir -p result && cd result && export TZ=Asia/Shanghai && /home/zhanghaokun/anaconda3/envs/mamba_install/envs/motif/bin/Rscript $submitr_path $species $da_peaks_ident 2>&1 | tee ../output.log
