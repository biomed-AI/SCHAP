#!/bin/bash

set -v on
set -x # debug 
set -e # exit script with any err
# set -u # exit script with any use of empty env

env > env.log

mkdir -p input

tar -xvf $samplestar -C input


cd input && mv * sample && cd ../

export submitr_path=/app/common/bio-platform/SCHAP/scatac/scATAC-seq-signac-code.R


export species=$species
echo "$species"

if [[ $species = Human ]]; then
export anno_data=/app/common/bio-platform/SCHAP/scatac/hsa_annotation.Rdata
elif [[ $species = Mouse ]]; then
export anno_data=/app/common/bio-platform/SCHAP/scatac/mmu_annotation.Rdata
fi


export FragmentHistogram_region=$FragmentHistogram_region

export http_proxy="http://10.20.18.21:3128" && export https_proxy="http://10.20.18.21:3128"

mkdir -p /home/zhanghaokun/anaconda3/envs/mamba_install/envs/
ln -s /app/common/bio-platform/SCHAP/atacenv  /home/zhanghaokun/anaconda3/envs/mamba_install/envs/atacenv

export PATH=/home/zhanghaokun/anaconda3/envs/mamba_install/envs/atacenv/bin:$PATH

# 运行
mkdir -p result && cd result && export TZ=Asia/Shanghai && Rscript $submitr_path $species $anno_data $FragmentHistogram_region 2>&1 | tee ../output.log
