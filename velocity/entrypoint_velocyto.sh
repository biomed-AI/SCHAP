#!/bin/bash

#set -v on
set -x # debug 
set -e # exit script with any err
# set -u # exit script with any use of empty env

env > env.log

mkdir input


cp $loom input/

export OMP_NUM_THREADS=1


export submitr_path=/app/common/bio-platform/SCHAP/velocyto/velocity.R

export http_proxy="http://10.20.18.21:3128" && export https_proxy="http://10.20.18.21:3128"

mkdir -p /home/zhanghaokun/anaconda3/envs/mamba_install/envs/
ln -s /app/common/bio-platform/SCHAP/R4velocyto /home/zhanghaokun/anaconda3/envs/mamba_install/envs/R4velocyto

mkdir -p result && cd result && export TZ=Asia/Shanghai && /home/zhanghaokun/anaconda3/envs/mamba_install/envs/R4velocyto/bin/Rscript $submitr_path 2>&1 | tee ../output.log
