#!/bin/bash

set -v on
set -x # debug 
set -e # exit script with any err
# set -u # exit script with any use of empty env

# 记录环境变量
env > env.log

mkdir -p input

# 指定samples的路径
tar -xvf $samplestar -C input

# 改名

cd input && mv * sample && cd ../

cp $celltype_pred input/

# 提交脚本的路径
export submitr_path=/app/common/bio-platform/SCHAP/SANGO/MotifPlot.R

# 设置R程序需要load的数据环境变量
# 配置物种
export species=$species
echo "$species"

#if [[ $species = Mouse ]]; then
#export anno_data=/app/common/bio-platform/SCHAP/SANGO/mmu_annotation.Rdata
#fi

## (Default 'Microglia')
export da_peaks_ident=$da_peaks_ident

# 设置联网proxy
export http_proxy="http://10.20.18.21:3128" && export https_proxy="http://10.20.18.21:3128"

# 链接程序到指定位置
############链接指定库#############

mkdir -p /home/zhanghaokun/anaconda3/envs/mamba_install/envs/
ln -s /app/common/bio-platform/SCHAP/motif  /home/zhanghaokun/anaconda3/envs/mamba_install/envs/motif

#export LD_LIBRARY_PATH=/home/zhanghaokun/anaconda3/envs/mamba_install/envs/motif/lib:$LD_LIBRARY_PATH
#echo "$LD_LIBRARY_PATH"

export PATH=/home/zhanghaokun/anaconda3/envs/mamba_install/envs/motif/bin:$PATH

# 运行
mkdir -p result && cd result && export TZ=Asia/Shanghai && /home/zhanghaokun/anaconda3/envs/mamba_install/envs/motif/bin/Rscript $submitr_path $species $da_peaks_ident 2>&1 | tee ../output.log