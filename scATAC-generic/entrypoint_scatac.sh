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

# 提交脚本的路径
export submitr_path=/app/common/bio-platform/SCHAP/scatac/scATAC-seq-signac-code.R

# 设置R程序需要load的数据环境变量
# 配置物种
export species=$species
echo "$species"

if [[ $species = Human ]]; then
export anno_data=/app/common/bio-platform/SCHAP/scatac/hsa_annotation.Rdata
elif [[ $species = Mouse ]]; then
export anno_data=/app/common/bio-platform/SCHAP/scatac/mmu_annotation.Rdata
fi

##配置参数

export FragmentHistogram_region=$FragmentHistogram_region


# 设置联网proxy
export http_proxy="http://10.20.18.21:3128" && export https_proxy="http://10.20.18.21:3128"

# 链接程序到指定位置
mkdir -p /home/zhanghaokun/anaconda3/envs/mamba_install/envs/
ln -s /app/common/bio-platform/SCHAP/atacenv  /home/zhanghaokun/anaconda3/envs/mamba_install/envs/atacenv

export PATH=/home/zhanghaokun/anaconda3/envs/mamba_install/envs/atacenv/bin:$PATH

# 运行
mkdir -p result && cd result && export TZ=Asia/Shanghai && Rscript $submitr_path $species $anno_data $FragmentHistogram_region 2>&1 | tee ../output.log