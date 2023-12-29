#!/bin/bash

set -v on
set -x # debug 
set -e # exit script with any err
# set -u # exit script with any use of empty env

# 记录环境变量
env > env.log

mkdir input

# 指定sample的路径
cp $loom input/

#设定单核运行
export OMP_NUM_THREADS=1

#指定执行代码
export submitr_path=/app/common/bio-platform/SCHAP/velocyto/velocity.R

# 设置联网proxy
export http_proxy="http://10.20.18.21:3128" && export https_proxy="http://10.20.18.21:3128"

# 链接程序到指定位置
mkdir -p /home/zhanghaokun/anaconda3/envs/mamba_install/envs/
ln -s /app/common/bio-platform/SCHAP/R4velocyto /home/zhanghaokun/anaconda3/envs/mamba_install/envs/R4velocyto

# 运行
mkdir -p result && cd result && export TZ=Asia/Shanghai && /home/zhanghaokun/anaconda3/envs/mamba_install/envs/R4velocyto/bin/Rscript $submitr_path 2>&1 | tee ../output.log