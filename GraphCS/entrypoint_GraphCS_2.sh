#!/bin/bash

set -v on
set -x # debug 
set -e # exit script with any err
# set -u # exit script with any use of empty env

# 记录环境变量
env > env.log

mkdir input

# 指定samples的路径
cp $Rdata input/
cp $celltype_pred input/
export Rdata=$(realpath input/*.Rdata)
export celltype_pred=$(realpath input/*.csv)

# export submitr_path=/home/zhanghaokun/GraphCS_2.R
export submitr_path=/app/common/bio-platform/SCHAP/GraphCS/GraphCS_2.R

# 设置画图参数


export pvalueCutoff=$pvalueCutoff
export qvalueCutoff=$qvalueCutoff
export heatmap_width=$heatmap_width
export heatmap_height=$heatmap_height
export bubble_width=$bubble_width
export bubble_height=$bubble_height


# 设置联网proxy
export http_proxy="http://10.20.18.21:3128" && export https_proxy="http://10.20.18.21:3128"

# 链接程序到指定位置
mkdir -p /home/zhanghaokun/anaconda3/envs/
ln -s /app/common/bio-platform/SCHAP/singler /home/zhanghaokun/anaconda3/envs/singler

# 运行
mkdir -p result && cd result && export TZ=Asia/Shanghai && /home/zhanghaokun/anaconda3/envs/singler/bin/Rscript $submitr_path $Rdata $celltype_pred $pvalueCutoff $qvalueCutoff $heatmap_width $heatmap_height $bubble_width $bubble_height 2>&1 | tee ../output.log