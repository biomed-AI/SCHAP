#!/bin/bash

#set -v on
set -x # debug 
set -e # exit script with any err
# set -u # exit script with any use of empty env

mkdir input

cp $Rdata input/
cp $celltype_pred input/
export Rdata=$(realpath input/*.Rdata)
export celltype_pred=$(realpath input/*.csv)

export submitr_path=/app/common/bio-platform/SCHAP/GraphCS/GraphCS_2.R

export pvalueCutoff=$pvalueCutoff
export qvalueCutoff=$qvalueCutoff
export heatmap_width=$heatmap_width
export heatmap_height=$heatmap_height
export bubble_width=$bubble_width
export bubble_height=$bubble_height

export http_proxy="http://10.20.18.21:3128" && export https_proxy="http://10.20.18.21:3128"

mkdir -p /home/zhanghaokun/anaconda3/envs/
ln -s /app/common/bio-platform/SCHAP/singler /home/zhanghaokun/anaconda3/envs/singler

mkdir -p result && cd result && export TZ=Asia/Shanghai && /home/zhanghaokun/anaconda3/envs/singler/bin/Rscript $submitr_path $Rdata $celltype_pred $pvalueCutoff $qvalueCutoff $heatmap_width $heatmap_height $bubble_width $bubble_height 2>&1 | tee ../output.log
