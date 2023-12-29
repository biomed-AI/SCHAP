#!/bin/bash

set -v on
set -x # debug 
set -e # exit script with any err
# set -u # exit script with any use of empty env

# 记录环境变量
env > env.log

mkdir input

# 指定samples的路径
tar -xvf $samplestar -C input

# 获取 input 文件夹中的文件夹数量
folder_count=$(find input -mindepth 1 -maxdepth 1 -type d | wc -l)
if [ "$folder_count" -eq 1 ]; then
    # 只有一个文件夹，将其路径赋给 sample1 环境变量
    export sample1=$(realpath input/*)
    mv $sample1 $(realpath input/sample1)
    echo -e "$sample1 --> $(realpath input/sample1)" >> mapper.txt
elif [ "$folder_count" > 1 ]; then
    # 有两个文件夹，将第一个文件夹路径赋给 sample1 环境变量，以此类推
 find input -mindepth 1 -maxdepth 1 -type d >> tmp.txt
 for i in `seq $folder_count`
  do
    export folders=$(cat tmp.txt|sed -n $i\p)
    export samples=$(realpath $folders)
    mv $samples $(realpath input/samples$i)
    # 将samples文件的路径环境变量写入文件
    echo -e "$samples --> $(realpath input/samples$i)"  >> mapper.txt
  done
else
    # 没有输入，或者解压后结构不规范
    echo "Error: You should check your input file for compliance."
    exit 1
fi
rm -rf tmp.txt

# export submitr_path=/home/zhanghaokun/submit_3.R
export submitr_path=/app/common/bio-platform/SCHAP/scRNA-generic.R


# 设置R程序需要load的数据环境变量
# 配置物种
export species=$species
echo "$species"
if [[ $species = Human ]]; then
export transcriptome=/app/common/bio-platform/SCHAP/reference/refdata-gex-GRCh38-2020-A
export Immu_Rdata_path=/app/common/bio-platform/SCHAP/humanImmu.Rdata
export RNA_Rdata_path=/app/common/bio-platform/SCHAP/humanRNA.Rdata
elif [[ $species = Mouse ]]; then
export transcriptome=/app/common/bio-platform/SCHAP/reference/refdata-gex-mm10-2020-A
export Immu_Rdata_path=/app/common/bio-platform/SCHAP/mouseImmu.Rdata
export RNA_Rdata_path=/app/common/bio-platform/SCHAP/mouseRNA.Rdata
fi

#设置画图参数

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
mkdir -p result && mv mapper.txt ./result/ && cd result && export TZ=Asia/Shanghai && /home/zhanghaokun/anaconda3/envs/singler/bin/Rscript $submitr_path $Immu_Rdata_path $RNA_Rdata_path $species $pvalueCutoff $qvalueCutoff $heatmap_width $heatmap_height $bubble_width $bubble_height 2>&1 | tee ../output.log