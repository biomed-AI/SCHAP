#!/bin/bash

set -v on
set -x # debug 
set -e # exit script with any err
# set -u # exit script with any use of empty env


env > env.log

mkdir input


tar -xvf $samplestar -C input


folder_count=$(find input -mindepth 1 -maxdepth 1 -type d | wc -l)
if [ "$folder_count" -eq 1 ]; then
   
    export sample1=$(realpath input/*)
    mv $sample1 $(realpath input/sample1)
    echo -e "$sample1 --> $(realpath input/sample1)" >> mapper.txt
elif [ "$folder_count" > 1 ]; then
    
 find input -mindepth 1 -maxdepth 1 -type d >> tmp.txt
 for i in `seq $folder_count`
  do
    export folders=$(cat tmp.txt|sed -n $i\p)
    export samples=$(realpath $folders)
    mv $samples $(realpath input/samples$i)
    
    echo -e "$samples --> $(realpath input/samples$i)"  >> mapper.txt
  done
else

    echo "Error: You should check your input file for compliance."
    exit 1
fi
rm -rf tmp.txt


export submitr_path=/app/common/bio-platform/SCHAP/scRNA-generic.R



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



export pvalueCutoff=$pvalueCutoff
export qvalueCutoff=$qvalueCutoff
export heatmap_width=$heatmap_width
export heatmap_height=$heatmap_height
export bubble_width=$bubble_width
export bubble_height=$bubble_height


export http_proxy="http://10.20.18.21:3128" && export https_proxy="http://10.20.18.21:3128"


mkdir -p /home/zhanghaokun/anaconda3/envs/
ln -s /app/common/bio-platform/SCHAP/singler /home/zhanghaokun/anaconda3/envs/singler


mkdir -p result && mv mapper.txt ./result/ && cd result && export TZ=Asia/Shanghai && /home/zhanghaokun/anaconda3/envs/singler/bin/Rscript $submitr_path $Immu_Rdata_path $RNA_Rdata_path $species $pvalueCutoff $qvalueCutoff $heatmap_width $heatmap_height $bubble_width $bubble_height 2>&1 | tee ../output.log
