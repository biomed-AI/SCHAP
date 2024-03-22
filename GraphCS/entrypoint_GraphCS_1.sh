#!/bin/bash

#set -v on
set -x # debug 
set -e # exit script with any err
# set -u # exit script with any use of empty env

mkdir input

export ref=$ref
export rds=$rds

if [ -z "$rds"]; then
 tar -xvf $samplestar -C input
else
 cp $rds input/
fi

folder_2_count=$(find input -mindepth 2 -maxdepth 2 -type d | wc -l)
if [ "$folder_2_count" -eq 0 ] || [ -z "$rds"]; then

  folder_count=$(find input -mindepth 1 -maxdepth 1 -type d | wc -l)

elif [ "$folder_2_count" > 0 ] || [ -z "$rds"]; then
  
  mv input/* input/folder_useless
  mv input/folder_useless/* input/
  rm input/folder_useless -rf
  folder_count=$(find input -mindepth 1 -maxdepth 1 -type d | wc -l)

else
    echo "Error: You should check your input file for compliance."
    exit 1
fi

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
  rm -rf tmp.txt
else
    
    echo "Error: You should check your input file for compliance."
    exit 1
fi


export submitr_path=/app/common/bio-platform/SCHAP/GraphCS/GraphCS_1.R


export species=$species
echo "$species"
if [[ $species = Human ]]; then
export transcriptome=/app/common/bio-platform/SCHAP/reference/refdata-gex-GRCh38-2020-A
export Immu_Rdata_path=/app/common/bio-platform/SCHAP/humanImmu.Rdata
export RNA_Rdata_path=/app/common/bio-platform/SCHAP/humanRNA.Rdata
export ref=/app/common/bio-platform/SCHAP/reference/GraphCS_human.Rdata
elif [[ $species = Mouse ]]; then
export transcriptome=/app/common/bio-platform/SCHAP/reference/refdata-gex-mm10-2020-A
export Immu_Rdata_path=/app/common/bio-platform/SCHAP/mouseImmu.Rdata
export RNA_Rdata_path=/app/common/bio-platform/SCHAP/mouseRNA.Rdata
export ref=/app/common/bio-platform/SCHAP/reference/GraphCS_mouse.Rdata
fi

mkdir -p /home/zhanghaokun/anaconda3/envs/
ln -s /app/common/bio-platform/SCHAP/singler /home/zhanghaokun/anaconda3/envs/singler

if [ -z "$ref"]; then
 echo "Using the default reference!"
else
 cp $ref input/
fi

mkdir -p result && mv mapper.txt ./result/ && cd result && export TZ=Asia/Shanghai && /home/zhanghaokun/anaconda3/envs/singler/bin/Rscript $submitr_path $Immu_Rdata_path $RNA_Rdata_path $species $ref $rds 2>&1 | tee ../output.log
