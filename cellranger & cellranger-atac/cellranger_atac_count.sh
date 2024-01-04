#!/bin/bash
export transcriptome=$1
export fastq=$(realpath cellranger_atac_input)
export PATH=/app/common/bio-platform/SCHAP/cellranger-atac-2.1.0:$PATH
# 获取 cellranger_input 文件夹中的文件夹数量
cd cellranger_atac_input
file_count=$(ls *.fastq.gz|cut -f 1 -d "_"|uniq|wc -l)
echo "`ls *.fastq.gz`" >>samples_fullname.txt
echo "`ls *.fastq.gz|cut -f 1 -d "_"|uniq`" >>samples_name.txt
if [ "$file_count" -eq 1 ]; then
  for i in `cat samples_name.txt`
  do
    for n in `cat samples_fullname.txt`
      do
      echo -e "$n --> $(realpath sample1_$(echo "`(echo $n)|cut -f 2- -d "_"`"))" >> mapper.txt
      rename $i sample1 $n
      done
  done
elif [ "$file_count" > 1 ]; then
    # 有两个文件夹，将第一个文件夹路径赋给 sample1 环境变量，以此类推
  for i in `seq $file_count`
  do
    for n in `cat samples_fullname.txt`
    do
      export samples=$(cat samples_name.txt|sed -n $i\p)
      export original=$(echo $n|grep "$samples")
      [ ! -z "$original" ] && echo -e "$original --> $(realpath sample$i\_$(echo "`(echo $original)|cut -f 2- -d "_"`"))" >> mapper.txt
      [ ! -z "$original" ] && rename $samples sample$i $original
    done
  done
else
    # 没有输入，或者解压后结构不规范
    echo "Error: You should check your input file for compliance."
    exit 1
fi
rm -rf samples_fullname.txt

mkdir ../result && mv *.txt ../result && cd ../result

for i in `seq $(cat samples_name.txt|wc -l)`
do
    mkdir sample$i
    cellranger-atac count --id=atac_count_sample$i \
		--fastqs=$fastq \
		--sample=sample$i \
		--reference=$transcriptome

    mv atac_count_sample$i/outs/filtered_peak_bc_matrix.h5 sample$i
    #zcat atac_count_sample$i/outs/fragments.tsv.gz |  awk '!/^#/' | gzip > atac_count_sample$i/outs/fragments2.tsv.gz
    #cp atac_count_sample$i/outs/fragments2.tsv.gz sample$i/fragments.tsv.gz
    mv atac_count_sample$i/outs/fragments.tsv.gz sample$i
    mv atac_count_sample$i/outs/fragments.tsv.gz.tbi sample$i
    mv atac_count_sample$i/outs/singlecell.csv sample$i
done
rm -rf ../cellranger_atac_input
rm -rf samples_name.txt atac_count_*
tar -zcvf results.tar sample*
rm -rf sample*
