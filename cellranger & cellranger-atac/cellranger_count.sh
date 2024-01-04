#!/bin/bash
export transcriptome=$1
export fastq=$(realpath cellranger_input)
export PATH=/app/common/bio-platform/SCHAP/cellranger-7.1.0:$PATH

cd cellranger_input
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
    echo "Error: You should check your input file for compliance."
    exit 1
fi
rm -rf samples_fullname.txt

mkdir ../result && mv *.txt ../result && cd ../result

for i in `seq $(cat samples_name.txt|wc -l)`
do
    mkdir sample$i
    cellranger count --id=run_count_sample$i \
		--fastqs=$fastq \
		--sample=sample$i \
		--transcriptome=$transcriptome \
		--nosecondary
    mv run_count_sample$i/outs/filtered_feature_bc_matrix/* sample$i
done
rm -rf ../cellranger_input/
rm -rf samples_name.txt run_count_*
tar -zcvf results.tar sample*
rm -rf sample*
