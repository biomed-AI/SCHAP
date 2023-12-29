#!/bin/bash
cd bg_peaks_txt

echo "`ls *.txt|cut -f 1 -d "."`" >samples_name.txt
mkdir ../bg_peaks_bed
for i in `cat samples_name.txt`
do
     echo $i
     /app/common/bio-platform/SCHAP/SANGO/liftOver $i.txt /app/common/bio-platform/SCHAP/SANGO/mm9ToMm10.over.chain $i.bed unmapped
     rm -f unmapped
     mv $i.bed ../bg_peaks_bed
done
rm samples_name.txt -f

cd ../cell_peaks_txt

echo "`ls *.txt|cut -f 1 -d "."`" >>samples_name.txt
mkdir ../cell_peaks_bed
for i in `cat samples_name.txt`
do
    /app/common/bio-platform/SCHAP/SANGO/liftOver $i.txt /app/common/bio-platform/SCHAP/SANGO/mm9ToMm10.over.chain $i.bed unmapped
    rm -f unmapped
    mv $i.bed ../cell_peaks_bed
done

rm samples_name.txt -f

## Run SNPsea ##

cd ../ && mkdir bedfile && cp cell_peaks_bed/* bedfile && cp bg_peaks_bed/* bedfile

cd bedfile
mkdir ../SNPsea_result

for i in `ls *.bed`
do
export bed=$i
export name=$(echo $bed |cut -f 1 -d ".")
awk '{print $0, "rs" $NF}' $bed >tmp.bed
awk '{$3=""; print $0}' tmp.bed > $bed
sed -i '1i\CHR POS SNP' $bed
sed -i 's/\s\+/\t/g' $bed
rm tmp.bed -f

process_do(){
options=(
    --snps              $bed
    --gene-matrix       /app/common/bio-platform/SCHAP/SANGO/snpsea_dataset/GeneAtlas2004.gct.gz
    --gene-intervals    /app/common/bio-platform/SCHAP/SANGO/snpsea_dataset/NCBIgenes2013.bed.gz
    --snp-intervals     /app/common/bio-platform/SCHAP/SANGO/snpsea_dataset/TGP2011.bed.gz
    --null-snps         /app/common/bio-platform/SCHAP/SANGO/snpsea_dataset/Lango2010.txt.gz
    --out               $name\_cell_out
    --slop              10e3
    --threads           2
    --null-snpsets      0
    --min-observations  100
    --max-iterations    1e7
)

# Run SNPsea.
/app/common/bio-platform/SCHAP/SANGO/snpsea-1.0.3/bin/snpsea-linux64 ${options[*]}

# Create a horizontal bar plot of condition p-values.
/app/common/bio-platform/SCHAP/SANGO/snpsea-1.0.3/bin/snpsea-barplot $name\_cell_out

cp $name\_cell_out/condition_pvalues_barplot.pdf ../SNPsea_result/$name\_barplot.pdf

mv $name\_cell_out ../SNPsea_result/

}

process_do &

done

wait

echo "SNPsea barplot have completed"