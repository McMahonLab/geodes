#!/bin/bash

mkdir input
mv $1.fna input/

checkm lineage_wf input output

grep ">" input/$1.fna > contigs.txt
sed -i 's/>//g' contigs.txt
while read line; do
	echo $line $1;
	done < contigs.txt > $1-contigs.txt

tax=$(awk -F', ' '{for(i=1;i<=NF;i++){if ($i ~ /marker lineage/){print $i}}}' output/storage/bin_stats_ext.tsv | awk -F': ' '{print $2}')
length=$(awk -F', ' '{for(i=1;i<=NF;i++){if ($i ~ /Genome size/){print $i}}}' output/storage/bin_stats_ext.tsv | awk -F': ' '{print $2}')
complete=$(awk -F', ' '{for(i=1;i<=NF;i++){if ($i ~ /Completeness/){print $i}}}' output/storage/bin_stats_ext.tsv | awk -F': ' '{print $2}')
contamination=$(awk -F', ' '{for(i=1;i<=NF;i++){if ($i ~ /Contamination/){print $i}}}' output/storage/bin_stats_ext.tsv | awk -F': ' '{print $2}')

echo $1 $tax $length $complete $contamination > $1-checkm.txt

rm -r input/
rm -r output/
