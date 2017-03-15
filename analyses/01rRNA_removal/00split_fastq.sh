#!/bin/bash
tar -xvzf BBMap_36.99.tar.gz

cp /mnt/gluster/amlinz/GEODES_metaT/$1.filter-MTF.fastq.gz .
gzip -d $1.filter-MTF.fastq.gz

sed -i '/^$/d' $1.filter-MTF.fastq

maxreads=$((`wc -l < $1.filter-MTF.fastq` / 8 - 1))
startpoints=$(seq 0 500000 $maxreads)

for num in $startpoints;
  do endpoint=$(($num + 499999));
  bbmap/getreads.sh in=$1.filter-MTF.fastq id=$num-$endpoint out=$1_filter-MTF_$endpoint.fastq overwrite=T;
  done

rm $1.filter-MTF.fastq 
gzip $1*
cp $1* /mnt/gluster/amlinz/GEODES_metaT_split
rm $1*
rm BBMap_36.99.tar.gz
rm -r bbmap

