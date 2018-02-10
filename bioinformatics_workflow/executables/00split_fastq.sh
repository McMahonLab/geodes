#!/bin/bash
tar -xvzf BBMap_36.99.tar.gz

cp /mnt/gluster/amlinz/filtered/$1.filter-MTF.fastq.gz .
gzip -d $1.filter-MTF.fastq.gz

sed -i '/^$/d' $1.filter-MTF.fastq

maxreads=$((`wc -l < $1.filter-MTF.fastq` / 8 - 1))
startpoints=$(seq 0 500000 $maxreads)

for num in $startpoints;
  do endpoint=$(($num + 499999));
  bbmap/getreads.sh in=$1.filter-MTF.fastq id=$num-$endpoint out=$1-$endpoint.fastq overwrite=T;
  done

rm $1.filter-MTF.fastq
mkdir $1-splitfiles
mv $1*.fastq $1-splitfiles
tar cvzf $1-splitfiles.tar.gz $1-splitfiles
rm BBMap_36.99.tar.gz
rm -r bbmap
