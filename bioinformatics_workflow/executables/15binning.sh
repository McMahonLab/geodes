#!/bin/bash
# Test binning on one metagenome assembly
# Both the assembled data and the reads have already been slimmed down
tar xvzf MaxBin.tar.gz

export PATH=$(pwd)/MaxBin-2.2.4/auxiliary/FragGeneScan1.30:$PATH
export PATH=$(pwd)/MaxBin-2.2.4/auxiliary/bowtie2-2.2.3:$PATH
export PATH=$(pwd)/MaxBin-2.2.4/auxiliary/hmmer-3.1b1/src:$PATH
export PATH=$(pwd)/MaxBin-2.2.4/auxiliary/idba-1.1.1/bin:$PATH

cp /mnt/gluster/amlinz/filtered_assemblies/$1-filtered.assembled.fna.gz .
cp /mnt/gluster/amlinz/downsampled_metagenomes/$1-sampled.fastq.gz .

gzip -d $1-filtered.assembled.fna.gz
gzip -d $1-sampled.fastq.gz

./MaxBin-2.2.4/run_MaxBin.pl -contig $1-filtered.assembled.fna -out $1-binned -reads $1-sampled.fastq

mkdir $1-binning
mv $1-binned* $1-binning/
tar cvzf $1-binning.tar.gz $1-binning/
mv $1-binning.tar.gz /mnt/gluster/amlinz/

rm *assembled.fna
rm *fastq
