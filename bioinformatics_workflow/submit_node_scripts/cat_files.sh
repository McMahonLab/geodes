#!/bin/bash
#Concatenate sortmerna output
mv /home/amlinz/*-rRNA/* /home/amlinz/GEODES_rRNA_split
mv /home/amlinz/*-nonrRNA/* /home/amlinz/GEODES_nonrRNA_split
rmdir *RNA

gzip -d /home/amlinz/GEODES_nonrRNA_split/*
gzip -d /home/amlinz/GEODES_rRNA_split/*

cat /home/amlinz/samplenames.txt | while read line;
  do cat /home/amlinz/GEODES_nonrRNA_split/$line*nonrRNA.fastq > /home/amlinz/GEODES_nonrRNA/$line-nonrRNA.fastq;
  gzip /home/amlinz/GEODES_nonrRNA/$line-nonrRNA.fastq
  cat /home/amlinz/GEODES_rRNA_split/$line*rRNA.fastq > /home/amlinz/GEODES_rRNA/$line-rRNA.fastq;
  gzip /home/amlinz/GEODES_rRNA/$line-rRNA.fastq
done
  
