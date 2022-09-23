#!/bin/bash

#wget ftp://ftp.wormbase.org/pub/wormbase/releases/WS283/species/c_elegans/PRJNA13758/c_elegans.PRJNA13758.WS283.genomic.fa.gz

bwa index c_elegans.PRJNA13758.WS283.genomic.fa.gz

for SAMPLE in SRR16356854 SRR16356855 SRR16356856
do
	bwa mem \
	-R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}" \
	-t 4 \
	-o ${SAMPLE}.sam \
	c_elegans.PRJNA13758.WS283.genomic.fa.gz ${SAMPLE}_1.fastq.gz ${SAMPLE}_2.fastq.gz
done

for SAMPLE in SRR16356854 SRR16356855 SRR16356856
do
	samtools sort -@4 -O bam -o ${SAMPLE}.bam ${SAMPLE}.sam
done

for SAMPLE in SRR16356854 SRR16356855 SRR16356856
do
	samtools index ${SAMPLE}.bam
done

gzcat c_elegans.PRJNA13758.WS283.genomic.fa.gz > c_elegans.PRJNA13758.WS283.genomic.fa

freebayes -f c_elegans.PRJNA13758.WS283.genomic.fa \
	SRR16356854.bam SRR16356855.bam SRR16356856.bam > c_elegans.vcf
