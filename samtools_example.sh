#!/bin/bash

#PBS -N samtools_1
#PBS -q darwin
#PBS -l nodes=1:ppn=2,naccesspolicy=shared
#PBS -l walltime=300:00:00
#PBS -m ae
#PBS -M avrilharder@gmail.com

cd $PBS_O_WORKDIR

module purge
module load samtools/1.8

samtools view -@ 2 -Su hisat2_7-b.sam | samtools sort -@ 2 -o hisat2_7-b.sorted.bam 
samtools view -@ 2 -Su hisat2_9-a2.sam | samtools sort -@ 2 -o hisat2_9-a2.sorted.bam 
samtools view -@ 2 -Su hisat2_11b.sam | samtools sort -@ 2 -o hisat2_11b.sorted.bam 
samtools view -@ 2 -Su hisat2_1-b.sam | samtools sort -@ 2 -o hisat2_1-b.sorted.bam