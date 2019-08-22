#!/bin/bash

#PBS -N stringtie
#PBS -q darwin
#PBS -l nodes=1:ppn=4,naccesspolicy=singleuser
#PBS -l walltime=300:00:00

cd $PBS_O_WORKDIR

module purge
module load stringtie/1.3.3
module load gffcompare/0.10.2


# create assembly per sample
stringtie /scratch/snyder/h/harder/analyses_ssalar_rna_seq/hisat2/hisat2_11-a.sorted.bam -l 11-a -p 4 -G /scratch/snyder/h/harder/ssalar_genome_and_annotation_files/use_this_one_ref_ICSASG_v2_top_level.gtf -o hisat2_11-a.sorted.bam.gtf
stringtie /scratch/snyder/h/harder/analyses_ssalar_rna_seq/hisat2/hisat2_11a.sorted.bam -l 11a -p 4 -G /scratch/snyder/h/harder/ssalar_genome_and_annotation_files/use_this_one_ref_ICSASG_v2_top_level.gtf -o hisat2_11a.sorted.bam.gtf
stringtie /scratch/snyder/h/harder/analyses_ssalar_rna_seq/hisat2/hisat2_11-b.sorted.bam -l 11-b -p 4 -G /scratch/snyder/h/harder/ssalar_genome_and_annotation_files/use_this_one_ref_ICSASG_v2_top_level.gtf -o hisat2_11-b.sorted.bam.gtf
stringtie /scratch/snyder/h/harder/analyses_ssalar_rna_seq/hisat2/hisat2_11b.sorted.bam -l 11b -p 4 -G /scratch/snyder/h/harder/ssalar_genome_and_annotation_files/use_this_one_ref_ICSASG_v2_top_level.gtf -o hisat2_11b.sorted.bam.gtf
stringtie /scratch/snyder/h/harder/analyses_ssalar_rna_seq/hisat2/hisat2_13a2.sorted.bam -l 13a2 -p 4 -G /scratch/snyder/h/harder/ssalar_genome_and_annotation_files/use_this_one_ref_ICSASG_v2_top_level.gtf -o hisat2_13a2.sorted.bam.gtf
stringtie /scratch/snyder/h/harder/analyses_ssalar_rna_seq/hisat2/hisat2_13-a.sorted.bam -l 13-a -p 4 -G /scratch/snyder/h/harder/ssalar_genome_and_annotation_files/use_this_one_ref_ICSASG_v2_top_level.gtf -o hisat2_13-a.sorted.bam.gtf
stringtie /scratch/snyder/h/harder/analyses_ssalar_rna_seq/hisat2/hisat2_13-b.sorted.bam -l 13-b -p 4 -G /scratch/snyder/h/harder/ssalar_genome_and_annotation_files/use_this_one_ref_ICSASG_v2_top_level.gtf -o hisat2_13-b.sorted.bam.gtf
stringtie /scratch/snyder/h/harder/analyses_ssalar_rna_seq/hisat2/hisat2_13b.sorted.bam -l 13b -p 4 -G /scratch/snyder/h/harder/ssalar_genome_and_annotation_files/use_this_one_ref_ICSASG_v2_top_level.gtf -o hisat2_13b.sorted.bam.gtf


## once all sample assemblies have been generated:

# merge all transcripts from all samples + reference annotation
stringtie --merge -p 4 \
-G /scratch/snyder/h/harder/ssalar_genome_and_annotation_files/use_this_one_ref_ICSASG_v2_top_level.gtf \
-o stringtie_all_merged.gtf mergelist.txt


# compare assembled transcripts to known transcripts
gffcompare -r \
/scratch/snyder/h/harder/ssalar_genome_and_annotation_files/use_this_one_ref_ICSASG_v2_top_level.gtf \
-G -o stringtie_all_merged stringtie_all_merged.gtf
