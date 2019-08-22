#!/bin/bash

#PBS -N all_hisat2
#PBS -q standby
#PBS -l nodes=1:ppn=20,naccesspolicy=singleuser
#PBS -l walltime=4:00:00

cd $PBS_O_WORKDIR

module purge
module load hisat2/2.1.0

## build index file for reference genome
hisat2-build -p 20 /scratch/snyder/h/harder/ssalar_genome_and_annotation_files/simple_names_chroms_scaf_MT.fa simple_names_chroms_scaf_MT

## map reads to reference genome
hisat2 --dta --time --threads 20 \
-S hisat2_7-b.sam \
-x simple_names_chroms_scaf_MT \
-1 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033774_7-b_S86_R1_filtered.fastq.gz \
-2 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033774_7-b_S86_R2_filtered.fastq.gz

hisat2 --dta --time --threads 20 \
-S hisat2_9-a2.sam \
-x simple_names_chroms_scaf_MT \
-1 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033781_9-a2_S93_R1_filtered.fastq.gz \
-2 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033781_9-a2_S93_R2_filtered.fastq.gz

hisat2 --dta --time --threads 20 \
-S hisat2_11b.sam \
-x simple_names_chroms_scaf_MT \
-1 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033752_11b_S64_R1_filtered.fastq.gz \
-2 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033752_11b_S64_R2_filtered.fastq.gz

hisat2 --dta --time --threads 20 \
-S hisat2_1-b.sam \
-x simple_names_chroms_scaf_MT \
-1 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033746_1-b_S58_R1_filtered.fastq.gz \
-2 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033746_1-b_S58_R2_filtered.fastq.gz

hisat2 --dta --time --threads 20 \
-S hisat2_3-b.sam \
-x simple_names_chroms_scaf_MT \
-1 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033760_3-b_S72_R1_filtered.fastq.gz \
-2 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033760_3-b_S72_R2_filtered.fastq.gz

hisat2 --dta --time --threads 20 \
-S hisat2_8a.sam \
-x simple_names_chroms_scaf_MT \
-1 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033779_8a_S91_R1_filtered.fastq.gz \
-2 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033779_8a_S91_R2_filtered.fastq.gz

hisat2 --dta --time --threads 20 \
-S hisat2_1-a.sam \
-x simple_names_chroms_scaf_MT \
-1 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033745_1-a_S57_R1_filtered.fastq.gz \
-2 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033745_1-a_S57_R2_filtered.fastq.gz

hisat2 --dta --time --threads 20 \
-S hisat2_11a.sam \
-x simple_names_chroms_scaf_MT \
-1 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033751_11a_S63_R1_filtered.fastq.gz \
-2 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033751_11a_S63_R2_filtered.fastq.gz

hisat2 --dta --time --threads 20 \
-S hisat2_4b.sam \
-x simple_names_chroms_scaf_MT \
-1 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033768_4b_S80_R1_filtered.fastq.gz \
-2 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033768_4b_S80_R2_filtered.fastq.gz

hisat2 --dta --time --threads 20 \
-S hisat2_4a2.sam \
-x simple_names_chroms_scaf_MT \
-1 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033767_4a2_S79_R1_filtered.fastq.gz \
-2 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033767_4a2_S79_R2_filtered.fastq.gz

hisat2 --dta --time --threads 20 \
-S hisat2_8b.sam \
-x simple_names_chroms_scaf_MT \
-1 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033780_8b_S92_R1_filtered.fastq.gz \
-2 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033780_8b_S92_R2_filtered.fastq.gz

hisat2 --dta --time --threads 20 \
-S hisat2_9-b.sam \
-x simple_names_chroms_scaf_MT \
-1 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033782_9-b_S94_R1_filtered.fastq.gz \
-2 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033782_9-b_S94_R2_filtered.fastq.gz

hisat2 --dta --time --threads 20 \
-S hisat2_13a2.sam \
-x simple_names_chroms_scaf_MT \
-1 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033755_13a2_S67_R1_filtered.fastq.gz \
-2 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033755_13a2_S67_R2_filtered.fastq.gz

hisat2 --dta --time --threads 20 \
-S hisat2_1b.sam \
-x simple_names_chroms_scaf_MT \
-1 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033748_1b_S60_R1_filtered.fastq.gz \
-2 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033748_1b_S60_R2_filtered.fastq.gz

hisat2 --dta --time --threads 20 \
-S hisat2_7-a.sam \
-x simple_names_chroms_scaf_MT \
-1 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033773_7-a_S85_R1_filtered.fastq.gz \
-2 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033773_7-a_S85_R2_filtered.fastq.gz

hisat2 --dta --time --threads 20 \
-S hisat2_3a.sam \
-x simple_names_chroms_scaf_MT \
-1 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033763_3a_S75_R1_filtered.fastq.gz \
-2 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033763_3a_S75_R2_filtered.fastq.gz

hisat2 --dta --time --threads 20 \
-S hisat2_3-a.sam \
-x simple_names_chroms_scaf_MT \
-1 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033759_3-a_S71_R1_filtered.fastq.gz \
-2 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033759_3-a_S71_R2_filtered.fastq.gz

hisat2 --dta --time --threads 20 \
-S hisat2_3b.sam \
-x simple_names_chroms_scaf_MT \
-1 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033764_3b_S76_R1_filtered.fastq.gz \
-2 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033764_3b_S76_R2_filtered.fastq.gz

hisat2 --dta --time --threads 20 \
-S hisat2_8-b.sam \
-x simple_names_chroms_scaf_MT \
-1 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033778_8-b_S90_R1_filtered.fastq.gz \
-2 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033778_8-b_S90_R2_filtered.fastq.gz

hisat2 --dta --time --threads 20 \
-S hisat2_31b2.sam \
-x simple_names_chroms_scaf_MT \
-1 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033762_31b2_S74_R1_filtered.fastq.gz \
-2 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033762_31b2_S74_R2_filtered.fastq.gz

hisat2 --dta --time --threads 20 \
-S hisat2_4-b2.sam \
-x simple_names_chroms_scaf_MT \
-1 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033766_4-b2_S78_R1_filtered.fastq.gz \
-2 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033766_4-b2_S78_R2_filtered.fastq.gz

hisat2 --dta --time --threads 20 \
-S hisat2_3-1b.sam \
-x simple_names_chroms_scaf_MT \
-1 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033758_3-1b_S70_R1_filtered.fastq.gz \
-2 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033758_3-1b_S70_R2_filtered.fastq.gz

hisat2 --dta --time --threads 20 \
-S hisat2_13-a.sam \
-x simple_names_chroms_scaf_MT \
-1 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033753_13-a_S65_R1_filtered.fastq.gz \
-2 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033753_13-a_S65_R2_filtered.fastq.gz

hisat2 --dta --time --threads 20 \
-S hisat2_1c.sam \
-x simple_names_chroms_scaf_MT \
-1 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033747_1c_S59_R1_filtered.fastq.gz \
-2 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033747_1c_S59_R2_filtered.fastq.gz

hisat2 --dta --time --threads 20 \
-S hisat2_4-a.sam \
-x simple_names_chroms_scaf_MT \
-1 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033765_4-a_S77_R1_filtered.fastq.gz \
-2 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033765_4-a_S77_R2_filtered.fastq.gz

hisat2 --dta --time --threads 20 \
-S hisat2_7b.sam \
-x simple_names_chroms_scaf_MT \
-1 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033776_7b_S88_R1_filtered.fastq.gz \
-2 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033776_7b_S88_R2_filtered.fastq.gz

hisat2 --dta --time --threads 20 \
-S hisat2_8-a.sam \
-x simple_names_chroms_scaf_MT \
-1 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033777_8-a_S89_R1_filtered.fastq.gz \
-2 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033777_8-a_S89_R2_filtered.fastq.gz

hisat2 --dta --time --threads 20 \
-S hisat2_7a.sam \
-x simple_names_chroms_scaf_MT \
-1 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033775_7a_S87_R1_filtered.fastq.gz \
-2 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033775_7a_S87_R2_filtered.fastq.gz

hisat2 --dta --time --threads 20 \
-S hisat2_6-b.sam \
-x simple_names_chroms_scaf_MT \
-1 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033770_6-b_S82_R1_filtered.fastq.gz \
-2 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033770_6-b_S82_R2_filtered.fastq.gz

hisat2 --dta --time --threads 20 \
-S hisat2_9b.sam \
-x simple_names_chroms_scaf_MT \
-1 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033784_9b_S96_R1_filtered.fastq.gz \
-2 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033784_9b_S96_R2_filtered.fastq.gz 

hisat2 --dta --time --threads 20 \
-S hisat2_13b.sam \
-x simple_names_chroms_scaf_MT \
-1 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033756_13b_S68_R1_filtered.fastq.gz \
-2 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033756_13b_S68_R2_filtered.fastq.gz

hisat2 --dta --time --threads 20 \
-S hisat2_9a.sam \
-x simple_names_chroms_scaf_MT \
-1 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033783_9a_S95_R1_filtered.fastq.gz \
-2 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033783_9a_S95_R2_filtered.fastq.gz

hisat2 --dta --time --threads 20 \
-S hisat2_11-a.sam \
-x simple_names_chroms_scaf_MT \
-1 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033749_11-a_S61_R1_filtered.fastq.gz \
-2 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033749_11-a_S61_R2_filtered.fastq.gz

hisat2 --dta --time --threads 20 \
-S hisat2_11-b.sam \
-x simple_names_chroms_scaf_MT \
-1 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033750_11-b_S62_R1_filtered.fastq.gz \
-2 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033750_11-b_S62_R2_filtered.fastq.gz

hisat2 --dta --time --threads 20 \
-S hisat2_13-b.sam \
-x simple_names_chroms_scaf_MT \
-1 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033754_13-b_S66_R1_filtered.fastq.gz \
-2 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033754_13-b_S66_R2_filtered.fastq.gz

hisat2 --dta --time --threads 20 \
-S hisat2_6a.sam \
-x simple_names_chroms_scaf_MT \
-1 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033771_6a_S83_R1_filtered.fastq.gz \
-2 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033771_6a_S83_R2_filtered.fastq.gz

hisat2 --dta --time --threads 20 \
-S hisat2_6-a.sam \
-x simple_names_chroms_scaf_MT \
-1 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033769_6-a_S81_R1_filtered.fastq.gz \
-2 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033769_6-a_S81_R2_filtered.fastq.gz

hisat2 --dta --time --threads 20 \
-S hisat2_31a.sam \
-x simple_names_chroms_scaf_MT \
-1 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033761_31a_S73_R1_filtered.fastq.gz \
-2 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033761_31a_S73_R2_filtered.fastq.gz

hisat2 --dta --time --threads 20 \
-S hisat2_6b.sam \
-x simple_names_chroms_scaf_MT \
-1 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033772_6b_S84_R1_filtered.fastq.gz \
-2 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033772_6b_S84_R2_filtered.fastq.gz

hisat2 --dta --time --threads 20 \
-S hisat2_3-1a2.sam \
-x simple_names_chroms_scaf_MT \
-1 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033757_3-1a2_S69_R1_filtered.fastq.gz \
-2 /scratch/snyder/h/harder/ssalar_rna_seq_data/downloaded_reads/033757_3-1a2_S69_R2_filtered.fastq.gz
