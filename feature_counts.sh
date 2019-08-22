#!/bin/bash

#PBS -N featureCounts
#PBS -q standby
#PBS -l nodes=1:ppn=20,naccesspolicy=singleuser
#PBS -l walltime=4:00:00

cd $PBS_O_WORKDIR

module purge
module load Subread/1.6.1

# -p = is paired-end; -B = require that fragments have both ends successfully aligned; -C
# = chimeric fragments (ends mapped to 2 different chroms) not counted

featureCounts -p -B -C -T 20 -t exon -g gene_id \
-a /scratch/snyder/h/harder/analyses_ssalar_rna_seq/stringtie/stringtie_all_merged.gtf \
-o featureCount_counts.txt \
/scratch/snyder/h/harder/analyses_ssalar_rna_seq/hisat2/hisat2_11-a.sorted.bam \
/scratch/snyder/h/harder/analyses_ssalar_rna_seq/hisat2/hisat2_11a.sorted.bam \
/scratch/snyder/h/harder/analyses_ssalar_rna_seq/hisat2/hisat2_11-b.sorted.bam \
/scratch/snyder/h/harder/analyses_ssalar_rna_seq/hisat2/hisat2_11b.sorted.bam \
/scratch/snyder/h/harder/analyses_ssalar_rna_seq/hisat2/hisat2_13a2.sorted.bam \
/scratch/snyder/h/harder/analyses_ssalar_rna_seq/hisat2/hisat2_13-a.sorted.bam \
/scratch/snyder/h/harder/analyses_ssalar_rna_seq/hisat2/hisat2_13-b.sorted.bam \
/scratch/snyder/h/harder/analyses_ssalar_rna_seq/hisat2/hisat2_13b.sorted.bam \
/scratch/snyder/h/harder/analyses_ssalar_rna_seq/hisat2/hisat2_1-a.sorted.bam \
/scratch/snyder/h/harder/analyses_ssalar_rna_seq/hisat2/hisat2_1-b.sorted.bam \
/scratch/snyder/h/harder/analyses_ssalar_rna_seq/hisat2/hisat2_1b.sorted.bam \
/scratch/snyder/h/harder/analyses_ssalar_rna_seq/hisat2/hisat2_1c.sorted.bam \
/scratch/snyder/h/harder/analyses_ssalar_rna_seq/hisat2/hisat2_3-1a2.sorted.bam \
/scratch/snyder/h/harder/analyses_ssalar_rna_seq/hisat2/hisat2_31a.sorted.bam \
/scratch/snyder/h/harder/analyses_ssalar_rna_seq/hisat2/hisat2_31b2.sorted.bam \
/scratch/snyder/h/harder/analyses_ssalar_rna_seq/hisat2/hisat2_3-1b.sorted.bam \
/scratch/snyder/h/harder/analyses_ssalar_rna_seq/hisat2/hisat2_3-a.sorted.bam \
/scratch/snyder/h/harder/analyses_ssalar_rna_seq/hisat2/hisat2_3a.sorted.bam \
/scratch/snyder/h/harder/analyses_ssalar_rna_seq/hisat2/hisat2_3-b.sorted.bam \
/scratch/snyder/h/harder/analyses_ssalar_rna_seq/hisat2/hisat2_3b.sorted.bam \
/scratch/snyder/h/harder/analyses_ssalar_rna_seq/hisat2/hisat2_4a2.sorted.bam \
/scratch/snyder/h/harder/analyses_ssalar_rna_seq/hisat2/hisat2_4-a.sorted.bam \
/scratch/snyder/h/harder/analyses_ssalar_rna_seq/hisat2/hisat2_4-b2.sorted.bam \
/scratch/snyder/h/harder/analyses_ssalar_rna_seq/hisat2/hisat2_4b.sorted.bam \
/scratch/snyder/h/harder/analyses_ssalar_rna_seq/hisat2/hisat2_6-a.sorted.bam \
/scratch/snyder/h/harder/analyses_ssalar_rna_seq/hisat2/hisat2_6a.sorted.bam \
/scratch/snyder/h/harder/analyses_ssalar_rna_seq/hisat2/hisat2_6-b.sorted.bam \
/scratch/snyder/h/harder/analyses_ssalar_rna_seq/hisat2/hisat2_6b.sorted.bam \
/scratch/snyder/h/harder/analyses_ssalar_rna_seq/hisat2/hisat2_7-a.sorted.bam \
/scratch/snyder/h/harder/analyses_ssalar_rna_seq/hisat2/hisat2_7a.sorted.bam \
/scratch/snyder/h/harder/analyses_ssalar_rna_seq/hisat2/hisat2_7-b.sorted.bam \
/scratch/snyder/h/harder/analyses_ssalar_rna_seq/hisat2/hisat2_7b.sorted.bam \
/scratch/snyder/h/harder/analyses_ssalar_rna_seq/hisat2/hisat2_8-a.sorted.bam \
/scratch/snyder/h/harder/analyses_ssalar_rna_seq/hisat2/hisat2_8a.sorted.bam \
/scratch/snyder/h/harder/analyses_ssalar_rna_seq/hisat2/hisat2_8-b.sorted.bam \
/scratch/snyder/h/harder/analyses_ssalar_rna_seq/hisat2/hisat2_8b.sorted.bam \
/scratch/snyder/h/harder/analyses_ssalar_rna_seq/hisat2/hisat2_9-a2.sorted.bam \
/scratch/snyder/h/harder/analyses_ssalar_rna_seq/hisat2/hisat2_9a.sorted.bam \
/scratch/snyder/h/harder/analyses_ssalar_rna_seq/hisat2/hisat2_9-b.sorted.bam \
/scratch/snyder/h/harder/analyses_ssalar_rna_seq/hisat2/hisat2_9b.sorted.bam
