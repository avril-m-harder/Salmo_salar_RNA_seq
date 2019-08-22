# Salmo_salar_RNA_seq
In run order:

1. hisat2.sh: Builds HISAT2 index for reference genome and maps reads to reference genome. Outputs SAM file per sample.

2. samtools_example.sh: Sorts SAM files and writes sorted BAM files.

3. stringtie_example.sh: Assembles transcripts for each sample, merges all transcripts from all samples and the reference annotation, and compares assembled sample transcripts to known transcripts. Writes merged GTF file.

4. feature_counts.sh: Generates text file with count matrix for all samples using sorted BAM files and merged GTF file.
