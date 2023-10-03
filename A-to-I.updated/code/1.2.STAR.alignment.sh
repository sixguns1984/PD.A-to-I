#!/usr/bin/bash

STAR --runThreadN 4 \
--genomeDir genome_index \
--genomeLoad NoSharedMemory \
--outFileNamePrefix XXX.Aligned.sortedByCoord.out.bam \
--outReadsUnmapped Fastx \
--outSAMtype BAM SortedByCoordinate \
--outSAMstrandField intronMotif \
--outSAMattributes All \
--readFilesCommand zcat \
--outFilterType BySJout \
--outFilterMultimapNmax 1 \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 1 \
--outFilterMismatchNmax 999 \
--outFilterMismatchNoverLmax 0.04 \
--alignIntronMin 20 \
--alignIntronMax 1000000 \
--alignMatesGapMax 1000000 \
--readFilesIn clean_data.R1.fastq.gz clean_data.R2.fastq.gz
