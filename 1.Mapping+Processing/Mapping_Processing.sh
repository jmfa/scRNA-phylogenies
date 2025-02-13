#!/bin/sh

# define paths for resources and data

reference=/path/to/reference_genome_hg38
fastqs=/path/to/fastqs

# Map FASTQs to hg38 using STAR

STAR \
--runThreadN 5 \
--genomeDir ${reference} \
--outFileNamePrefix ${cellID} \
--outSAMtype BAM SortedByCoordinate \
--readFilesIn ${fastqs}/${sraID}_1.fastq ${fastqs}/${sraID}_2.fastq


# Dedup with Picard

temp=/path/to/temp_folder

java -jar picard.jar MarkDuplicates \
	I=${cellID}Aligned.sortedByCoord.out.bam \
	OUTPUT=${cellID}.dedup.bam \
	CREATE_INDEX=true \
	REMOVE_DUPLICATES=false \
	TMP_DIR=${temp} \
	M=${cellID}_Metrics.txt \
	VALIDATION_STRINGENCY=LENIENT


