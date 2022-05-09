#!/bin/bash

##### Perform bwa alignment of raw ChIP-seq sequencing files #####

#Define key files/directories
GENOME=[.../mm10.fa] #mm10 genome fasta file from UCSC
INDEX_DIR=[...] #bwa index directory
FASTQ_DIR=[...] #directory containing fastq files
ALIGN_DIR=[...] #alignment output directory

SAMPLE_ARRAY=('CGTACG' 'GAGTGG' 'ACTGAT' 'ATTCCT')
NAME_ARRAY=('TKO_input' 'TKO_H3K4me3' 'TKO_H3K27me3' 'TKO_H3K36me2')
INDEX=0

bwa index -p $INDEX_DIR/mm10_bwa -a bwtsw $GENOME

for bc in "${SAMPLE_ARRAY[@]}"; do
	echo $bc
	echo $INDEX
	echo ${NAME_ARRAY[$INDEX]}

	bwa aln -t 16 $INDEX_DIR/mm10_bwa $FASTQ_DIR/1_HMVW3DRXY.1.${bc}.unmapped.1.fastq.gz > $ALIGN_DIR/${bc}.aln_sa.1.sai

	bwa aln -t 16 $INDEX_DIR/mm10_bwa $FASTQ_DIR/1_HMVW3DRXY.1.${bc}.unmapped.2.fastq.gz > $ALIGN_DIR/${bc}.aln_sa.2.sai

	bwa sampe $INDEX_DIR/mm10_bwa $ALIGN_DIR/${bc}.aln_sa.1.sai $ALIGN_DIR/${bc}.aln_sa.2.sai \
	$FASTQ_DIR/1_HMVW3DRXY.1.${bc}.unmapped.1.fastq.gz \
	$FASTQ_DIR/1_HMVW3DRXY.1.${bc}.unmapped.2.fastq.gz \
	> $ALIGN_DIR/${bc}.${NAME_ARRAY[$INDEX]}.aln_pe.sam

	((++INDEX))

done

##### Sort and deduplicate aligned reads #####

#Define key directories
PICARD_DIR=[...] #directory containing picard.jar file
SAM_DIR=$ALIGN_DIR
OUT_DIR=[...] #output directory for sorted, deduplicated bam files

SAM_ARRAY=('CGTACG.TKO_input.aln_pe.sam' 'GAGTGG.TKO_H3K4me3.aln_pe.sam' 'ACTGAT.TKO_H3K27me3.aln_pe.sam' 'ATTCCT.TKO_H3K36me2.aln_pe.sam')

for file in "${SAM_ARRAY[@]}"; do
	echo $file

	#Use samtools view function to convert .sam to .bam (compress data)
	samtools view -S -b $SAM_DIR/${file} > $SAM_DIR/${file//pe.sam/pe.bam}

	#Use picard SortSam function to sort by coordinate
	java -jar $PICARD_DIR/picard.jar SortSam \
		I=$SAM_DIR/${file//pe.sam/pe.bam} \
		O=$OUT_DIR/${file//aln_pe.sam/sorted.bam} \
		SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT

	#Use picard MarkDuplicates to mark and remove duplicate alignments
	java -jar $PICARD_DIR/picard.jar MarkDuplicates \
		I=$OUT_DIR/${file//aln_pe.sam/sorted.bam} \
		O=$OUT_DIR/${file//aln_pe.sam/deduplicated.bam} \
		M=$OUT_DIR/${file//aln_pe.sam/marked_dup_metrics.txt} \
		REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT
done

##### Generate input-normalized bigwig files #####

#Define directories
CHIP_DIR=[...] #output directory for bigwig files
BAM_DIR=$OUT_DIR

#Use DeepTools bamCompare to generate input-normalized bigwig files
cd $CHIP_DIR
bamCompare -b1 $BAM_DIR/GAGTGG.TKO_H3K4me3.deduplicated.bam -b2 $BAM_DIR/CGTACG.TKO_input.deduplicated.bam -o $CHIP_DIR/H3K4me3_input_log2fc_200bp_bins_1kb_smoothed.bw -p max --extendReads --operation log2 --skipZeroOverZero -bs 200 --smoothLength 1000 --scaleFactorsMethod SES
bamCompare -b1 $BAM_DIR/ACTGAT.TKO_H3K27me3.deduplicated.bam -b2 $BAM_DIR/CGTACG.TKO_input.deduplicated.bam -o $CHIP_DIR/H3K27me3_input_log2fc_200bp_bins_1kb_smoothed.bw -p max --extendReads --operation log2 --skipZeroOverZero -bs 200 --smoothLength 1000 --scaleFactorsMethod SES
bamCompare -b1 $BAM_DIR/ATTCCT.TKO_H3K36me2.deduplicated.bam -b2 $BAM_DIR/CGTACG.TKO_input.deduplicated.bam -o $CHIP_DIR/H3K36me2_input_log2fc_200bp_bins_1kb_smoothed.bw -p max --extendReads --operation log2 --skipZeroOverZero -bs 200 --smoothLength 1000 --scaleFactorsMethod SES


