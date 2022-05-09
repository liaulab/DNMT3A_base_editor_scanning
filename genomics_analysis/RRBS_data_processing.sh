#!/bin/bash

##### Trim and align fastq files #####
#Representative commands for trimming of raw sequencing files. These were repeated for the other replicates/samples.

#Define key directories
FASTQ_DIR=[...] #directory containing raw sequencing reads
TRIM_DIR=[...] #directory for trimmed reads
GENOME_DIR=[...] #directory containing mm10 genome prepared with bismark
LAMBDA_DIR=[...] #directory containing lambda phage genome prepared with bismark
ALIGN_DIR=[...] #directory for alignments to mm10 genome
ALIGN_LAMBDA_DIR=[...] #directory for alignments to lambda phage genome

#Define arrays/variables for barcodes/names/etc
SAMPLE_ARRAY=('A1' 'B1' 'C1')
LANE_2=('ACTCTC' 'TGAGCT' 'GAGACG')
NAMES_2=('A1' 'B1' 'C1')
INDEX_2=0

#Trim using trim_galore and also run fastQC on trimmed files
for bc in "${LANE_2[@]}"; do
	trim_galore --illumina --rrbs --paired --length 21 --fastqc \
				--output_dir $TRIM_DIR --basename ${NAMES_2[$INDEX_2]} \
				$FASTQ_DIR/2_HTJ3FDRXY.2.${bc}.unmapped.2.fastq.gz $FASTQ_DIR/2_HTJ3FDRXY.2.${bc}.unmapped.1.fastq.gz
	((++INDEX_2))
done

#Align using bismark
for samp in "${SAMPLE_ARRAY[@]}"; do
	#Align first to mm10
	bismark --parallel 4 --bowtie2 --un --ambiguous --output_dir $ALIGN_DIR \
	--genome $GENOME_DIR -1 $TRIM_DIR/${samp}_val_1.fq.gz -2 $TRIM_DIR/${samp}_val_2.fq.gz
	#Then align to lambda phage genome
	bismark --parallel 4 --bowtie2 --un --ambiguous --output_dir $ALIGN_LAMBDA_DIR \
	--genome $LAMBDA_DIR -1 $TRIM_DIR/${samp}_val_1.fq.gz -2 $TRIM_DIR/${samp}_val_2.fq.gz
done



##### Extract methylation using MethylDackel #####

#Define key directories
ALIGN_DIR=[...] #same directory as above containing mm10 alignments
SORT_DIR=[...] #directory for sorted bam files
OUT_DIR=[...] #directory for extracted methylation information
GENOME=[.../mm10.fa] #mm10 genome file

#Define arrays/variables for barcodes/names/etc
SAMPLE_ARRAY=('A1' 'B1' 'C1' 'D1' 'E1' 'A2' 'B2' 'C2' 'D2' 'E2')

#Sort and index
for samp in "${SAMPLE_ARRAY[@]}"; do
	samtools sort -o $SORT_DIR/${samp}_resorted.bam $ALIGN_DIR/${samp}_val_1_bismark_bt2_pe.bam
	samtools index $SORT_DIR/${samp}_resorted.bam
done

#Index mm10 with faidx
samtools faidx $GENOME

#Extract methylation with methyldackel
for samp in "${SAMPLE_ARRAY[@]}"; do
	MethylDackel extract -d 5 --mergeContext --keepDupes -o $OUT_DIR/${samp}_5_coverage $GENOME $SORT_DIR/${samp}_resorted.bam
done



##### Process bedGraph files and make bigWigs #####

BEDGRAPH_DIR=$OUT_DIR #directory containing extracted methylation information
BIGWIG_DIR=[...] #directory for bigwig files
SCRIPT_DIR=[...] #directory containing bedGraphToBigWig script and mm10.chrom.sizes file

#Define arrays/variables for barcodes/names/etc
SAMPLE_ARRAY=('A1' 'B1' 'C1' 'D1' 'E1' 'A2' 'B2' 'C2' 'D2' 'E2')

#Process bedgraph files
for samp in "${SAMPLE_ARRAY[@]}"; do cat $BEDGRAPH_DIR/${samp}_5_coverage_CpG.bedGraph | sed '1d' | awk '{print $1 "\t" $2 "\t" $3 "\t" $5/($5+$6) "\t" $5 "\t" $6}' > $BIGWIG_DIR/${samp}_5_coverage_CpG.meth.bedGraph; done

#Make bigwigs
cd $BIGWIG_DIR
for FILE in *.bedGraph; do LC_COLLATE=C sort -k1,1 -k2,2n $FILE | awk '{print $1 "\t" $2 "\t" $3 "\t" $4}' > ${FILE//CpG.meth.bedGraph/sorted_4col_CpG.bedGraph}; done
for FILE in *4col_CpG.bedGraph; do $SCRIPT_DIR/bedGraphToBigWig $FILE $SCRIPT_DIR/mm10.chrom.sizes ${FILE//sorted_4col_CpG.bedGraph/CpG.bw}; done



##### Intersect bedGraph files to retain only CpGs with 5x coverage in all samples #####

#Intersection
cd $BIGWIG_DIR
bedtools intersect -a A1_5_coverage_CpG.meth.bedGraph -b A2_5_coverage_CpG.meth.bedGraph > intersections/A12_5x_intersection.bedGraph
bedtools intersect -a B1_5_coverage_CpG.meth.bedGraph -b B2_5_coverage_CpG.meth.bedGraph > intersections/B12_5x_intersection.bedGraph
bedtools intersect -a C1_5_coverage_CpG.meth.bedGraph -b C2_5_coverage_CpG.meth.bedGraph > intersections/C12_5x_intersection.bedGraph
bedtools intersect -a D1_5_coverage_CpG.meth.bedGraph -b D2_5_coverage_CpG.meth.bedGraph > intersections/D12_5x_intersection.bedGraph
bedtools intersect -a E1_5_coverage_CpG.meth.bedGraph -b E2_5_coverage_CpG.meth.bedGraph > intersections/E12_5x_intersection.bedGraph
bedtools intersect -a intersections/A12_5x_intersection.bedGraph -b intersections/B12_5x_intersection.bedGraph | bedtools intersect -a stdin -b intersections/C12_5x_intersection.bedGraph | bedtools intersect -a stdin -b intersections/D12_5x_intersection.bedGraph | bedtools intersect -a stdin -b intersections/E12_5x_intersection.bedGraph > intersections/ABCDE_5x_intersection.bedGraph

#Make bed file containing only CpGs with 5x coverage in all samples/replicates
cat intersections/ABCDE_5x_intersection.bedGraph | awk '{print $1 "\t" $2 "\t" $3}' > intersections/ABCDE_5x_intersection.bed

#Make intersected bedgraphs
for FILE in *.meth.bedGraph; do bedtools intersect -a $FILE -b intersections/ABCDE_5x_intersection.bedGraph > intersections/${FILE//CpG.meth.bedGraph/ABCDE_intersected.bedGraph}; done

#Make bigwigs of intersected files
cd intersections
for FILE in *intersected.bedGraph; do LC_COLLATE=C sort -k1,1 -k2,2n $FILE | awk '{print $1 "\t" $2 "\t" $3 "\t" $4}' > ${FILE//.bedGraph/_sorted_4col_CpG.bedGraph}; done
for FILE in *4col_CpG.bedGraph; do $SCRIPT_DIR/bedGraphToBigWig $FILE $SCRIPT_DIR/mm10.chrom.sizes ${FILE//_sorted_4col_CpG.bedGraph/.bw}; done



##### Use DeepTools multiBigwigSummary to generate .tsv files to use in downstream analysis #####

#Define key directories
CHIP_DIR=/n/holystore01/LABS/liau_lab/Users/nlue/NZL10273/220323_bigwigs
OUT_DIR_2=/n/holystore01/LABS/liau_lab/Users/nlue/NZL10280/220323_misc_analysis/CpG_level_global_analysis
BIGWIG_DIR=/n/holystore01/LABS/liau_lab/Users/nlue/NZL10280/220323_misc_analysis/rrbs_bigwigs

cd $BIGWIG_DIR

# Methylation at CpG-level. Only CpGs with 5x coverage across all samples.
multiBigwigSummary BED-file -b A1_5_coverage_CpG.bw A2_5_coverage_CpG.bw B1_5_coverage_CpG.bw B2_5_coverage_CpG.bw C1_5_coverage_CpG.bw C2_5_coverage_CpG.bw D1_5_coverage_CpG.bw D2_5_coverage_CpG.bw E1_5_coverage_CpG.bw E2_5_coverage_CpG.bw \
	-l WT_1 WT_2 R297W_1 R297W_2 W326R_1 W326R_2 E338K_1 E338K_2 E752K_1 E752K_2 --BED intersections/ABCDE_5x_intersection.bed -o $OUT_DIR_2/220326_ABCDE_global_methylation.npz --outRawCounts $OUT_DIR_2/220326_ABCDE_global_methylation.tsv -p max

cd intersections/

# Methylation in 10 kb bins. Only CpGs with 5x coverage across all samples.
multiBigwigSummary bins --binSize 10000 -p max --smartLabels -b A1_5_coverage_ABCDE_intersected.bw \
	A2_5_coverage_ABCDE_intersected.bw \
	B1_5_coverage_ABCDE_intersected.bw \
	B2_5_coverage_ABCDE_intersected.bw \
	C1_5_coverage_ABCDE_intersected.bw \
	C2_5_coverage_ABCDE_intersected.bw \
	D1_5_coverage_ABCDE_intersected.bw \
	D2_5_coverage_ABCDE_intersected.bw \
	E1_5_coverage_ABCDE_intersected.bw \
	E2_5_coverage_ABCDE_intersected.bw \
	$CHIP_DIR/H3K4me3_input_log2fc_200bp_bins_1kb_smoothed.bw $CHIP_DIR/H3K27me3_input_log2fc_200bp_bins_1kb_smoothed.bw $CHIP_DIR/H3K36me2_input_log2fc_200bp_bins_1kb_smoothed.bw \
	-o 220326_ABCDE_intersection_global_methylation.npz --outRawCounts 220326_ABCDE_intersection_global_methylation.tsv

# Methylation in 500 kb bins. Only CpGs with 5x coverage across all samples.
multiBigwigSummary bins --binSize 500000 -p max --smartLabels -b A1_5_coverage_ABCDE_intersected.bw A2_5_coverage_ABCDE_intersected.bw \
	B1_5_coverage_ABCDE_intersected.bw B2_5_coverage_ABCDE_intersected.bw C1_5_coverage_ABCDE_intersected.bw C2_5_coverage_ABCDE_intersected.bw \
	D1_5_coverage_ABCDE_intersected.bw D2_5_coverage_ABCDE_intersected.bw E1_5_coverage_ABCDE_intersected.bw E2_5_coverage_ABCDE_intersected.bw \
	-o 220408_ABCDE_intersection_500kb_bins_global_methylation.npz --outRawCounts 220408_ABCDE_intersection_500kb_bins_global_methylation.tsv























