#!/bin/bash

##### Use Defiant to call DMRs #####

#Define key directories
OUT_DIR=[...] #directory to output to
IN_DIR=[...] #directory containing bedGraph files with methylation (after intersection to retain only CpGs at 5x coverage across all samples)
DEFIANT_DIR=[...] #directory containing defiant

cd $OUT_DIR

$DEFIANT_DIR/defiant -b -N -v -c 5 -p 0.05 -s 5 -CpN 3 -d 2 -P 4 -S 2 -G 5000 -l R297W -L WT,R297W -i $IN_DIR/A1_5_coverage_ABCDE_intersected.bedGraph,$IN_DIR/A2_5_coverage_ABCDE_intersected.bedGraph $IN_DIR/B1_5_coverage_ABCDE_intersected.bedGraph,$IN_DIR/B2_5_coverage_ABCDE_intersected.bedGraph
$DEFIANT_DIR/defiant -b -N -v -c 5 -p 0.05 -s 5 -CpN 3 -d 2 -P 4 -S 2 -G 5000 -l W326R -L WT,W326R -i $IN_DIR/A1_5_coverage_ABCDE_intersected.bedGraph,$IN_DIR/A2_5_coverage_ABCDE_intersected.bedGraph $IN_DIR/C1_5_coverage_ABCDE_intersected.bedGraph,$IN_DIR/C2_5_coverage_ABCDE_intersected.bedGraph
$DEFIANT_DIR/defiant -b -N -v -c 5 -p 0.05 -s 5 -CpN 3 -d 2 -P 4 -S 2 -G 5000 -l E338K -L WT,E338K -i $IN_DIR/A1_5_coverage_ABCDE_intersected.bedGraph,$IN_DIR/A2_5_coverage_ABCDE_intersected.bedGraph $IN_DIR/D1_5_coverage_ABCDE_intersected.bedGraph,$IN_DIR/D2_5_coverage_ABCDE_intersected.bedGraph
$DEFIANT_DIR/defiant -b -N -v -c 5 -p 0.05 -s 5 -CpN 3 -d 2 -P 4 -S 2 -G 5000 -l E752K -L WT,E752K -i $IN_DIR/A1_5_coverage_ABCDE_intersected.bedGraph,$IN_DIR/A2_5_coverage_ABCDE_intersected.bedGraph $IN_DIR/E1_5_coverage_ABCDE_intersected.bedGraph,$IN_DIR/E2_5_coverage_ABCDE_intersected.bedGraph

#Further process the output bed files and separate up and down DMRs
#First remove header line of the bed files
for FILE in *.bed; do sed -i '1d' $FILE; done
#Filter out increasing or decreasing DMRs
for FILE in *S2.tsv; do
	sed '1d' $FILE | awk '$6 > 0 {print "chr"$1 "\t" $2 "\t" $3}' > categorized_DMRs/${FILE//c5_CpN3_d2_G5000_p0.05_P4_S2.tsv/higher.bed}
	sed '1d' $FILE | awk '$6 < 0 {print "chr"$1 "\t" $2 "\t" $3}' > categorized_DMRs/${FILE//c5_CpN3_d2_G5000_p0.05_P4_S2.tsv/lower.bed}
done

#Note: renamed files manually to '220731_R297W_lower.bed', etc.


##### Use BedTools to perform intersections with genomic annotations #####

cd $OUT_DIR/genome_annotations

#Make bed file with CpG islands
#220731_mm10_cpgIslandExt.txt downloaded from UCSC
sed '1d' 220731_mm10_cpgIslandExt.txt | awk '{print $2 "\t" $3 "\t" $4 "\t" "CpG_"$6 "\t" $7}' > 220731_mm10_CGIs.bed

#Make bed file with genes
#220731_mm10_gencode_vm23_knownGene.txt downloaded from UCSC
sed '1d' 220731_mm10_gencode_vm23_knownGene.txt | awk '{print $2 "\t" $4 "\t" $5}' > 220731_mm10_gencode23_genes.bed
#Sort and merge
LC_COLLATE=C sort -k1,1 -k2,2n 220731_mm10_gencode23_genes.bed > 220731_mm10_gencode23_genes_sorted.bed
bedtools merge -i 220731_mm10_gencode23_genes_sorted.bed > 220731_mm10_gencode23_genes_merged.bed

#Take complement of genic regions to get intergenic regions
LC_COLLATE=C sort -k1,1 mm10.chrom.sizes > mm10.sorted.chrom.sizes
LC_COLLATE=C sort -k1,1 -k2,2n 220731_mm10_gencode23_genes_merged.bed | bedtools complement -i stdin -g mm10.sorted.chrom.sizes > 220731_mm10_intergenic.bed

#Use bedtools intersect to count # DMRs that overlap each annotation
cd $OUT_DIR/categorized_DMRs

#First, filter out noncanonical chromosomes
for FILE in *.bed; do
	awk '$1 !~ /_/' $FILE | grep -vw 'chrM' > ${FILE//.bed/_filtered.bed}
done

#Intersection
bedtools intersect -u -a 220731_R297W_higher_filtered.bed -b $OUT_DIR/genome_annotations/220731_mm10_CGIs.bed > overlaps/R297W_higher_CGI_unique.bed
bedtools intersect -u -a 220731_R297W_higher_filtered.bed -b $OUT_DIR/genome_annotations/220731_mm10_gencode23_genes_merged.bed > overlaps/R297W_higher_genic_unique.bed
bedtools intersect -u -a 220731_R297W_higher_filtered.bed -b $OUT_DIR/genome_annotations/220731_mm10_intergenic.bed > overlaps/R297W_higher_intergenic_unique.bed

bedtools intersect -u -a 220731_R297W_lower_filtered.bed -b $OUT_DIR/genome_annotations/220731_mm10_CGIs.bed > overlaps/R297W_lower_CGI_unique.bed
bedtools intersect -u -a 220731_R297W_lower_filtered.bed -b $OUT_DIR/genome_annotations/220731_mm10_gencode23_genes_merged.bed > overlaps/R297W_lower_genic_unique.bed
bedtools intersect -u -a 220731_R297W_lower_filtered.bed -b $OUT_DIR/genome_annotations/220731_mm10_intergenic.bed > overlaps/R297W_lower_intergenic_unique.bed


bedtools intersect -u -a 220731_W326R_higher_filtered.bed -b $OUT_DIR/genome_annotations/220731_mm10_CGIs.bed > overlaps/W326R_higher_CGI_unique.bed
bedtools intersect -u -a 220731_W326R_higher_filtered.bed -b $OUT_DIR/genome_annotations/220731_mm10_gencode23_genes_merged.bed > overlaps/W326R_higher_genic_unique.bed
bedtools intersect -u -a 220731_W326R_higher_filtered.bed -b $OUT_DIR/genome_annotations/220731_mm10_intergenic.bed > overlaps/W326R_higher_intergenic_unique.bed

bedtools intersect -u -a 220731_W326R_lower_filtered.bed -b $OUT_DIR/genome_annotations/220731_mm10_CGIs.bed > overlaps/W326R_lower_CGI_unique.bed
bedtools intersect -u -a 220731_W326R_lower_filtered.bed -b $OUT_DIR/genome_annotations/220731_mm10_gencode23_genes_merged.bed > overlaps/W326R_lower_genic_unique.bed
bedtools intersect -u -a 220731_W326R_lower_filtered.bed -b $OUT_DIR/genome_annotations/220731_mm10_intergenic.bed > overlaps/W326R_lower_intergenic_unique.bed


bedtools intersect -u -a 220731_E338K_higher_filtered.bed -b $OUT_DIR/genome_annotations/220731_mm10_CGIs.bed > overlaps/E338K_higher_CGI_unique.bed
bedtools intersect -u -a 220731_E338K_higher_filtered.bed -b $OUT_DIR/genome_annotations/220731_mm10_gencode23_genes_merged.bed > overlaps/E338K_higher_genic_unique.bed
bedtools intersect -u -a 220731_E338K_higher_filtered.bed -b $OUT_DIR/genome_annotations/220731_mm10_intergenic.bed > overlaps/E338K_higher_intergenic_unique.bed

bedtools intersect -u -a 220731_E338K_lower_filtered.bed -b $OUT_DIR/genome_annotations/220731_mm10_CGIs.bed > overlaps/E338K_lower_CGI_unique.bed
bedtools intersect -u -a 220731_E338K_lower_filtered.bed -b $OUT_DIR/genome_annotations/220731_mm10_gencode23_genes_merged.bed > overlaps/E338K_lower_genic_unique.bed
bedtools intersect -u -a 220731_E338K_lower_filtered.bed -b $OUT_DIR/genome_annotations/220731_mm10_intergenic.bed > overlaps/E338K_lower_intergenic_unique.bed


bedtools intersect -u -a 220731_E752K_higher_filtered.bed -b $OUT_DIR/genome_annotations/220731_mm10_CGIs.bed > overlaps/E752K_higher_CGI_unique.bed
bedtools intersect -u -a 220731_E752K_higher_filtered.bed -b $OUT_DIR/genome_annotations/220731_mm10_gencode23_genes_merged.bed > overlaps/E752K_higher_genic_unique.bed
bedtools intersect -u -a 220731_E752K_higher_filtered.bed -b $OUT_DIR/genome_annotations/220731_mm10_intergenic.bed > overlaps/E752K_higher_intergenic_unique.bed

bedtools intersect -u -a 220731_E752K_lower_filtered.bed -b $OUT_DIR/genome_annotations/220731_mm10_CGIs.bed > overlaps/E752K_lower_CGI_unique.bed
bedtools intersect -u -a 220731_E752K_lower_filtered.bed -b $OUT_DIR/genome_annotations/220731_mm10_gencode23_genes_merged.bed > overlaps/E752K_lower_genic_unique.bed
bedtools intersect -u -a 220731_E752K_lower_filtered.bed -b $OUT_DIR/genome_annotations/220731_mm10_intergenic.bed > overlaps/E752K_lower_intergenic_unique.bed















