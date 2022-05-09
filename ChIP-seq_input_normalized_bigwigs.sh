#!/bin/bash

#Load relevant modules
module load python
source activate /n/holystore01/LABS/liau_lab/Users/nlue/conda_envs/deeptools_env

#Define directories
CHIP_DIR=[...]
BAM_DIR=[...]

#Use DeepTools bamCompare to generate input-normalized bigwig files
cd $CHIP_DIR
bamCompare -b1 $BAM_DIR/GAGTGG.TKO_H3K4me3.deduplicated.bam -b2 $BAM_DIR/CGTACG.TKO_input.deduplicated.bam -o $CHIP_DIR/H3K4me3_input_log2fc_200bp_bins_1kb_smoothed.bw -p max --extendReads --operation log2 --skipZeroOverZero -bs 200 --smoothLength 1000 --scaleFactorsMethod SES
bamCompare -b1 $BAM_DIR/ACTGAT.TKO_H3K27me3.deduplicated.bam -b2 $BAM_DIR/CGTACG.TKO_input.deduplicated.bam -o $CHIP_DIR/H3K27me3_input_log2fc_200bp_bins_1kb_smoothed.bw -p max --extendReads --operation log2 --skipZeroOverZero -bs 200 --smoothLength 1000 --scaleFactorsMethod SES
bamCompare -b1 $BAM_DIR/ATTCCT.TKO_H3K36me2.deduplicated.bam -b2 $BAM_DIR/CGTACG.TKO_input.deduplicated.bam -o $CHIP_DIR/H3K36me2_input_log2fc_200bp_bins_1kb_smoothed.bw -p max --extendReads --operation log2 --skipZeroOverZero -bs 200 --smoothLength 1000 --scaleFactorsMethod SES

