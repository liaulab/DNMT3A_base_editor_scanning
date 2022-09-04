#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Process fastq files and calculate normalized sgRNA counts in each replicate
and condition. Further calculate 'sgRNA scores.'

@author: Nicholas Lue and Kevin Ngan
"""

#%% Import Packages

import os
import pandas as pd

from ScreenAnalysisFunctions_v2 import (
    batch_count, batch_process, compare_conditions
)

#%% Key Inputs

path = os.getcwd()
in_ref = 'Annotated_Library_Full.csv' #Output of annotation script
in_batch = 'Batch_info.csv' #Spreadsheet with sample/replicate ID information


#%% Read and process fastq files

# Define paths
subdir_counts = os.path.join('Count_Results', 'counts')
subdir_np = os.path.join('Count_Results', 'npcounts')
subdir_stats = os.path.join('Count_Results', 'stats')

# Count sgRNAs in each sample
# Input fastq files placed in directory 'Raw_Data'
batch_count(in_batch=in_batch, in_ref=in_ref, subdir_counts=subdir_counts,
            subdir_np=subdir_np, subdir_stats=subdir_stats, dir_fastq='Raw_Data')

#%% Preprocessing

# Process counts (merge files, normalize and log2 transform counts, normalize to day 0)
df_ref = batch_process(in_batch=in_batch, in_ref=in_ref,
                       out_prefix='DNMT3A_scanning_', out_folder='Processed',
                       counts_path=os.path.join(path,subdir_counts),
                       stats_path=os.path.join(path,subdir_stats),
                       pDNA_ID='NL-P', day0_ID='NL0')

#%% Compare conditions

# Read in batch info file
df_batch = pd.read_csv(in_batch)

# List of conditions
list_conds = df_batch['condition'].unique().tolist()
list_conds.remove('t0')

# Define comparisons (sorted vs. unsorted)
# Example: 'd9-pos' comparison refers to d9 citrine-positive condition normalized
# to d9 citrine-unsorted condition
list_comps = [
    (list_conds[2],list_conds[1]), (list_conds[3],list_conds[1]),
    (list_conds[5],list_conds[4]), (list_conds[6],list_conds[4]),
    (list_conds[8],list_conds[7]), (list_conds[9],list_conds[7])
    ]
list_compnames = [
    'd3-pos', 'd3-neg',
    'd6-pos', 'd6-neg',
    'd9-pos', 'd9-neg'
    ]
list_comparisons = [(name, cond[0], cond[1]) for name,cond in zip(list_compnames, list_comps)]

compare_conditions(comparisons=list_comparisons,
                   in_t0norm=os.path.join(path,'Processed','DNMT3A_scanning_t0norm_conds.csv'),
                   df_ref=df_ref, out_comps='DNMT3A_scanning_comparisons.csv')



#%% Create output file for paper
# Export .csv file with processed read counts (at replicate level)

df_reps = pd.read_csv('Processed/DNMT3A_scanning_t0norm_reps.csv')

# Dictionary mapping sample names
samples_dict = {'NL0':'day0',
                'NL1':'day3_unsorted_rep1',
                'NL10':'day3_unsorted_rep2',
                'NL19':'day3_unsorted_rep3',
                'NL2':'day3_citrine_pos_rep1',
                'NL11':'day3_citrine_pos_rep2',
                'NL20':'day3_citrine_pos_rep3',
                'NL3':'day3_citrine_neg_rep1',
                'NL12':'day3_citrine_neg_rep2',
                'NL21':'day3_citrine_neg_rep3',
                'NL4':'day6_unsorted_rep1',
                'NL13':'day6_unsorted_rep2',
                'NL22':'day6_unsorted_rep3',
                'NL5':'day6_citrine_pos_rep1',
                'NL14':'day6_citrine_pos_rep2',
                'NL23':'day6_citrine_pos_rep3',
                'NL6':'day6_citrine_neg_rep1',
                'NL15':'day6_citrine_neg_rep2',
                'NL24':'day6_citrine_neg_rep3',
                'NL7':'day9_unsorted_rep1',
                'NL16':'day9_unsorted_rep2',
                'NL25':'day9_unsorted_rep3',
                'NL8':'day9_citrine_pos_rep1',
                'NL17':'day9_citrine_pos_rep2',
                'NL26':'day9_citrine_pos_rep3',
                'NL9':'day9_citrine_neg_rep1',
                'NL18':'day9_citrine_neg_rep2',
                'NL27':'day9_citrine_neg_rep3'}

# Rename columns
df_reps = df_reps.rename(mapper=samples_dict, axis=1)

# Columns to keep
key_cols = ['sgRNA_ID', 'sgRNA_seq', 'Gene', 'Mut_type', 'Domain',
            'Mut_list_3A1', 'Edit_site_3A1'] + list(samples_dict.values())

# Export processed data for Supplementary Data 2
df_reps[key_cols].to_csv('Processed/DNMT3A_scanning_t0norm_reps_keycols.csv', index=False)



#%% Calculation of sgRNA scores
# Here, calculate sgRNA scores and prepare data for plotting

# Read in comparisons
df_data = pd.read_csv('DNMT3A_scanning_comparisons.csv')

# Prep data for plotting and define key variables

# Define list of annotations
list_annocols = ['sgRNA_ID', 'Gene', 'Targeted_exon', 'C_count', 'is_C',
                 'Splice_check', 'Mut_type', 'Edit_site_3A1', 'Domain', 'Is_domain',
                 'Submut_type']

# Define list of comparisons
list_compnames = ['d3-pos', 'd3-neg', 'd6-pos', 'd6-neg', 'd9-pos', 'd9-neg']

# Normalize data to intergenic controls to calculate 'sgRNA scores'
# First, use negative controls to set cutoffs
df_negctrl = df_data.loc[df_data['Gene'] == 'NON-GENE'].copy()
list_negctrlstats = [] # list of tups of (comp, avg, sd, avg+2sd, avg-2sd)
cutoff_dict = {} # dictionary of comp: (avg+2sd, avg-2sd)
avg_dict = {} # dictionary of comp: avg
for comp in list_compnames:
    temp = (df_negctrl[comp].mean(), df_negctrl[comp].std())
    tup_comp = (comp, temp[0], temp[1], temp[0] + (2*temp[1]), temp[0] - (2*temp[1]))
    list_negctrlstats.append(tup_comp)
    cutoff_dict[comp] = (tup_comp[3], tup_comp[4])
    avg_dict[comp] = tup_comp[1]
# Second, perform normalization
df_not_normalized = df_data.copy()
for comp in list_compnames:
    df_data[comp] = df_not_normalized[comp].sub(avg_dict[comp])


# Export sgRNA scores
df_data.to_csv('DNMT3A_scanning_sgRNA_scores.csv', index=False)
