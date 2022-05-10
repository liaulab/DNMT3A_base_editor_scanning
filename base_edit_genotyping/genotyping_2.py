#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Analysis of genotyping for sgRNA validation experiments shown in Fig. 1i, 2f,
and 4c and Extended Data Fig. 4b,d, as well as construction of allele table in
Supplementary Table 4.

@author: Nicholas Lue
"""

#%% Import Packages

import os
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

#%% Inputs

# Define directory where CRISPResso2 analysis is stored
data_dir = 'CRISPResso2_output'

# Define output directory for this analysis
out_dir= 'Analysis'

# Sample names
sample_list = ['sg118', 'sg128', 'sg135', 'sg136', 'sg137', 'sg149', 'sg151',
               'sg164', 'sg165', 'sg168', 'sg175', 'sg336', 'sg428', 'sg456',
               'sg458', 'sg522', 'sg566', 'sg588', 'sg651']
# Guide names
guide_list = ['DNMT3A2_'+sample_list[i][-3:] for i in range(len(sample_list))]

# Read in mapping statistics
df_map = pd.read_csv(
    os.path.join(data_dir,'CRISPRessoBatch_mapping_statistics.txt'), sep='\t')

# Read in quantification of editing frequency file
df_edits = pd.read_csv(
    os.path.join(
        data_dir,'CRISPRessoBatch_quantification_of_editing_frequency.txt'),
    sep='\t', index_col=0)

# Read in sample-specific information (allele and nucleotide tables)
# Initialize dictionaries to hold sample-specific dataframes
allele_tables_dict = {}
nt_tables_dict = {}

for sample in sample_list:
    
    sample_dir = os.path.join(data_dir, '_'.join(['CRISPResso_on', sample]))
    guide_id = 'DNMT3A2_' + sample[-3:]    
    file = os.path.join(
        sample_dir, ('Alleles_frequency_table_around_sgRNA_' + guide_id + '.txt'))
    
    # Read .txt file and store in dictionary of data frames
    allele_tables_dict[sample] = pd.read_csv(file, sep='\t')
    temp = pd.read_csv(
        os.path.join(sample_dir, 'Nucleotide_frequency_table.txt'),
        sep='\t', header=None, index_col=0)
    temp.index = ['nt'] + temp.index.tolist()[1:] # Replace NaN index with 'nt'
    temp.iloc[1:] = temp.iloc[1:].astype(float) # Make values float, not str
    nt_tables_dict[sample] = temp.rename(index={})

# Import batch settings
batch_settings = pd.read_csv('Batch_settings_2.txt', sep='\t', index_col=0)



#%% Define functions for next sections

# Function to translate codon
def get_aa(codon):
    codon_dict = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
        'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'
    }
    return codon_dict[codon]

# Function to reverse complement sequence
def reverse_complement(seq):
    rev_seq = seq[::-1].upper()
    new_seq = ''
    for base in rev_seq:
        if base == 'A':
            new_seq += 'T'
        elif base == 'T':
            new_seq += 'A'
        elif base == 'G':
            new_seq += 'C'
        elif base == 'C':
            new_seq += 'G'
    return new_seq

# Function to make certain parts lowercase
def make_lower_case(seq, start, end):
    subseq = seq[start:end+1].lower()
    newseq = seq[:start] + subseq + seq[end+1:]
    return newseq

# Function to take allele and turn it into sense direction and mark introns
def annot_allele(inseq, strand, intron_start, intron_end):
    if np.isnan(intron_start):
        if strand == 'sense':
            return inseq
        else:
            return reverse_complement(inseq)
        
    else:
        if strand == 'sense':
            return make_lower_case(inseq, int(intron_start), int(intron_end))
        else:
            temp = reverse_complement(inseq)
            return make_lower_case(temp, int(intron_start), int(intron_end))


# Translate
# Define function to translate allele
def translate_allele(inseq, tr_i, tr_f):
    
    # Isolate the part to be translated
    translate_seq = inseq[tr_i:tr_f+1]
    
    # Split into codons
    codon_list = [translate_seq[i:i+3] for i in range(0,len(translate_seq),3)]
    
    # Initialize aa sequence
    out_aa = ''
    # Translate
    for codon in codon_list:
        out_aa += get_aa(codon)
    return out_aa



#%% Allele table analysis
# This section presents code used to construct Supplementary Table 4.

# Filter to keep only alleles with >1% allele frequence
allele_tables_filtered = {}
for key,value in allele_tables_dict.items():
    
    df_temp = value.copy()
    df_temp = df_temp.loc[df_temp['%Reads']>=1]
    allele_tables_filtered[key] = df_temp.copy()

# Key samples info
df_in = pd.read_csv('Analysis/input.csv')
key_samples = df_in['sgRNA'].tolist()
df_in = df_in.set_index('sgRNA')
            
# Make corrected alleles and translate
df_processed = pd.DataFrame(columns=['sgRNA','Allele','Product','%Reads'])

for key,df in allele_tables_filtered.items():
    if key not in key_samples:
        continue
    
    df_temp = df.copy()
    
    # Pull out info from df_in
    sg_strand = df_in.loc[key,'strand']
    tr_i = df_in.loc[key,'translate_start']
    tr_f = df_in.loc[key,'translate_end']
    in_i = df_in.loc[key,'intron_start']
    in_f = df_in.loc[key,'intron_end']
    
    # Correct allele
    df_temp = df_temp.assign(Allele=df_temp.apply(
        lambda x: annot_allele(x.Aligned_Sequence,sg_strand,in_i,in_f), axis=1))
    
    # Make translation
    df_temp = df_temp.assign(Product=df_temp.apply(
        lambda x: translate_allele(x.Allele,tr_i,tr_f), axis=1))
    
    # Annotate guide
    df_temp['sgRNA'] = key
    
    # Merge
    df_temp = df_temp[['sgRNA','Allele','Product','%Reads']]
    df_processed = pd.concat([df_processed,df_temp],ignore_index=True)
    

# Export to csv
df_processed.to_csv('Analysis/Allele_table_analysis.csv',index=False)


#%% Extract C to T conversion data
# This section is for calculating C to T conversion for use in heatmaps
# (e.g., Fig. 4c) or dotplots (e.g., Fig. 1i).

# Define the positions in nt_frequency dataframes where each guide is found
# Number corresponds to protospacer position +1 in the amplicon, 1-indexed
# Remember that for each guide, the amplicon is indexed in the same orientation
# as the guide.
sg_positions = {'sg118':82, 'sg128':78, 'sg135':109, 'sg136':110, 'sg137':114,
                'sg149':77, 'sg151':79, 'sg164':94, 'sg165':93, 'sg168':90,
                'sg175':116, 'sg336':103, 'sg428':41, 'sg456':57, 'sg458':46,
                'sg522':57, 'sg566':66, 'sg588':59, 'sg651':46
                }

# Define total number of aligned reads
reads_dict = {}
for sample in sample_list:
    reads_dict[sample] = df_edits.loc[sample,'Reads_aligned']

# Initialize dataframe to store data
df_data = pd.DataFrame(columns=[i for i in range(-5,21)])

def calculate_edits(sample_name, df_nts, tot_reads, pos):
    
    temp = []
    
    # Take protospacer positions -6 to 20, inclusive
    for col in range(pos-6,pos+20):
        
        if df_nts.loc['nt',col] != 'C':
            temp.append(0.0)
        else:
            num_T = df_nts.loc['T',col]
            percent_T = num_T / tot_reads * 100
            temp.append(percent_T)
    
    return pd.Series(data=temp, index=[i for i in range(-5,21)],
                     name=sample_name)

# Assemble dataframe
for sample in sample_list:
    df_temp = calculate_edits(sample_name=sample,
                              df_nts=nt_tables_dict[sample],
                              tot_reads=reads_dict[sample],
                              pos=sg_positions[sample]
                              )
    df_data = df_data.append(df_temp)

#%% Extract sequence around each protospacer

# Initialize dataframe to store sequences
df_annot = pd.DataFrame(columns=[i for i in range(-5,21)])

def extract_seqs(sample_name, df_nts, pos):
    
    temp = df_nts.loc['nt',pos-6:pos+19].copy()
    temp.index = [i for i in range(-5,21)]
    temp.name = sample_name
    
    return temp

# Assemble dataframe
for sample in sample_list:
    df_temp = extract_seqs(sample_name=sample,
                           df_nts=nt_tables_dict[sample],
                           pos=sg_positions[sample]
                           )
    df_annot = df_annot.append(df_temp)


#%% Extract indel data and store in df_data

# Create a new column in df_edits containing total # of reads with indels
df_edits['Contains_indels'] = df_edits['Modified'].subtract(
    df_edits['Only Substitutions'])

# Calculate % of aligned reads that contain indels
df_edits['Percent_indel'] = df_edits['Contains_indels'].divide(
    df_edits['Reads_aligned'])*100

# Add new column to df_data with indel info
df_data['Indel'] = df_edits['Percent_indel']

# Add new column to df_annot that is blank
df_annot['Indel'] = ' '

#%% Make heatmaps

def plot_heatmap(df_data, out_prefix, sns_context='talk',
                 cmap=sns.light_palette('#ED202D',as_cmap=True),
                 annot=None, v_min=0, v_max=100):
    
    # Plotting parameters and variables
    sns.set_context(sns_context)
    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['ps.fonttype'] = 42    
    mpl.rcParams['font.sans-serif'] = ['Arial']
    
    # Generate heatmap
    fig, ax = plt.subplots(figsize=(5,4))
    sns.heatmap(data=df_data, ax=ax, square=True, cmap=cmap, center=0,
                linewidths=0.5, linecolor='black',
                annot=annot, fmt='', annot_kws={"fontsize":6,"color":'black'},
                xticklabels=True, yticklabels=True, vmin=v_min, vmax=v_max,
                cbar_kws={"shrink": .70})
    for edge,spine in ax.spines.items():
        spine.set_visible(True)
        spine.set_color('k')
        spine.set_linewidth(0.5)
    ax.set_xlabel('Protospacer position')
    ax.set_ylabel('Guide')
    ax.set_xticklabels(ax.get_xticklabels(), fontsize=8)
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=8)
    ax.tick_params(width=0.5, length=2.5)
    plt.tight_layout()
    plt.savefig(out_prefix+'_editing_heatmap.pdf', format='pdf', bbox_inches='tight')
    plt.close()


# Plot heatmap for all guides
plot_heatmap(df_data=df_data, out_prefix='All',
             cmap='RdBu_r', annot=df_annot.to_numpy())

# Plot heatmap for only hits
df_onlyhits = df_data.copy()
df_onlyhits = df_onlyhits.drop(index=['sg136', 'sg165'])
df_annot_onlyhits = df_annot.copy()
df_annot_onlyhits = df_annot_onlyhits.drop(index=['sg136', 'sg165'])
plot_heatmap(df_data=df_onlyhits, out_prefix='Onlyhits',
             cmap='RdBu_r', annot=df_annot_onlyhits.to_numpy())


