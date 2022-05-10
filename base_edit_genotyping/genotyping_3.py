#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Analysis of sgE342.1, .2, and .3 genotyping data (Extended Data Fig. 6b)
@author: nlclarinet8
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
sample_list = ['sg164_uns', 'sg164_pos',
               'sg165_uns', 'sg165_pos',
               'sg168_uns', 'sg168_pos']

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
    guide_id = 'DNMT3A2_' + sample[2:5]
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
batch_settings = pd.read_csv('Batch_settings_3.txt', sep='\t', index_col=0)

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

# Define a function to truncate allele
def truncate_allele(in_allele, positions):
    
    pos_i, pos_f = positions
    return in_allele[pos_i:pos_f]

# Define function to parse allele into triplets
def parse_allele(in_allele, revcomp=False, codon_tup=None):
    # Reverse complement input if necessary
    if revcomp:
        in_allele = reverse_complement(in_allele)
    # Break input into codons
    codon_list = [in_allele[i:i+3] for i in range(0,len(in_allele),3)]
    # Truncate if necessary
    if codon_tup == None:
        codon_i, codon_f = (0, len(codon_list))
    else:
        codon_i, codon_f = codon_tup
    codon_list = codon_list[codon_i:codon_f]
    # Make output string
    out_seq = ''
    for codon in codon_list:
        out_seq += ' '
        out_seq += codon
    return out_seq[1:]

# Define function to translate allele
def translate_allele(in_parsed):
    # Split into codons
    codon_list = in_parsed.split(' ')
    # Initialize aa sequence
    out_aa = ''
    # Translate
    for codon in codon_list:
        if len(codon) != 3:
            return 'Frameshift'
        else:
            out_aa += get_aa(codon)
    return out_aa

def check_be(in_seq, poslist):
    out_str = ''
    for pos in poslist:
        if in_seq[pos] == 'C':
            out_str += '.'
        elif in_seq[pos] == '-':
            out_str += '-'
        else:
            out_str += 'X'
    return out_str



#%% Analysis of alleles around E153

dict_positions = {'sg164_uns':(1,36), 'sg164_pos':(1,36),
                  'sg165_uns':(2,37), 'sg165_pos':(2,37),
                  'sg168_uns':(5,40), 'sg168_pos':(5,40)}

# Initialize
top_alleles = [] #List of top alleles (>=1% in at least one sample)
merged_allele_tables_dict = {} #Dictionary of processed allele tables

for sample in sample_list:
    
    # Retrieve allele table and make new column with just common subsequences
    df_temp = allele_tables_dict[sample].copy()
    df_temp = df_temp.assign(Allele=df_temp.apply(
        lambda x: truncate_allele(x.Aligned_Sequence, dict_positions[sample]),
        axis=1))
    
    # Combine alleles that are the same and sort the table by %Reads
    df_temp = df_temp.groupby('Allele')['%Reads'].sum().reset_index()
    df_temp = df_temp.sort_values(by='%Reads', ascending=False)    
    
    # Save into dict
    merged_allele_tables_dict[sample] = df_temp.copy()

    # Save top alleles
    df_top = df_temp.loc[df_temp['%Reads']>=1].copy()
    list_temp = df_top['Allele'].tolist()
    
    # Add alleles that have not already been added to top_alleles
    top_alleles.extend([allele for allele in list_temp if allele not in top_alleles])

# Initialize dataframe to store allele data
df_E153K = pd.DataFrame(top_alleles, columns=['Allele'])

# Filter and merge allele data
for sample in sample_list:
    
    # Get dataframe and rename columns
    df_temp = merged_allele_tables_dict[sample].copy()
    df_temp.columns = ['Allele', sample]
    
    # Merge
    df_E153K = df_E153K.merge(df_temp, how='left', on='Allele', validate='1:1')
    
# Change NaN values to 0
df_E153K = df_E153K.fillna(0)

# Sort
df_E153K['avg'] = df_E153K[sample_list].mean(axis=1)
df_E153K = df_E153K.sort_values(by='avg', ascending=False, ignore_index=True)

#%% Annotate new columns

# Get coding subsequence of allele and put spaces between triplets
df_E153K = df_E153K.assign(Allele_parsed=df_E153K.apply(
    lambda x: parse_allele(x.Allele[:24], revcomp=True), axis=1))

# Translate
df_E153K = df_E153K.assign(aa_seq=df_E153K.apply(
    lambda x: translate_allele(x.Allele_parsed), axis=1))

# Base editing summary
be_pos = [0,3,6,9,12,14,17,19,21,23,24,33] #Indices of base edited Cs in Allele string
df_E153K = df_E153K.assign(be_summary=df_E153K.apply(
    lambda x:check_be(x.Allele,be_pos), axis=1))

# Calculate fold-change and log2 fold-change
df_E153K['fc_164'] = df_E153K['sg164_pos'].divide(df_E153K['sg164_uns'])
df_E153K = df_E153K.assign(lfc_164=df_E153K.apply(
    lambda x: np.log2(x.fc_164), axis=1))
df_E153K['fc_165'] = df_E153K['sg165_pos'].divide(df_E153K['sg165_uns'])
df_E153K = df_E153K.assign(lfc_165=df_E153K.apply(
    lambda x: np.log2(x.fc_165), axis=1))
df_E153K['fc_168'] = df_E153K['sg168_pos'].divide(df_E153K['sg168_uns'])
df_E153K = df_E153K.assign(lfc_168=df_E153K.apply(
    lambda x: np.log2(x.fc_168), axis=1))

# Prepare file and export
list_columns = ['Allele', 'Allele_parsed', 'aa_seq', 'be_summary',
                'sg164_uns', 'sg164_pos', 'lfc_164',
                'sg165_uns', 'sg165_pos', 'lfc_165',
                'sg168_uns', 'sg168_pos', 'lfc_168']
df_E153K = df_E153K[list_columns]
df_E153K.to_csv('E153K_allele_freq.csv', index=False)

#%% Barplot function

def plot_bars(df_data, out_prefix):
        
    # Plotting parameters and variables
    sns.set_context('talk')
    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['ps.fonttype'] = 42    
    mpl.rcParams['font.sans-serif'] = ['Arial']
    
    # Generate heatmap
    fig, ax = plt.subplots(figsize=(6,6))
    sns.barplot(data=df_data, ax=ax,
                y='be_summary', x='value', hue='variable')
    plt.tight_layout()
    plt.savefig(out_prefix+'_allele_barplot.pdf', format='pdf', bbox_inches='tight')
    plt.close()

#%% Make barplot

# Put data in long form separately for each guide
df_164_melt = df_E153K.melt(id_vars=['Allele','be_summary','aa_seq'],
                            value_vars=['sg164_uns','sg164_pos'])
df_165_melt = df_E153K.melt(id_vars=['Allele','be_summary','aa_seq'],
                            value_vars=['sg165_uns','sg165_pos'])
df_168_melt = df_E153K.melt(id_vars=['Allele','be_summary','aa_seq'],
                            value_vars=['sg168_uns','sg168_pos'])


plot_bars(df_data=df_164_melt, out_prefix=out_dir+'/sg164')
plot_bars(df_data=df_165_melt, out_prefix=out_dir+'/sg165')
plot_bars(df_data=df_168_melt, out_prefix=out_dir+'/sg168')

#%% Make log2 fold-change heatmap

def plot_heatmap(df_data, out_prefix, v_min=-2, v_max=+2, sns_context='talk',
                 cmap=sns.light_palette('#ED202D',as_cmap=True)):
    
    # Plotting parameters and variables
    sns.set_context(sns_context)
    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['ps.fonttype'] = 42    
    mpl.rcParams['font.sans-serif'] = ['Arial']
    
    # Generate heatmap
    fig, ax = plt.subplots(figsize=(6,6))
    sns.heatmap(data=df_data, ax=ax, square=True, cmap=cmap, center=0,
                linewidths=0.5, linecolor='black', vmin=v_min, vmax=v_max,
                annot=True, fmt='.3f', annot_kws={"fontsize":6,"color":'black'},
                xticklabels=True, yticklabels=True,
                cbar_kws={"shrink":.70})
    for edge,spine in ax.spines.items():
        spine.set_visible(True)
        spine.set_color('k')
        spine.set_linewidth(0.5)
    ax.tick_params(width=0.5, length=2.5)
    plt.tight_layout()
    plt.savefig(out_prefix+'_lfc_heatmap.pdf', format='pdf', bbox_inches='tight')
    plt.close()

df_lfc = df_E153K[['be_summary', 'lfc_164', 'lfc_168', 'lfc_165']].copy()
df_lfc = df_lfc.set_index('be_summary')
plot_heatmap(df_data=df_lfc, out_prefix='sgE153K', v_min=-3, v_max=3,
             cmap='coolwarm_r')



