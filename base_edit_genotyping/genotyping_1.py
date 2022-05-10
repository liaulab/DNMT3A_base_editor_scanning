#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Analysis of sgW698, sgE756, and sgG532 genotyping data.

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
sample_list = ['NL28', 'NL29', 'NL30', 'NL31',
               'NL32', 'NL33', 'NL34', 'NL35',
               'NL36', 'NL37', 'NL38', 'NL39']

# Make dictionary for guide names
guide_dict = {'NL28':'GGTAGCCGTCGTCGTCGTAC', 'NL29':'GGTAGCCGTCGTCGTCGTAC',
              'NL30':'GGTAGCCGTCGTCGTCGTAC', 'NL31':'GGTAGCCGTCGTCGTCGTAC',
              'NL32':'ATTCTCAAAGAGCCAGAAGA', 'NL33':'ATTCTCAAAGAGCCAGAAGA',
              'NL34':'ATTCTCAAAGAGCCAGAAGA', 'NL35':'ATTCTCAAAGAGCCAGAAGA',
              'NL36':'GCCCCACTCCTGGATCTGGG', 'NL37':'GCCCCACTCCTGGATCTGGG',
              'NL38':'GCCCCACTCCTGGATCTGGG', 'NL39':'GCCCCACTCCTGGATCTGGG'}

# Make dictionary for sample IDs
samp_dict = {'NL30':'sgG532_uns', 'NL31':'sgG532_neg',
             'NL32':'sgE756_uns', 'NL33':'sgE756_pos',
             'NL36':'sgW698_uns', 'NL37':'sgW698_pos'}

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
    guide_id = guide_dict[sample]
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
batch_settings = pd.read_csv('Batch_settings_1.txt', sep='\t', index_col=0)

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
    
def plot_bars(df_data, out_prefix):
        
    # Plotting parameters and variables
    sns.set_context('talk')
    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['ps.fonttype'] = 42    
    mpl.rcParams['font.sans-serif'] = ['Arial']
    
    # Generate heatmap
    fig, ax = plt.subplots(figsize=(5,6))
    sns.barplot(data=df_data, ax=ax,
                y='be_summary', x='value', hue='variable')
    plt.tight_layout()
    plt.savefig(out_prefix+'_allele_barplot.pdf', format='pdf', bbox_inches='tight')
    plt.close()


#%% Analysis of alleles around W698

# Define samples to analyze
list_W698 = ['NL36', 'NL37']

# Initialize
top_W698 = [] #List of top alleles (>=1% in at least one sample)
merged_W698_alleles_dict = {} #Dictionary of processed allele tables

for sample in list_W698:
    
    # Retrieve allele table and make new column with just common subsequences
    df_temp = allele_tables_dict[sample].copy()
    df_temp = df_temp.assign(Allele=df_temp.apply(
        lambda x: truncate_allele(x.Aligned_Sequence, (10,28)), axis=1))
    
    # Combine alleles that are the same and sort the table by %Reads
    df_temp = df_temp.groupby('Allele')['%Reads'].sum().reset_index()
    df_temp = df_temp.sort_values(by='%Reads', ascending=False)    
    
    # Save into dict
    merged_W698_alleles_dict[sample] = df_temp.copy()

    # Save top alleles
    df_top = df_temp.loc[df_temp['%Reads']>=1].copy()
    list_temp = df_top['Allele'].tolist()
    
    # Add alleles that have not already been added to top_alleles
    top_W698.extend([allele for allele in list_temp if allele not in top_W698])

# Initialize dataframe to store allele data
df_W698 = pd.DataFrame(top_W698, columns=['Allele'])

# Filter and merge allele data
for sample in list_W698:
    
    # Get dataframe and rename columns
    df_temp = merged_W698_alleles_dict[sample].copy()
    df_temp.columns = ['Allele', samp_dict[sample]]
    
    # Merge
    df_W698 = df_W698.merge(df_temp, how='left', on='Allele', validate='1:1')
    
# Change NaN values to 0
df_W698 = df_W698.fillna(0)

# Sort
df_W698 = df_W698.sort_values(by='sgW698_uns', ascending=False, ignore_index=True)

# Annotate new columns
# Get coding subsequence, flip to coding sense, and put spaces between triplets
df_W698 = df_W698.assign(Allele_parsed=df_W698.apply(
    lambda x: parse_allele(x.Allele,revcomp=True),axis=1))

# Translate
df_W698 = df_W698.assign(aa_seq=df_W698.apply(
    lambda x: translate_allele(x.Allele_parsed), axis=1))

# Base editing summary
be_pos = [1,2,3,4,6,8,9,15] #indices of base edited Cs in Allele
df_W698 = df_W698.assign(be_summary=df_W698.apply(
    lambda x:check_be(x.Allele,be_pos), axis=1))

# Calculate fold-change and then log2 fold-change
df_W698['fc'] = df_W698['sgW698_pos'].divide(df_W698['sgW698_uns'])
df_W698 = df_W698.assign(lfc=df_W698.apply(
    lambda x: np.log2(x.fc), axis=1))

# Prepare file and export
df_W698 = df_W698[['Allele','Allele_parsed','aa_seq','be_summary',
                   'sgW698_uns','sgW698_pos','lfc']]
df_W698.to_csv('W698_allele_freq.csv', index=False)

#%% Analysis of alleles around E756

# Define samples to analyze
list_E756 = ['NL32', 'NL33']

# Initialize
top_E756 = [] #List of top alleles (>=1% in at least one sample)
merged_E756_alleles_dict = {} #Dictionary of processed allele tables

for sample in list_E756:
    
    # Retrieve allele table and make new column with just common subsequences
    df_temp = allele_tables_dict[sample].copy()
    df_temp = df_temp.assign(Allele=df_temp.apply(
        lambda x: truncate_allele(x.Aligned_Sequence, (10,28)), axis=1))
    
    # Combine alleles that are the same and sort the table by %Reads
    df_temp = df_temp.groupby('Allele')['%Reads'].sum().reset_index()
    df_temp = df_temp.sort_values(by='%Reads', ascending=False)    
    
    # Save into dict
    merged_E756_alleles_dict[sample] = df_temp.copy()

    # Save top alleles
    df_top = df_temp.loc[df_temp['%Reads']>=1].copy()
    list_temp = df_top['Allele'].tolist()
    
    # Add alleles that have not already been added to top_alleles
    top_E756.extend([allele for allele in list_temp if allele not in top_E756])

# Initialize dataframe to store allele data
df_E756 = pd.DataFrame(top_E756, columns=['Allele'])

# Filter and merge allele data
for sample in list_E756:
    
    # Get dataframe and rename columns
    df_temp = merged_E756_alleles_dict[sample].copy()
    df_temp.columns = ['Allele', samp_dict[sample]]
    
    # Merge
    df_E756 = df_E756.merge(df_temp, how='left', on='Allele', validate='1:1')
    
# Change NaN values to 0
df_E756 = df_E756.fillna(0)

# Sort
df_E756 = df_E756.sort_values(by='sgE756_uns', ascending=False, ignore_index=True)

# Annotate new columns
# Get coding subsequence, flip to coding sense, and put spaces between triplets
df_E756 = df_E756.assign(Allele_parsed=df_E756.apply(
    lambda x: parse_allele(x.Allele,revcomp=True),axis=1))

# Translate
df_E756 = df_E756.assign(aa_seq=df_E756.apply(
    lambda x: translate_allele(x.Allele_parsed), axis=1))

# Base editing summary
be_pos = [3,5,12,13] #indices of base edited Cs in Allele
df_E756 = df_E756.assign(be_summary=df_E756.apply(
    lambda x:check_be(x.Allele,be_pos), axis=1))

# Calculate fold-change and then log2 fold-change
df_E756['fc'] = df_E756['sgE756_pos'].divide(df_E756['sgE756_uns'])
df_E756 = df_E756.assign(lfc=df_E756.apply(
    lambda x: np.log2(x.fc), axis=1))

# Prepare file and export
df_E756 = df_E756[['Allele','Allele_parsed','aa_seq','be_summary',
                   'sgE756_uns','sgE756_pos','lfc']]
df_E756.to_csv('E756_allele_freq.csv', index=False)

#%% Analysis of alleles around G532

# Define samples to analyze
list_G532 = ['NL30', 'NL31']

# Initialize
top_G532 = [] #List of top alleles (>=1% in at least one sample)
merged_G532_alleles_dict = {} #Dictionary of processed allele tables

for sample in list_G532:
    
    # Retrieve allele table and make new column with just common subsequences
    df_temp = allele_tables_dict[sample].copy()
    df_temp = df_temp.assign(Allele=df_temp.apply(
        lambda x: truncate_allele(x.Aligned_Sequence, (10,28)), axis=1))
    
    # Combine alleles that are the same and sort the table by %Reads
    df_temp = df_temp.groupby('Allele')['%Reads'].sum().reset_index()
    df_temp = df_temp.sort_values(by='%Reads', ascending=False)    
    
    # Save into dict
    merged_G532_alleles_dict[sample] = df_temp.copy()

    # Save top alleles
    df_top = df_temp.loc[df_temp['%Reads']>=1].copy()
    list_temp = df_top['Allele'].tolist()
    
    # Add alleles that have not already been added to top_alleles
    top_G532.extend([allele for allele in list_temp if allele not in top_G532])

# Initialize dataframe to store allele data
df_G532 = pd.DataFrame(top_G532, columns=['Allele'])

# Filter and merge allele data
for sample in list_G532:
    
    # Get dataframe and rename columns
    df_temp = merged_G532_alleles_dict[sample].copy()
    df_temp.columns = ['Allele', samp_dict[sample]]
    
    # Merge
    df_G532 = df_G532.merge(df_temp, how='left', on='Allele', validate='1:1')
    
# Change NaN values to 0
df_G532 = df_G532.fillna(0)

# Sort
df_G532 = df_G532.sort_values(by='sgG532_uns', ascending=False, ignore_index=True)

# Annotate new columns
# Get coding subsequence, flip to coding sense, and put spaces between triplets
df_G532 = df_G532.assign(Allele_parsed=df_G532.apply(
    lambda x: parse_allele(x.Allele,revcomp=True),axis=1))

# Translate
df_G532 = df_G532.assign(aa_seq=df_G532.apply(
    lambda x: translate_allele(x.Allele_parsed), axis=1))

# Base editing summary
be_pos = [5,6,9,12,15] #indices of base edited Cs in Allele
df_G532 = df_G532.assign(be_summary=df_G532.apply(
    lambda x:check_be(x.Allele,be_pos), axis=1))

# Calculate fold-change and then log2 fold-change
df_G532['fc'] = df_G532['sgG532_neg'].divide(df_G532['sgG532_uns'])
df_G532 = df_G532.assign(lfc=df_G532.apply(
    lambda x: np.log2(x.fc), axis=1))

# Prepare file and export
df_G532 = df_G532[['Allele','Allele_parsed','aa_seq','be_summary',
                   'sgG532_uns','sgG532_neg','lfc']]
df_G532.to_csv('G532_allele_freq.csv', index=False)



#%% Make barplot

# Put data in long form
df_W698_melt = df_W698.melt(id_vars=['Allele', 'be_summary','aa_seq'],
                            value_vars=['sgW698_uns','sgW698_pos'])
df_E756_melt = df_E756.melt(id_vars=['Allele', 'be_summary','aa_seq'],
                            value_vars=['sgE756_uns','sgE756_pos'])
df_G532_melt = df_G532.melt(id_vars=['Allele', 'be_summary','aa_seq'],
                            value_vars=['sgG532_uns','sgG532_neg'])


plot_bars(df_data=df_W698_melt, out_prefix=out_dir+'/sgW698')
plot_bars(df_data=df_E756_melt, out_prefix=out_dir+'/sgE756')
plot_bars(df_data=df_G532_melt, out_prefix=out_dir+'/sgG532')




