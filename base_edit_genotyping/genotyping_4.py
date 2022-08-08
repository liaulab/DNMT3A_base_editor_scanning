#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Analysis of sequencing data for clones of sgE756-, sgG532-, sgR301-, and sgE907-
treated cells presented in Extended Data Fig. 3.

@author: Nicholas Lue
"""

#%% Import Packages

import os
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
list_28 = ['NL'+str(i).zfill(2) for i in range(1,25)] #samples for sgE756
list_33 = ['NL'+str(i).zfill(2) for i in range(25,49)] #samples for sgG532
list_128 = ['NL'+str(i).zfill(2) for i in range(49,73)] #samples for sgR301
list_651 = ['NL'+str(i).zfill(2) for i in range(73,97)] #samples for sgE907
list_all = list_28 + list_33 + list_128 + list_651

# Read in allele tables
# Initialize dictionaries to hold sample-specific dataframes
allele_tables_dict = {}

for sample in list_all:
    
    # Define directory with sample-specific CRISPResso output
    sample_dir = os.path.join(data_dir, '_'.join(['CRISPResso_on', sample]))
    
    # Identify which guide this sample corresponds to
    guide_id = '' #reset guide_id
    if sample in list_28:
        guide_id = 'sgNL28'
    elif sample in list_33:
        guide_id = 'sgNL33'
    elif sample in list_128:
        guide_id = 'sg3A.128'
    elif sample in list_651:
        guide_id = 'sg3A.651'
    
    # Define the file name of the allele frequency table
    file = os.path.join(
        sample_dir, ('Alleles_frequency_table_around_sgRNA_' + guide_id + '.txt'))
    
    # Read .txt file and store in dictionary of data frames
    allele_tables_dict[sample] = pd.read_csv(file, sep='\t')

# Import batch settings
batch_settings = pd.read_csv('Batch_settings_4.txt', sep='\t', index_col=0)

#%% Define functions for analysis

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


# For clonal analysis — name allele if it is one of the top alleles as defined
# by dictionary. If allele is <5% frequency OR is not in top alleles list,
# classify as Other.
def classify_allele(allele, top_dict, freq):
    if freq < 5:
        return 'Other'
    elif allele in top_dict.keys():
        return top_dict[allele]
    else:
        return 'Other'
    return

#%% Plotting functions

def plot_stacked_bar(df_data, list_alleles, color_dict, out_prefix,
                     legend_check=True):

    # Plotting parameters and variables
    sns.set_context('talk')
    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['ps.fonttype'] = 42    
    mpl.rcParams['font.sans-serif'] = ['Arial']
    
    # Initialize figure
    fig, ax = plt.subplots(figsize=(8,6))
    
    # For legend
    list_bars = []
    
    # Iterate through list of allele possibilities, and add bars
    for allele in list_alleles:
        sns.barplot(data=df_data, x='Sample_ID', y=allele, color=color_dict[allele],
                    ax=ax)
        list_bars.append(mpl.patches.Patch(color=color_dict[allele], label=allele))
    
    # Make dashed lines corresponding to theoretical allele frequencies
    plt.axhline(y=100/3, color='k', ls='--', lw=1)
    plt.axhline(y=100/3*2, color='k', ls='--', lw=1)
    plt.axhline(y=100, color='k', ls='--', lw=1)
    
    # Rotate xticks
    plt.xticks(rotation=45, horizontalalignment='right')
    
    # Final adjustments, then save figure
    plt.ylabel('% aligned reads')
    if legend_check:
        plt.legend(handles=list_bars)
    plt.tight_layout()
    plt.savefig(out_prefix+'_allele_barplot.pdf', format='pdf', bbox_inches='tight')
    plt.close()





#%% Analysis of sgNL28 (sgE756) clones

# (Prep work) Get a list of the top alleles
top_28 = []
for sample in list_28:
    temp = allele_tables_dict[sample].copy()
    temp = temp.loc[temp['%Reads']>=5]
    list_temp = temp['Aligned_Sequence'].tolist()
    top_28.extend([allele for allele in list_temp if allele not in top_28])

# Make dataframe with info on top alleles
df_top_28 = pd.DataFrame(top_28, columns=['Allele'])
# Parse and translate
df_top_28 = df_top_28.assign(Allele_parsed=df_top_28.apply(
    lambda x: parse_allele(x.Allele[1:], revcomp=True), axis=1))
df_top_28 = df_top_28.assign(Translation=df_top_28.apply(
    lambda x: translate_allele(x.Allele_parsed), axis=1))

# Make a dictionary assigning each top allele to a name
dict_top_28 = {
    'TGGCCACCACATTTTTAAAGAGCCAGAAGAAGGGGCGATC':'E756K',
    'TGGCCACCACATTCTCAAAGAGCCAGAAGAAGGGGCGATC':'WT',
    'TGGCCACCATATTTTTAAAGAGCCAGAAGAAGGGGCGATC':'E756K/V758M',
    'TGGCCACCACATTCTTAAAGAGCCAGAAGAAGGGGCGATC':'E756K'}

# Now, assemble final master dataframe with percent each allele for each clone
list_alleles_28 = ['WT','E756K','E756K/V758M','Other']
df_clones_28 = pd.DataFrame(columns=['Sample_ID']+list_alleles_28)

# Iterate through clones
for sample in list_28:
    temp = allele_tables_dict[sample].copy()
    temp['Percent_reads'] = temp['%Reads'] #Rename column
    
    # Classify/name each allele
    temp = temp.assign(Allele_classified=temp.apply(
        lambda x: classify_allele(x.Aligned_Sequence,dict_top_28,x.Percent_reads), axis=1))
    
    # Merge Other alleles (and repeats)
    temp_merged = temp.groupby('Allele_classified')['Percent_reads'].sum().to_frame().transpose()
    
    # Add sample ID
    temp_merged['Sample_ID'] = sample
    
    # Add information for this clone to the master dataframe
    df_clones_28 = pd.concat([df_clones_28,temp_merged], ignore_index=True)

# Change NaN values to 0
df_clones_28 = df_clones_28.fillna(0)

# Prep dataframe for plotting
# Sort dataframe
df_clones_28 = df_clones_28.sort_values(by=list_alleles_28,
                                        ascending=[False,False,False,False],
                                        ignore_index=True)
# Make sample ID the index
df_clones_28 = df_clones_28.set_index(keys='Sample_ID', verify_integrity=True)

# Calculate cumulative sum of allele frequences. This is necessary to make
# stacked barplot
df_clones_28_summed = df_clones_28[list_alleles_28].cumsum(axis=1)
df_clones_28_summed = df_clones_28_summed.reset_index()

# Make plot
color_dict = {'WT':'royalblue',
              'E756K':'coral',
              'E756K/V758M':'firebrick',
              'Other':'lightgray'}
plot_stacked_bar(df_data=df_clones_28_summed,
                 list_alleles=list_alleles_28[::-1],
                 color_dict=color_dict,
                 out_prefix=out_dir+'/sgE756_stacked_simplified')

#%% Analysis of sgNL33 (sgG532) clones

# (Prep work) Get a list of the top alleles
top_33 = []
for sample in list_33:
    temp = allele_tables_dict[sample].copy()
    temp = temp.loc[temp['%Reads']>=5]
    list_temp = temp['Aligned_Sequence'].tolist()
    top_33.extend([allele for allele in list_temp if allele not in top_33])

# Make dataframe with info on top alleles
df_top_33 = pd.DataFrame(top_33, columns=['Allele'])
# Parse and translate
df_top_33 = df_top_33.assign(Allele_parsed=df_top_33.apply(
    lambda x: parse_allele(x.Allele[2:-2], revcomp=True), axis=1))
df_top_33 = df_top_33.assign(Translation=df_top_33.apply(
    lambda x: translate_allele(x.Allele_parsed), axis=1))

# Make a dictionary assigning each top allele to a name
dict_top_33 = {
    'CAGTAGGACTGGTAGCCGTCGTCGTCGTACTGGTACGCAC':'WT',
    'CAGTAGGACTGGTAGTTGTCGTCGTCGTACTGGTACGCAC':'G532N/S/D',
    'CAGTAGGACTGGTAGCTGTCGTCGTCGTACTGGTACGCAC':'G532N/S/D',
    'CAGTAGAACTGGTAGTCGTCGTCGTCGTACTGGTACGCAC':'G532N/D_S535F',
    'CAGTAGAACTGGTAGTTGTCGTCGTCGTACTGGTACGCAC':'G532N/D_S535F',
    'CAGTAGGACTGGTAGTCGTCGTCGTCGTACTGGTACGCAC':'G532N/S/D'}

# Now, assemble final master dataframe with percent each allele for each clone
list_alleles_33 = ['WT','G532N/S/D','G532N/D_S535F','Other']
df_clones_33 = pd.DataFrame(columns=['Sample_ID']+list_alleles_33)

# Iterate through clones
for sample in list_33:
    temp = allele_tables_dict[sample].copy()
    temp['Percent_reads'] = temp['%Reads'] #Rename column
    
    # Classify/name each allele
    temp = temp.assign(Allele_classified=temp.apply(
        lambda x: classify_allele(x.Aligned_Sequence,dict_top_33,x.Percent_reads), axis=1))
    
    # Merge Other alleles (and repeats)
    temp_merged = temp.groupby('Allele_classified')['Percent_reads'].sum().to_frame().transpose()
    
    # Add sample ID
    temp_merged['Sample_ID'] = sample
    
    # Add information for this clone to the master dataframe
    df_clones_33 = pd.concat([df_clones_33,temp_merged], ignore_index=True)

# Change NaN values to 0
df_clones_33 = df_clones_33.fillna(0)

# Prep dataframe for plotting
# Sort dataframe
df_clones_33 = df_clones_33.sort_values(by=list_alleles_33,
                                        ascending=[False,False,False,False],
                                        ignore_index=True)
# Make sample ID the index
df_clones_33 = df_clones_33.set_index(keys='Sample_ID', verify_integrity=True)

# Calculate cumulative sum of allele frequences. This is necessary to make
# stacked barplot
df_clones_33_summed = df_clones_33[list_alleles_33].cumsum(axis=1)
df_clones_33_summed = df_clones_33_summed.reset_index()

# Make plot
color_dict = {'WT':'royalblue',
              'G532N/S/D':'coral',
              'G532N/D_S535F':'firebrick',
              'Other':'lightgray'}
plot_stacked_bar(df_data=df_clones_33_summed,
                 list_alleles=list_alleles_33[::-1],
                 color_dict=color_dict,
                 out_prefix=out_dir+'/sgG532_stacked_simplified')

#%% Analysis of sg3A.128 (sgR301) clones

# (Prep work) Get a list of the top alleles
top_128 = []
for sample in list_128:
    temp = allele_tables_dict[sample].copy()
    temp = temp.loc[temp['%Reads']>=5]
    list_temp = temp['Aligned_Sequence'].tolist()
    top_128.extend([allele for allele in list_temp if allele not in top_128])

# Make dataframe with info on top alleles
df_top_128 = pd.DataFrame(top_128, columns=['Allele'])
# Parse and translate
df_top_128 = df_top_128.assign(Allele_parsed=df_top_128.apply(
    lambda x: parse_allele(x.Allele[:-1], revcomp=False), axis=1))
df_top_128 = df_top_128.assign(Translation=df_top_128.apply(
    lambda x: translate_allele(x.Allele_parsed), axis=1))

# Make a dictionary assigning each top allele to a name
dict_top_128 = {
    'GTGTGGGGGAAATTGTGGGGCTTTTTCTGGTGGCCAGGCC':'R301W/S304F',
    'GTGTGGGGGAAACTGCGGGGCTTTTCCTGGTGGCCAGGCC':'silent',
    'GTGTGGGGGAAATTGTGGGGCTTCTCCTGGTGGCCAGGCC':'R301W',
    'GTGTGGGGGAAACTGCGGGGCTTCTCCTGGTGGCCAGGCC':'WT',
    'GTGTGGGGGAAACTGTGGGGCTTCTCCTGGTGGCCAGGCC':'R301W',
    'GTGTGGGGGAAACTGCGGGGCTTTTTCTGGTGGCCAGGCC':'S304F',
    'GTGTGGGGGAAAATGTGGGGCTTCTCCTGGTGGCCAGGCC':'L300M/R301W',
    'GTGTGGGGGAAATTGCGGGGCTTCTCCTGGTGGCCAGGCC':'silent',
    'GTGTGGGGGAAACTGCGGGGCTTCTTCTGGTGGCCAGGCC':'S304F',
    'GTGTGGGGGAAACTGTGGGGCTTCTTCTGGTGGCCAGGCC':'R301W/S304F',
    'GTGTGGGGGAAATTGCGGGGCTTCTTCTGGTGGCCAGGCC':'S304F',
    'GTGTGGGGGAAATTGTGGGGCTTCTTCTGGTGGCCAGGCC':'R301W/S304F'}

# Now, assemble final master dataframe with percent each allele for each clone
list_alleles_128 = ['WT','silent','R301W','S304F','R301W/S304F',
                    'L300M/R301W','Other']
df_clones_128 = pd.DataFrame(columns=['Sample_ID']+list_alleles_128)

# Iterate through clones
for sample in list_128:
    temp = allele_tables_dict[sample].copy()
    temp['Percent_reads'] = temp['%Reads'] #Rename column
    
    # Classify/name each allele
    temp = temp.assign(Allele_classified=temp.apply(
        lambda x: classify_allele(x.Aligned_Sequence,dict_top_128,x.Percent_reads), axis=1))
    
    # Merge Other alleles (and repeats)
    temp_merged = temp.groupby('Allele_classified')['Percent_reads'].sum().to_frame().transpose()
    
    # Add sample ID
    temp_merged['Sample_ID'] = sample
    
    # Add information for this clone to the master dataframe
    df_clones_128 = pd.concat([df_clones_128,temp_merged], ignore_index=True)

# Change NaN values to 0
df_clones_128 = df_clones_128.fillna(0)

# Prep dataframe for plotting
# Sort dataframe
# For purposes of sorting, group WT and silent alleles
df_clones_128['WT_silent'] = df_clones_128['WT'].add(df_clones_128['silent'])
df_clones_128 = df_clones_128.sort_values(by=['WT_silent']+list_alleles_128,
                                        ascending=[False,False,False,False,False,False,False,False],
                                        ignore_index=True)
# Make sample ID the index
df_clones_128 = df_clones_128.set_index(keys='Sample_ID', verify_integrity=True)

# Calculate cumulative sum of allele frequences. This is necessary to make
# stacked barplot
df_clones_128_summed = df_clones_128[list_alleles_128].cumsum(axis=1)
df_clones_128_summed = df_clones_128_summed.reset_index()

# Make plot
color_dict = {'WT':'royalblue',
              'silent':'cornflowerblue',
              'R301W':'coral',
              'S304F':'gold',
              'R301W/S304F':'orange',
              'L300M/R301W':'firebrick',
              'Other':'lightgray'}
plot_stacked_bar(df_data=df_clones_128_summed,
                 list_alleles=list_alleles_128[::-1],
                 color_dict=color_dict,
                 out_prefix=out_dir+'/sgR301_stacked_simplified')

#%% Analysis of sg3A.651 (sgE907) clones

# (Prep work) Get a list of the top alleles
top_651 = []
for sample in list_651:
    temp = allele_tables_dict[sample].copy()
    temp = temp.loc[temp['%Reads']>=5]
    list_temp = temp['Aligned_Sequence'].tolist()
    top_651.extend([allele for allele in list_temp if allele not in top_651])

# Make dataframe with info on top alleles
df_top_651 = pd.DataFrame(top_651, columns=['Allele'])
# Parse and translate
df_top_651 = df_top_651.assign(Allele_parsed=df_top_651.apply(
    lambda x: parse_allele(x.Allele[2:-2], revcomp=True), axis=1))
df_top_651 = df_top_651.assign(Translation=df_top_651.apply(
    lambda x: translate_allele(x.Allele_parsed), axis=1))

# Make a dictionary assigning each top allele to a name
dict_top_651 = {
    'CACGCAAAATACTTTTTTAGCGGAGCGAAGAGGTGGCGGA':'E907K',
    'CACGCAAAATACTCCTTCAGCGGAGCGAAGAGGTGGCGGA':'WT',
    'CACGCAAAATACTTTTTCAGCGGAGCGAAGAGGTGGCGGA':'E907K',
    'CACGCAAAATACTCCTTAAGCGGAGCGAAGAGGTGGCGGA':'silent',
    'CACGCAAAATACTTTTT-AGCGGAGCGAAGAGGTGGCGGA':'frameshift'}

# Now, assemble final master dataframe with percent each allele for each clone
list_alleles_651 = ['WT','silent','E907K','frameshift','Other']
df_clones_651 = pd.DataFrame(columns=['Sample_ID']+list_alleles_651)

# Iterate through clones
for sample in list_651:
    temp = allele_tables_dict[sample].copy()
    temp['Percent_reads'] = temp['%Reads'] #Rename column
    
    # Classify/name each allele
    temp = temp.assign(Allele_classified=temp.apply(
        lambda x: classify_allele(x.Aligned_Sequence,dict_top_651,x.Percent_reads), axis=1))
    
    # Merge Other alleles (and repeats)
    temp_merged = temp.groupby('Allele_classified')['Percent_reads'].sum().to_frame().transpose()
    
    # Add sample ID
    temp_merged['Sample_ID'] = sample
    
    # Add information for this clone to the master dataframe
    df_clones_651 = pd.concat([df_clones_651,temp_merged], ignore_index=True)

# Change NaN values to 0
df_clones_651 = df_clones_651.fillna(0)

# Prep dataframe for plotting
# Sort dataframe
# For purposes of sorting, group WT and silent alleles
df_clones_651['WT_silent'] = df_clones_651['WT'].add(df_clones_651['silent'])
df_clones_651 = df_clones_651.sort_values(by=['WT_silent']+list_alleles_651,
                                        ascending=[False,False,False,False,False,False],
                                        ignore_index=True)
# Make sample ID the index
df_clones_651 = df_clones_651.set_index(keys='Sample_ID', verify_integrity=True)

# Calculate cumulative sum of allele frequences. This is necessary to make
# stacked barplot
df_clones_651_summed = df_clones_651[list_alleles_651].cumsum(axis=1)
df_clones_651_summed = df_clones_651_summed.reset_index()

# Make plot
color_dict = {'WT':'royalblue',
              'silent':'cornflowerblue',
              'E907K':'coral',
              'frameshift':'firebrick',
              'Other':'lightgray'}
plot_stacked_bar(df_data=df_clones_651_summed,
                 list_alleles=list_alleles_651[::-1],
                 color_dict=color_dict,
                 out_prefix=out_dir+'/sgE907_stacked_simplified')

#%% Make supplementary table

# Initialize table
final_cols = ['Sample_ID','Allele','Translation','Annotation','Percent_reads']
df_supp = pd.DataFrame(columns=final_cols)

# sgNL28
list_28_reordered = df_clones_28_summed['Sample_ID'].tolist()
for sample in list_28_reordered:
    temp = allele_tables_dict[sample].copy()
    temp = temp.loc[temp['%Reads']>=5]
    temp['Sample_ID'] = sample
    temp = temp.rename(mapper={'Aligned_Sequence':'Allele','%Reads':'Percent_reads'},
                       axis=1)
    temp = temp.assign(Allele_parsed=temp.apply(
        lambda x: parse_allele(x.Allele[1:], revcomp=True), axis=1))
    temp = temp.assign(Translation=temp.apply(
        lambda x: translate_allele(x.Allele_parsed), axis=1))
    temp = temp.assign(Annotation=temp.apply(
        lambda x: classify_allele(x.Allele,dict_top_28,x.Percent_reads), axis=1))
    df_supp = pd.concat([df_supp,temp[final_cols]], ignore_index=True)

# sgNL33
list_33_reordered = df_clones_33_summed['Sample_ID'].tolist()
for sample in list_33_reordered:
    temp = allele_tables_dict[sample].copy()
    temp = temp.loc[temp['%Reads']>=5]
    temp['Sample_ID'] = sample
    temp = temp.rename(mapper={'Aligned_Sequence':'Allele','%Reads':'Percent_reads'},
                       axis=1)
    temp = temp.assign(Allele_parsed=temp.apply(
        lambda x: parse_allele(x.Allele[2:-2], revcomp=True), axis=1))
    temp = temp.assign(Translation=temp.apply(
        lambda x: translate_allele(x.Allele_parsed), axis=1))
    temp = temp.assign(Annotation=temp.apply(
        lambda x: classify_allele(x.Allele,dict_top_33,x.Percent_reads), axis=1))
    df_supp = pd.concat([df_supp,temp[final_cols]], ignore_index=True)

# sg3A.128
list_128_reordered = df_clones_128_summed['Sample_ID'].tolist()
for sample in list_128_reordered:
    temp = allele_tables_dict[sample].copy()
    temp = temp.loc[temp['%Reads']>=5]
    temp['Sample_ID'] = sample
    temp = temp.rename(mapper={'Aligned_Sequence':'Allele','%Reads':'Percent_reads'},
                       axis=1)
    temp = temp.assign(Allele_parsed=temp.apply(
        lambda x: parse_allele(x.Allele[:-1], revcomp=False), axis=1))
    temp = temp.assign(Translation=temp.apply(
        lambda x: translate_allele(x.Allele_parsed), axis=1))
    temp = temp.assign(Annotation=temp.apply(
        lambda x: classify_allele(x.Allele,dict_top_128,x.Percent_reads), axis=1))
    df_supp = pd.concat([df_supp,temp[final_cols]], ignore_index=True)

# sg3A.651
list_651_reordered = df_clones_651_summed['Sample_ID'].tolist()
for sample in list_651_reordered:
    temp = allele_tables_dict[sample].copy()
    temp = temp.loc[temp['%Reads']>=5]
    temp['Sample_ID'] = sample
    temp = temp.rename(mapper={'Aligned_Sequence':'Allele','%Reads':'Percent_reads'},
                       axis=1)
    temp = temp.assign(Allele_parsed=temp.apply(
        lambda x: parse_allele(x.Allele[2:-2], revcomp=True), axis=1))
    temp = temp.assign(Translation=temp.apply(
        lambda x: translate_allele(x.Allele_parsed), axis=1))
    temp = temp.assign(Annotation=temp.apply(
        lambda x: classify_allele(x.Allele,dict_top_651,x.Percent_reads), axis=1))
    df_supp = pd.concat([df_supp,temp[final_cols]], ignore_index=True)

# Export final table
df_supp.to_csv(out_dir+'/Clonal_allele_table.csv', index=False)




