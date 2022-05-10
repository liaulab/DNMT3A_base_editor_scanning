#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Analysis of reporter bisulfite amplicon sequencing to produce plot shown in
Extended Data Fig. 1c.

@author: Nicholas Lue
"""

#%% Import Packages

import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

#%% Inputs

# File names (CRISPResso2 output files)
in_1 = 'amp1_nt_freq_summary.csv'
in_2 = 'amp2_nt_freq_summary.csv'

# Read in nt frequency summary data
df_nt1 = pd.read_csv(in_1, header=None, index_col=[0,1])
df_nt2 = pd.read_csv(in_2, header=None, index_col=[0,1])

# Re-number columns
df_nt1.columns = [i+1 for i in range(146)]
df_nt2.columns = [i+1 for i in range(146)]

# Sample names
samples1 = ['d0_1', 'd3_1', 'd6_1', 'd9_1', 'd12_1']
samples2 = ['d0_2', 'd3_2', 'd6_2', 'd9_2', 'd12_2']

#%% Function to quantify methylation
# Here, we consider % methylation to be C/(C+T)*100 or G/(G+A)*100,
# depending on which strand is sequenced (i.e., how primers were designed).
# We will calculate % methylation at every C (or G), regardless of whether it
# is in a CpG context or not. For the other bases, A/T/G (or A/T/C), this
# function will just return 0 since it doesn't make sense to define a
# methylation for that position.

def quant_methyl(df_in, list_cols, list_ignore, edit_tup=('C','T')):
    
    # Initialize dataframe for output
    df_out = pd.DataFrame(columns=['position','me_level'])
    
    for col in df_in.columns.tolist():
        
        # If column (position) is in list_ignore, return 0
        if col in list_ignore:
            df_out = df_out.append({'position':col,'me_level':0.0},
                                   ignore_index=True)
            continue
        
        # Otherwise, slice relevant column corresponding to base/CpG site of interest
        df_sel = df_in[col].astype(int)
        
        # Define base before/after base edit
        base_i, base_f = edit_tup
        
        # Calculate
        n_converted = df_sel.loc[base_f]
        n_unconverted = df_sel.loc[base_i]
        frac = n_unconverted / (n_converted + n_unconverted) * 100
        
        df_out = df_out.append({'position':col,'me_level':frac},
                               ignore_index=True)
    
    return df_out


#%% Amplicon 1

# CpG positions, 1-indexed
# For amplicon 1 these are the G in CpG, which turn to As
amp1 = 'AGGCCAGGCCTCAACTCAAGCACGAGGCGAAGGGGCTCCTTAAGCGCAAGGCCTCGAACTCTCCCACCCACTTCCAACCCGAAGCTCGGGATCAAGAATCACGTACTGCAGCCAGGTGGAAGTAATTCAAGGCACGCAAGGGCCAT'
pos1_CpG = [24, 29, 46, 56, 81, 88, 103, 136]
pos1_allG = [i+1 for i in range(len(amp1)) if amp1[i]=='G']
pos1_other = [i for i in range(1,len(amp1)+1) if i not in (pos1_CpG + pos1_allG)]

list_df1 = []

for sample in samples1:
    
    # CpGs
    df_samp = quant_methyl(df_in=df_nt1.loc[sample], list_cols=pos1_allG,
                           list_ignore=pos1_other, edit_tup=('G','A'))
    
    df_samp['sample'] = sample
    list_df1.append(df_samp)    

df_quants1 = pd.concat(list_df1)

# Correct positions so that they're along the opposite strand
df_quants1['corr_position'] = 147 - df_quants1['position']
# Add 253 so that position counting begins at pEF
df_quants1['corr_position'] = df_quants1['corr_position'] + 253


#%% Amplicon 2

# CpG positions, 1-indexed
# For amplicon 2 these are the C in CpG, which turn to Ts
amp2 = 'TCTGGCATGGTGGCAAGCTTGATATCGAATTCCTGCAGCCCGGGGGATCCGCCCGGGCTAGAGCGGCCGCCACCGCGGTGGAGCTCCAGCTTTTGTTCCCTTTAGTGAGGGTTAATTTCGAGCTTGGCGTAATCGCTAGCTCACGA'
pos2_CpG = [26, 41, 50, 54, 64, 68, 74, 76, 119, 128, 134, 144]
pos2_allC = [i+1 for i in range(len(amp2)) if amp2[i]=='C']
pos2_other = [i for i in range(1,len(amp2)+1) if i not in (pos2_CpG + pos2_allC)]

list_df2 = []

for sample in samples2:
    
    # CpGs
    df_samp = quant_methyl(df_in=df_nt2.loc[sample], list_cols=pos2_allC,
                           list_ignore=pos2_other, edit_tup=('C','T'))
    
    df_samp['sample'] = sample
    list_df2.append(df_samp)    
    
df_quants2 = pd.concat(list_df2)

# Correct positions so that they're along the opposite strand
df_quants2['corr_position'] = 147 - df_quants2['position']
# Add 1149 so that position counting begins at pEF
df_quants2['corr_position'] = df_quants2['corr_position'] + 1149


#%% Plotting function

# Define plotting function
def plot_methylation(df_input, x_col='position', y_col='me_level',
                     hue_col='sample', hue_order=None, palette=None,
                     list_CpG=None, 
                     out_title='plot.pdf', sns_context='paper'):
    '''
    Make lineplot showing methylation over time.

    '''
    # Plotting parameters and variables
    sns.set_context(sns_context)
    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['ps.fonttype'] = 42    
    mpl.rcParams['font.sans-serif'] = ['Arial']
    
    # Plot lineplot
    fig, ax = plt.subplots(figsize=(5,4))
    sns.lineplot(data=df_input, ax=ax, x=x_col, y=y_col, hue=hue_col,
                 hue_order=hue_order, palette=palette)
    # Plot markers only for CpG positions
    df_CpG = df_input.loc[df_input['position'].isin(list_CpG)]
    sns.scatterplot(data=df_CpG, ax=ax, x=x_col, y=y_col, hue=hue_col,
                    hue_order=hue_order, palette=palette, linewidth=1,
                    edgecolor='black', s=50, zorder=3)
    
    # Save figure
    plt.tight_layout()
    plt.savefig(out_title, format='pdf', bbox_inches='tight')
    plt.close()

#%% Calls

# Plot amplicon 1
rev_samplist1 = samples1.copy()
rev_samplist1.reverse()
plot_methylation(df_input=df_quants1, x_col='corr_position',
                 hue_order=rev_samplist1, palette='rocket',
                 list_CpG=pos1_CpG,
                 out_title='amp1_methylation.pdf')

# Plot amplicon 2
rev_samplist2 = samples2.copy()
rev_samplist2.reverse()
plot_methylation(df_input=df_quants2, x_col='corr_position',
                 hue_order=rev_samplist2, palette='rocket',
                 list_CpG=pos2_CpG,
                 out_title='amp2_methylation.pdf')



