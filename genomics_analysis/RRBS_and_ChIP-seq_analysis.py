#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Code for downstream analysis of RRBS and ChIP-seq data.

@author: Nicholas Lue
"""

#%% Import Packages

import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import wilcoxon


#%% Plotting functions

def plot_boxes(df_plot, x_col, y_col, out_name,
               hue_col=None, palette=sns.color_palette("Set2"), outliers=False,
               sns_context='paper', sns_palette='deep', width=0.8, size=(6,4)):
    
    # Plotting parameters and variables
    sns.set_context(sns_context)
    sns.set_palette(sns_palette)
    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['ps.fonttype'] = 42
    mpl.rcParams['font.sans-serif'] = ['Arial']
    
    # Plot
    fig, ax = plt.subplots(figsize=size)
    
    sns.boxplot(data=df_plot, x=x_col, y=y_col, width=width, palette=palette,
                saturation=1, showfliers=outliers, hue=hue_col, ax=ax)
    
    # Other adjustments
    plt.setp(ax.artists, edgecolor='black')
    plt.setp(ax.lines, color='black')
    plt.xticks(rotation=45, horizontalalignment='right')
    plt.tight_layout()
    
    # Save to pdf
    plt.savefig(out_name, format='pdf', bbox_inches='tight')
    plt.close()


def plot_scatter(df_plot, x_col, y_col, hue_col, out_name, palette='coolwarm',
                 s=20, hue_norm=None, alpha=0.75, xlim=None, ylim=None, legend='auto',
                 sns_context='paper', sns_palette='deep', size=(6,6)):
    
    # Plotting parameters and variables
    sns.set_context(sns_context)
    sns.set_palette(sns_palette)
    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['ps.fonttype'] = 42
    mpl.rcParams['font.sans-serif'] = ['Arial']
    
    # Plot
    fig, ax = plt.subplots(figsize=size)
    sns.scatterplot(data=df_plot, x=x_col, y=y_col, hue=hue_col, alpha=alpha,
                    legend=legend, hue_norm=hue_norm, linewidth=0, s=s,
                    palette=palette, ax=ax)
    
    # Adjust x and y axis limits
    if xlim is not None:
        plt.xlim(xlim[0],xlim[1])
    if ylim is not None:
        plt.ylim(ylim[0],ylim[1])
    
    # Other adjustments
    plt.tight_layout()
    
    # Save to pdf
    plt.savefig(out_name, format='pdf', bbox_inches='tight')
    plt.close()



#%% CpG-level analysis
# Code to make Fig. 5f, left.

df_ABCDE = pd.read_csv('ABCDE_global_methylation.tsv', sep='\t')

# Define sample names
rrbs_list = ['WT_1', 'WT_2', 'R297W_1', 'R297W_2', 'W326R_1', 'W326R_2',
             'E338K_1', 'E338K_2', 'E752K_1', 'E752K_2']
avg_list = ['WT', 'R297W', 'W326R', 'E338K', 'E752K']
col_list = ['chr', 'start', 'end'] + rrbs_list

# Rename columns
df_ABCDE.columns = col_list

# Perform additional processing
chr_list = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7',
            'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14',
            'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chrX', 'chrY']
# Retain only CpGs within canonical chromosomes
df_ABCDE = df_ABCDE.loc[df_ABCDE['chr'].isin(chr_list)]
# Drop CpGs with NA values, if any
df_ABCDE = df_ABCDE.dropna(axis=0)
# Average methylation of replicates
df_ABCDE['WT'] = df_ABCDE[['WT_1','WT_2']].mean(axis=1)
df_ABCDE['R297W'] = df_ABCDE[['R297W_1','R297W_2']].mean(axis=1)
df_ABCDE['W326R'] = df_ABCDE[['W326R_1','W326R_2']].mean(axis=1)
df_ABCDE['E338K'] = df_ABCDE[['E338K_1','E338K_2']].mean(axis=1)
df_ABCDE['E752K'] = df_ABCDE[['E752K_1','E752K_2']].mean(axis=1)

# Remove any CpGs that are 0% methylated across all samples
df_ABCDE_nonzero = df_ABCDE.loc[~(df_ABCDE[rrbs_list]==0).all(axis=1)].copy()

# Transform data to long form
df_ABCDE_nonzero_melt = df_ABCDE_nonzero.melt(id_vars=['chr','start','end'],
                                              value_vars=rrbs_list+avg_list).copy()

# Make boxplot of replicate-averaged methylation
plot_boxes(df_plot=df_ABCDE_nonzero_melt.loc[df_ABCDE_nonzero_melt['variable'].isin(avg_list)],
           x_col='variable', y_col='value', size=(2,4),
           out_name='global_nonzero_CpG_ABCDE_boxes_resized.pdf')



#%% Analysis of methylation in 500 kb bins
# Code to make Fig. 5f, right.

df_500 = pd.read_csv('ABCDE_intersection_500kb_bins_global_methylation.tsv', sep='\t')

# Rename columns
df_500.columns = col_list

# Perform additional processing
# Retain only bins within canonical chromosomes
df_500 = df_500.loc[df_500['chr'].isin(chr_list)]
# Drop bins with NA values, if any
df_500 = df_500.dropna(axis=0)
# Average methylation of replicates
df_500['WT'] = df_500[['WT_1','WT_2']].mean(axis=1)
df_500['R297W'] = df_500[['R297W_1','R297W_2']].mean(axis=1)
df_500['W326R'] = df_500[['W326R_1','W326R_2']].mean(axis=1)
df_500['E338K'] = df_500[['E338K_1','E338K_2']].mean(axis=1)
df_500['E752K'] = df_500[['E752K_1','E752K_2']].mean(axis=1)

# Transform data to long form
df_500_melt = df_500.melt(id_vars=['chr','start','end'],
                          value_vars=rrbs_list+avg_list).copy()

# Make boxplot of replicate-averaged methylation in 500 kb bins
plot_boxes(df_plot=df_500_melt.loc[df_500_melt['variable'].isin(avg_list)],
           x_col='variable', y_col='value', size=(2,4),
           out_name='global_500kb_ABCDE_boxes_resized.pdf')




#%% Histone modification quartile analysis
# Code to make Fig. 5g and Extended Data Fig. 9c

df_10_ABCDE = pd.read_csv('ABCDE_intersection_global_methylation.tsv', sep='\t')

# Rename columns
df_10_ABCDE.columns = col_list + ['H3K4me3', 'H3K27me3', 'H3K36me2']

# Perform additional processing
# Retain only CpGs within canonical chromosomes
df_10_ABCDE = df_10_ABCDE.loc[df_10_ABCDE['chr'].isin(chr_list)]
# Drop any bins with NA values (i.e., no covered CpGs fall within bin)
df_10_ABCDE = df_10_ABCDE.dropna(axis=0)
# Average methylation of replicates
df_10_ABCDE['WT'] = df_10_ABCDE[['WT_1','WT_2']].mean(axis=1)
df_10_ABCDE['R297W'] = df_10_ABCDE[['R297W_1','R297W_2']].mean(axis=1)
df_10_ABCDE['W326R'] = df_10_ABCDE[['W326R_1','W326R_2']].mean(axis=1)
df_10_ABCDE['E338K'] = df_10_ABCDE[['E338K_1','E338K_2']].mean(axis=1)
df_10_ABCDE['E752K'] = df_10_ABCDE[['E752K_1','E752K_2']].mean(axis=1)

# Quartile analysis
df_quartiles_ABCDE = df_10_ABCDE.copy()
df_quartiles_ABCDE['K36_quartiles'] = pd.qcut(df_quartiles_ABCDE['H3K36me2'],q=4,labels=False)
df_quartiles_ABCDE['K27_quartiles'] = pd.qcut(df_quartiles_ABCDE['H3K27me3'],q=4,labels=False)
df_quartiles_ABCDE['K4_quartiles'] = pd.qcut(df_quartiles_ABCDE['H3K4me3'],q=4,labels=False)

# Transform data to long form
df_quartiles_ABCDE_melt = df_quartiles_ABCDE.melt(
    id_vars=['chr','start','end','K36_quartiles','K27_quartiles','K4_quartiles'],
    value_vars=rrbs_list+avg_list+['H3K4me3','H3K27me3','H3K36me2']).copy()

# Plots
plot_boxes(df_plot=df_quartiles_ABCDE_melt.loc[df_quartiles_ABCDE_melt['variable'].isin(avg_list)],
           x_col='variable', y_col='value', hue_col='K36_quartiles',
           out_name='K36_smoothed_quartile_methylation_10kb_bins.pdf')
plot_boxes(df_plot=df_quartiles_ABCDE_melt.loc[df_quartiles_ABCDE_melt['variable'].isin(avg_list)],
           x_col='variable', y_col='value', hue_col='K4_quartiles',
           out_name='K4_smoothed_quartile_methylation_10kb_bins.pdf')


#%% High-low quartile methylation difference analysis
# Analysis for Extended Data Fig. 9d

col_list_2 = rrbs_list+['K27_quartiles','K36_quartiles']
df_quartiles_collapsed = df_quartiles_ABCDE[col_list_2].groupby(by=['K36_quartiles']).median().reset_index()

# Initialize dataframe
df_diff = pd.DataFrame(columns=['variable','value'])

for samp in rrbs_list:
    diff = df_quartiles_collapsed.loc[3,samp] - df_quartiles_collapsed.loc[0,samp]
    df_diff = pd.concat([df_diff, pd.DataFrame([[samp[:-2],diff]],columns=['variable','value'])])

df_diff.to_csv('4th-1st_K36_quartile_delta_methylation.csv',index=False)

#%% Histone modification 1000-quantile analysis
# Code to make Extended Data Fig. 9e

# 1000 quantiles
df_fine = df_10_ABCDE[rrbs_list+avg_list+['H3K4me3', 'H3K27me3', 'H3K36me2']].copy()
df_fine['K36_quantiles'] = pd.qcut(df_fine['H3K36me2'],q=1000,labels=False)
df_fine['K27_quantiles'] = pd.qcut(df_fine['H3K27me3'],q=1000,labels=False)
df_fine['K4_quantiles'] = pd.qcut(df_fine['H3K4me3'],q=1000,labels=False)

df_K36_fine = df_fine.groupby(by=['K36_quantiles']).mean().reset_index()
df_K36_fine = df_K36_fine.melt(
    id_vars=['H3K4me3','H3K27me3','H3K36me2','K36_quantiles','K27_quantiles','K4_quantiles'],
    value_vars=rrbs_list+avg_list)
df_K27_fine = df_fine.groupby(by=['K27_quantiles']).mean().reset_index()
df_K27_fine = df_K27_fine.melt(
    id_vars=['H3K4me3','H3K27me3','H3K36me2','K36_quantiles','K27_quantiles','K4_quantiles'],
    value_vars=rrbs_list+avg_list)
df_K4_fine = df_fine.groupby(by=['K4_quantiles']).mean().reset_index()
df_K4_fine = df_K4_fine.melt(
    id_vars=['H3K4me3','H3K27me3','H3K36me2','K36_quantiles','K27_quantiles','K4_quantiles'],
    value_vars=rrbs_list+avg_list)

# Make scatterplots for K36 but only 1 variant at a time
plot_scatter(df_plot=df_K36_fine.loc[df_K36_fine['variable']=='WT'],
             x_col='H3K36me2', y_col='value', ylim=(0,0.16), xlim=(-3,4),
             hue_col='variable', palette=None, legend=False,
             out_name='1000_K36_quantiles_CpGme_scatter_WT_only.pdf')
plot_scatter(df_plot=df_K36_fine.loc[df_K36_fine['variable']=='R297W'],
             x_col='H3K36me2', y_col='value', ylim=(0,0.16), xlim=(-3,4),
             hue_col='variable', palette=None, legend=False,
             out_name='1000_K36_quantiles_CpGme_scatter_R297W_only.pdf')
plot_scatter(df_plot=df_K36_fine.loc[df_K36_fine['variable']=='W326R'],
             x_col='H3K36me2', y_col='value', ylim=(0,0.16), xlim=(-3,4),
             hue_col='variable', palette=None, legend=False,
             out_name='1000_K36_quantiles_CpGme_scatter_W326R_only.pdf')
plot_scatter(df_plot=df_K36_fine.loc[df_K36_fine['variable']=='E338K'],
             x_col='H3K36me2', y_col='value', ylim=(0,0.16), xlim=(-3,4),
             hue_col='variable', palette=None, legend=False,
             out_name='1000_K36_quantiles_CpGme_scatter_E338K_only.pdf')
plot_scatter(df_plot=df_K36_fine.loc[df_K36_fine['variable']=='E752K'],
             x_col='H3K36me2', y_col='value', ylim=(0,0.16), xlim=(-3,4),
             hue_col='variable', palette=None, legend=False,
             out_name='1000_K36_quantiles_CpGme_scatter_E752K_only.pdf')




