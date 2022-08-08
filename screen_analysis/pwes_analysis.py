#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Analyze screening results using proximity-weighted enrichment score
calculations and hierarchical clustering, including calculation of 'summed
delta PWES' for each sgRNA.

@author: Nicholas Lue and Kevin Ngan
"""

#%% Import Packages

import numpy as np
import pandas as pd

from pwes_functions import (
    get_centroids, get_pairwise_dist, chunks_scale, chunks_cluster
)


#%% Input parameters for script

# Structures (PDB IDs) and parameters to use for PWES calculations
struct_A = '4U7T' #Active DNMT3A ADD-MTase
struct_B = '4U7P' #Autoinhibited DNMT3A ADD-MTase
clust_num = 8 #Number of clusters
std, m, theta = (16,2,0.8) #(std, m, theta)

# Read in output from analysis/data processing script
df_comps = pd.read_csv('DNMT3A_scanning_sgRNA_scores.csv')


#%% Part 1: Read in screening data and prepare for calculations

# For PWES analysis, use d9 citrine-positive sgRNA scores (these are
# normalized to intergenic controls)
df_comps['d9pos_norm'] = df_comps['d9-pos'].copy() #Copy over

# Filter to retain only missense guides
df_missense = df_comps.loc[df_comps['Mut_type'] == 'Missense'].copy()

# Round edit sites to whole numbers (0.5 rounded to nearest even) and convert
# to long isoform numbering (note: 3A2 residues <25 don't map, but these
# will be filtered out later since no structure extends to this region).
df_missense['Clust_site'] = df_missense['Edit_site_3A2'].round(decimals=0)
df_missense['Clust_site'] = df_missense['Clust_site'].add(189)

# Get centroids
# 4U7T: chains A/C are both 3A. C has 1 less unresolved aa, though no guide maps to it.
df_centroids_A = get_centroids(pdb_id=struct_A, aa_sel='chain C and resi 476-912')
# 4U7P: chain A is 3A
df_centroids_B = get_centroids(pdb_id=struct_B, aa_sel='segi A')

# Calculate pairwise distances
df_pwdist_A = get_pairwise_dist(df_centroids_A)
df_pwdist_B = get_pairwise_dist(df_centroids_B)

# Filter guides to only those that map to residues resolved in
# both Structure A and Structure B
df_lfc_AB = df_missense.loc[df_missense['Clust_site'].isin(df_pwdist_A.index)].copy()
df_lfc_AB = df_lfc_AB.loc[df_lfc_AB['Clust_site'].isin(df_pwdist_B.index)]
df_lfc_AB = df_lfc_AB.sort_values(by=['Clust_site', 'sgRNA_ID']).reset_index()


#%% Part 2: Calculate PWES and perform clustering

# Calculate pairwise PWES using missense guides mapping to both structures
df_scaled_A = chunks_scale(df_lfc=df_lfc_AB, cond='d9pos_norm',
                           df_pwdist=df_pwdist_A, std=std, m=m, theta=theta)
df_scaled_B = chunks_scale(df_lfc=df_lfc_AB, cond='d9pos_norm',
                           df_pwdist=df_pwdist_B, std=std, m=m, theta=theta)

# Perform clustering
df_clusters_A, link_A = chunks_cluster(df_scaled=df_scaled_A, df_lfc=df_lfc_AB,
                                       num_clusters=clust_num)
df_clusters_B, link_B = chunks_cluster(df_scaled=df_scaled_B, df_lfc=df_lfc_AB,
                                       num_clusters=clust_num)

# Combine cluster IDs into one dataframe along with other data/annotations
df_clusters_A['4U7T_clusters'] = 'Cluster '
df_clusters_A['4U7T_clusters'] = df_clusters_A['4U7T_clusters'].add(
    df_clusters_A['Group'].astype(str))
df_lfc_AB = df_lfc_AB.merge(df_clusters_A[['sgRNA_ID', '4U7T_clusters']], on='sgRNA_ID')
df_clusters_B['4U7P_clusters'] = 'Cluster '
df_clusters_B['4U7P_clusters'] = df_clusters_B['4U7P_clusters'].add(
    df_clusters_B['Group'].astype(str))
df_lfc_AB = df_lfc_AB.merge(df_clusters_B[['sgRNA_ID', '4U7P_clusters']], on='sgRNA_ID')


#%% Calculate summed delta PWES
# First, find 'interaction scores' for each guide and compare between structures
# We define 'interaction score' as the sum of each guide's (absolute value)
# proximity-weighted enrichment scores against all other guides, not including
# itself. This is equivalent to the column sums of df_scaled.abs(), minus the
# diagonal entries. Then, to calculate 'summed delta PWES', take the difference
# between interaction scores for structure A vs. structure B.

# Initialize dataframe to hold info
df_sums = df_lfc_AB.copy()

# Calculate interaction scores based on struct A, append as column
df_abs_A = df_scaled_A.abs()
np.fill_diagonal(df_abs_A.values, 0) # Set diagonal to zero
df_sums_A = df_abs_A.sum(axis=0).to_frame(name='4U7T_scores').reset_index()
df_sums = df_sums.merge(df_sums_A, on='sgRNA_ID')

# Calculate interaction scores based on struct B, append as column
df_abs_B = df_scaled_B.abs()
np.fill_diagonal(df_abs_B.values, 0) # Set diagonal to zero
df_sums_B = df_abs_B.sum(axis=0).to_frame(name='4U7P_scores').reset_index()
df_sums = df_sums.merge(df_sums_B, on='sgRNA_ID')

# Calculate difference in interaction score between two structures
# This is the 'summed delta PWES'
df_sums['Summed_delta_PWES'] = df_sums['4U7T_scores'].sub(df_sums['4U7P_scores'])



#%% Export summary csv file

# Export calculations (i.e., for Supplementary Table)
df_sums.to_csv('DNMT3A_scanning_PWES_analysis_summary.csv', index=False)



