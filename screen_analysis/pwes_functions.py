#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Functions for PWES analysis.

@author: Kevin Ngan and Nicholas Lue
"""

#%% Import Packages

import time
import numpy as np
import pandas as pd
from pymol import cmd, stored
import scipy.spatial.distance as sp_sd
import scipy.cluster.hierarchy as sp_hca


#%% get_centroids()- get centroid xyz coords from structures w/ pymol

def get_centroids(pdb_id, aa_sel=None, list_aas=None, save_centroids=True,
                  out_csv=None, return_df=True):
    """
    This function gets all centroids from a pymol structure and outputs as csv

    Parameters
    ----------
    pdb_id : str
        pdb id of the desired structure
    aa_sel : str
        str in pymol syntax for the desired aa selection (NO PARENTHESES)
        The default is pdb_id + chain A + polymer.protein
    list_aas : list
        list of all amino acids to find centroids for.
        By default uses first to last AA in aa_sel w/ resolved alpha carbons
    save_centroids : TYPE, optional
        choose whether to save the centroid df as csv. The default is True.
    out_csv : str
        name of the output csv file. The default is pdb_id + _centroids.csv
    return_df : bool
        choose whether or not to return the centroid df (default is True)
    """

    # fetch pdb structure
    cmd.fetch(pdb_id)
    # default selection for centroids is the entire crystal structure
    # if a custom selection is defined, MUST be in pymol syntax
    if aa_sel is None:
        aa_sel = '%s and chain A and polymer.protein' % pdb_id
    if list_aas is None:
        # make list to hold all amino acids in the aa_sel
        stored.aa=[]
        # iterate thru AA (alpha Cs), append (resi num, resi name) to stored.list
        # pymol syntax iterate (selection), expression
        cmd.iterate('(' + aa_sel + ' and name ca)', 'stored.aa.append((int(resi),resn))')
        list_aas = [x for x in range(min(stored.aa)[0], max(stored.aa)[0] + 1)]

    # iterate through residues and add centroids to df as (x,y,z)
    # determine centroid by average of atoms(x,y,z) coords
    centroid_list = []
    for aa in list_aas:
        coords = cmd.get_coords('(' + aa_sel + ' and resi ' + str(aa) + ')')
        # if residue is unresolved, get_coords returns NoneType
        # therefore, skip unresolved residues (fill w/ None)
        if coords is None:
            centroid = (aa, None, None, None)
            centroid_list.append(centroid)
        else:
            xyz = np.mean(coords, axis=0)
            centroid = (aa, xyz[0], xyz[1], xyz[2])
            centroid_list.append(centroid)
    cmd.reinitialize()
    
    # turn list into pandas dataframe of xyz centroid coords
    df_centroids = pd.DataFrame(data=centroid_list,
                                columns=['aa_num', 'x', 'y', 'z'])
    # pymol sometimes makes numbers into strings -- make sure they are int
    df_centroids['aa_num'] = df_centroids['aa_num'].astype(int)
    num_resolved = df_centroids.loc[~df_centroids['x'].isna()].shape[0]

    if out_csv is None:
        out_csv = pdb_id + '_centroids.csv'
    # save centroid coordinates as csv
    if save_centroids:
        df_centroids.to_csv(out_csv, index=False)
    # print statistics
    print('# residues in list_aas: ' + str(len(list_aas)))
    print('# centroids in df_centroids: ' + str(df_centroids.shape[0]))
    print('# centroids resolved: ' + str(num_resolved))
    if return_df:
        print('get_centroids done for ' + pdb_id)
        return df_centroids
    else:
        return ('get_centroids done for ' + pdb_id)

#%% get_pairwise_dist() - calculates pairwise dist between resolved centroids

def get_pairwise_dist(df_centroids, aa_int=None):
    """
    This function calculates pairwise distances from centroid coordinates
    that were generated from get_centroids().
    Use pdb_id to call the correct centroids.csv output file
    Returns a pandas df of pairwise distances with columns/index as AA pos
    
    df_centroids: pandas dataframe containing 4 columns of
        ['aa_num', 'x', 'y', 'z'].
    aa_int: tuple of (aa_min, aa_max) defining the aas to calculate pdists
        the default is None, which takes the min/max of df_centroids
    """

    # check for correct columns in df_centroids
    list_cols = df_centroids.columns.tolist()
    if not all(col in list_cols for col in ['aa_num', 'x', 'y', 'z']):
        raise Exception('df_centroids is missing an essential column id')
    # make sure aa_num is the correct dtype (pymol uses strings)
    df_centroids['aa_num'] = df_centroids['aa_num'].astype('int64')

    # isolate the desired amino acid interval
    # default is all resolved residues in df_centroids
    if aa_int is None:
        # remove unresolved residues (xyz = NaN) before finding aa min/max
        df_aaint = df_centroids.loc[~df_centroids.isnull().any(axis=1)].copy()
        aa_min = df_aaint['aa_num'].min()
        aa_max = df_aaint['aa_num'].max()
    else:
        aa_min = aa_int[0]
        aa_max = aa_int[1]
        df_aaint = df_centroids.loc[df_centroids['aa_num'].between(aa_min, aa_max)].copy()
        df_aaint = df_aaint.loc[~df_aaint.isnull().any(axis=1)].copy()
        if df_aaint['aa_num'].min() != aa_min:
            print('Warning! User aa_min input was ' + str(aa_min))
            print('But first resolved AA was ' + str(df_aaint['aa_num'].min()))
        if df_aaint['aa_num'].max() != aa_max:
            print('Warning! User aa_max input was ' + str(aa_max))
            print('But last resolved AA was ' + str(df_aaint['aa_num'].max()))
    # calculate all pairwise distances in euclidean 3d space
    pairwise = sp_sd.pdist(df_aaint[['x','y','z']], 'euclidean')
    # turn condensed matrix into square-form matrix
    pairwise = sp_sd.squareform(pairwise)
    # convert to pandas df with index/col as aa numbers
    df_pwdist = pd.DataFrame(pairwise, index=df_aaint['aa_num'], columns=df_aaint['aa_num'])

    return df_pwdist

#%% gauss(), hill(), round() - gaussian/hill equations for CHUNKS scaling

# calculate distance component
def gauss(distance, std):
    """
    gaussian function to calculate the "distance component" for scaling
    input is distance (between residues i,j) and stdev (t)
    transforms the euclidean distance into a value between [0,1]
    shorter is closer to 1, longer closer to 0, relative to t
    t is arbitrary constant, represents a "soft distance threshold"
    LSD1 used 16 angstrom for t, CLUMPS used 6 angstrom
    see wikipedia gaussian function and the Gad Getz CLUMPS PNAS paper for ref
    """
    arg = -(distance * distance) / (2 * std * std)
    dist = np.exp(arg)
    return dist

# scale the summed LFC
def hill(lfc, m, theta):
    """
    hill function to scale the log2_fc numbers into range [0,1]
    sigmoidal func prevents highly enriched sgRNAs (e.g. jackpots) from having
    disproportionate influence vs. enriched but not jackpot sgRNAs
    input is the summed lfc (of sg1 and sg2), m, and theta
    theta controls the critical point (center), m controls steepness of function
    CLUMPS used m=3, theta=2, LSD1 used m=2, theta=3
    according to allison, m=3, theta=2 didn't work, so idk
    see the Gad Getz CLUMPS PNAS paper for reference
    """
    num = lfc**m
    denom = (lfc**m) + (theta**m)
    val = num/denom
    return val


#%% chunks_scale() - scale LFC with hill/gaussian functions (KCN, NZL modified)

def chunks_scale(df_lfc, cond, df_pwdist, std, m, theta):
    """
    Function to generate scaled LFC values using the "CLUMPS" method. The output
    is a pandas dataframe containing scaled LFC values corresponding
    to sgRNAs hitting resolved AAs in the pairwise distance matrix
    by KCN, modified by NZL

    Parameters
    ----------
    df_lfc : pandas dataframe
        The dataframe containing the sgRNA log2 fold change values, with
        each guide assigned to a residue in 'Clust_site' column
    cond : str
        The column name in the in_lfc file with the desired LFC values
    df_pwdist : pandas dataframe
        Pandas dataframe containing the pairwise distance values. this is
        generated from the get_pairwise_dist function
    std : int
        The standard deviation for the gaussian function.
    m : int
        Controls the steepness of the hill function.
    theta : int
        Controls the critical point (center) of the hill function.
    """
    
    # Set IDs
    ids = df_lfc['sgRNA_ID']

    # generate a pairwise matrix of summed LFC values with matrix operation
    df_sum = df_lfc[cond].values[:, None] + df_lfc[cond].values[None, :]
    df_sum = pd.DataFrame(index=ids, columns=ids, data=df_sum)
    # generate a matrix to remember the sign
    df_sign = df_sum.apply(lambda x: np.where(x > 0, 1, -1))
    # get abs value, scale w/ hill function, then adjust sign
    df_sum = df_sum.abs()
    df_sum = df_sum.apply(lambda x: hill(x, m, theta))
    df_hill = df_sum * df_sign
    # get AAs using sgRNA ids. not the best way but i'm overly anxious about the order
    aas = df_lfc.loc[df_lfc['sgRNA_ID'] == ids]['Clust_site']
    # retrieve pairwise distances from df_pwdist for the relevant sgRNAs
    df_gauss = df_pwdist.loc[aas, aas].copy()
    # reset index/column to match sgRNA ID rather than AA number
    df_gauss.set_index(keys=ids, drop=True, inplace=True)
    df_gauss.set_axis(labels=ids, axis=1, inplace=True)
    # scale the pairwise distances with the gaussian function
    df_gauss = df_gauss.apply(lambda x: gauss(x, std))
    # multiply gaussian component w/ hill component to get scaled LFCs
    df_scaled = df_gauss * df_hill
    # change data type to float
    df_scaled = df_scaled.astype(float)

    return df_scaled

#%% chunks_cluster() - generate clusters and annotate input csv file
# by NZL based on KCN code

def chunks_cluster(df_scaled, df_lfc, num_clusters, method='ward', metric='euclidean'):
    
    # Perform clustering
    link = sp_hca.linkage(df_scaled, method=method, metric=metric)
    clusters = sp_hca.fcluster(link, t=num_clusters, criterion='maxclust')
    df_clusters = pd.DataFrame(clusters, index=df_scaled.index, columns=['Group'])
    
    # Merge with df_lfc and return
    df_clusters = df_lfc.merge(df_clusters.reset_index(), on='sgRNA_ID')
    
    return df_clusters, link
    

