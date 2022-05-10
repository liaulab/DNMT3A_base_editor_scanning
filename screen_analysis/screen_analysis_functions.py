#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Functions for screen analysis.

@author: Kevin Ngan
"""

#%% Import Packages

import os
import pathlib
import time
import csv
from collections import OrderedDict, Counter
import numpy as np
import pandas as pd
from Bio import SeqIO

#%% Batch Count

def count_reads(in_fastq, in_seqinfo='libseq.csv', out_counts='lib_counts.csv',
                out_np='np_counts.csv', out_stats='stats.txt',
                KEY_INTERVAL=(10,80), KEY='CGAAACACCG', KEY_REV='GTTTTAGA',
                DIR='FWD'):
    """
    Creates dictionary with guide counts from fastq file, writes to out_counts

    Parameters
    ----------
    in_fastq : str
        fastq file with reads containing sgRNA sequences
    in_seqinfo : str, optional
        csv file with guide seqs and annotations. The default is 'libseq.csv'.
        in_seqinfo must have headers, with 'sgRNA_seq' as the column header
        for the sgRNA sequences
    out_counts : str, optional
        csv output for guide counts. The default is 'lib_counts.csv'.
    out_np : str, optional
        csv output for non-perfect matches. The default is 'np_counts.csv'.
    out_stats : str, optional
        txt output for QC statistics. The default is 'stats.txt'.
    KEY_INTERVAL : tuple, optional
        Tuple of (start, end) for the key region. The default is (10,80).
    KEY : str, optional
        string indicates start of sgRNA. The default is 'CGAAACACCG'.
    KEY_REV : str, optional
        string indicates sgRNA if using rev orientation.
        The default is 'GTTTTAGA'.
    DIR : {'FWD', 'REV'}, optional
        direction to identify sgRNA sequence. The default is 'FWD'.
        end of hU6 (fwd) or start of sgRNA scaffold (rev)
    """

    # STEP 1: RESET VARIABLES, OPEN INPUT FILES FOR PROCESSING
    num_reads = 0 # total number of reads processed
    perfect_matches = 0 # count of reads with a perfect match to library
    np_matches = 0 # count of reads without a perfect match to library
    key_not_found = 0 # count of reads where key was not found
    np_list = [] # placeholder list for non-perfect matches
    # try opening input files, raise exception if not possible
    try:
        df_input = pd.read_csv(in_seqinfo)
        # sgRNA seq as key, count as value, start at 0 count for all sgRNAs
        dict_perfects = {sgRNA:0 for sgRNA in df_input['sgRNA_seq']}
    except:
        print('could not process' + in_seqinfo)
        return
    try:
        handle = open(in_fastq)
    except:
        print('could not find fastq file')
        return

    # STEP 2: PROCESS FASTQ FILE READS AND ADD COUNTS TO DICT
    readiter = SeqIO.parse(handle, 'fastq') # process reads in fastq file
    KEY_START, KEY_END = KEY_INTERVAL[0:2] # set the key interval
    # find sgRNA using FORWARD direction (default)
    if DIR == 'FWD':
        for record in readiter: # contains the seq and Qscore etc.
            num_reads += 1
            read_sequence = str.upper(str(record.seq))
            key_region = read_sequence[KEY_START:KEY_END]
            key_index = key_region.find(KEY)
            if key_index >= 0:
                start_index = key_index + KEY_START + len(KEY)
                guide = read_sequence[start_index:(start_index + 20)]
                if guide in dict_perfects:
                    dict_perfects[guide] += 1
                    perfect_matches += 1
                else:
                    np_matches += 1
                    np_list.append(guide)
            else:
                key_not_found += 1
    # find sgRNA using REVERSE direction
    elif DIR == 'REV':
        for record in readiter: # contains the seq and Qscore etc.
            num_reads += 1
            read_sequence = str.upper(str(record.seq))
            key_region = read_sequence[KEY_START:KEY_END]
            key_index = key_region.find(KEY_REV)
            if key_index >= 0:
                start_index = key_index + KEY_START
                guide = read_sequence[(start_index - 20):(start_index)]
                if guide in dict_perfects:
                    dict_perfects[guide] += 1
                    perfect_matches += 1
                else:
                    np_matches += 1
                    np_list.append(guide)
            else:
                key_not_found += 1
    else:
        raise Exception('ERROR! Direction not valid')
    handle.close()

    # STEP 3: SORT DICTIONARIES AND GENERATE OUTPUT FILES
    # sort perf matches (A-Z) with guides,counts as k,v and output to csv
    perfdict_sorted = OrderedDict(sorted(dict_perfects.items(),
                                         key=lambda t: t[0]))
    with open(out_counts, 'w', newline='') as perfectmatches_csv:
        mywriter = csv.writer(perfectmatches_csv, delimiter=',')
        for guide in perfdict_sorted:
            count = perfdict_sorted[guide]
            mywriter.writerow([guide,count])
    # now sort non-perfect matches and output to csv
    dict_np = Counter(np_list) # create dict for non-perfect guide matches
    nonperfdict_sorted = OrderedDict(sorted(dict_np.items(), key=lambda t: t[0]))
    with open(out_np, 'w', newline='') as nonperfectmatches_csv:
        mywriter = csv.writer(nonperfectmatches_csv, delimiter=',')
        for guide in nonperfdict_sorted:
            count = nonperfdict_sorted[guide]
            mywriter.writerow([guide,count])

    # STEP 4: CALCULATE STATS AND GENERATE STAT OUTPUT FILE

    # percentage of guides that matched perfectly
    pct_perfmatch = round(perfect_matches/float(perfect_matches + np_matches) * 100, 1)
    # percentage of undetected guides (no read counts)
    guides_with_reads = np.count_nonzero(list(dict_perfects.values()))
    guides_no_reads = len(dict_perfects) - guides_with_reads
    pct_no_reads = round(guides_no_reads/float(len(dict_perfects.values())) * 100, 1)
    # skew ratio of top 10% to bottom 10% of guide counts
    top_10 = np.percentile(list(dict_perfects.values()), 90)
    bottom_10 = np.percentile(list(dict_perfects.values()), 10)
    if top_10 != 0 and bottom_10 != 0:
        skew_ratio = top_10/bottom_10
    else:
        skew_ratio = 'Not enough perfect matches to determine skew ratio'

    # write analysis statistics to stat_file
    with open(out_stats, 'w') as statfile:
        statfile.write('Number of reads processed: ' + str(num_reads) + '\n')
        statfile.write('Number of reads where key was not found: ' + str(key_not_found) + '\n')
        statfile.write('Number of perfect guide matches: ' + str(perfect_matches) + '\n')
        statfile.write('Number of nonperfect guide matches: ' + str(np_matches) + '\n')
        statfile.write('Number of undetected guides: ' + str(guides_no_reads) + '\n')
        statfile.write('Percentage of guides that matched perfectly: ' + str(pct_perfmatch) + '\n')
        statfile.write('Percentage of undetected guides: ' + str(pct_no_reads) + '\n')
        statfile.write('Skew ratio of top 10% to bottom 10%: ' + str(skew_ratio))
        statfile.close()

    return (in_fastq + ' processed')


def batch_count(in_batch, in_ref, subdir_counts=None, subdir_np=None,
                subdir_stats=None, dir_fastq=None, KEY_INTERVAL=(10,80),
                KEY='CGAAACACCG', KEY_REV='GTTTTAGA', DIR='FWD'):
    """
    Function to batch process FASTQ files based on input csv files. Requires
    two csv files: in_batch and in_ref. in_batch should have samplenames, fastq
    filenames, sample descriptions, and conditions (e.g. drug/DMSO)

    Parameters
    ----------
    in_batch : csv
        csv file with the sample names, fastq filenames, description, condition
    in_ref : csv
        CSV file with guide seqs and annotations.
    subdir_counts : TYPE
        DESCRIPTION.
    subdir_np : TYPE
        DESCRIPTION.
    subdir_stats : TYPE
        DESCRIPTION.
    dir_fastq : TYPE, optional
        DESCRIPTION. The default is None.
    """

    # Define directory paths
    batch_st = time.perf_counter()
    path = os.getcwd()
    if subdir_counts is None:
        subdir_counts = path
    else:
        subdir_counts = os.path.join(path, subdir_counts)
        if not os.path.exists(subdir_counts):
            pathlib.Path(subdir_counts).mkdir()
    if subdir_np is None:
        subdir_np = path
    else:
        subdir_np = os.path.join(path, subdir_np)
        if not os.path.exists(subdir_np):
            pathlib.Path(subdir_np).mkdir()
    if subdir_stats is None:
        subdir_stats = path
    else:
        subdir_stats = os.path.join(path, subdir_stats)
        if not os.path.exists(subdir_stats):
            pathlib.Path(subdir_stats).mkdir()
    if dir_fastq is None:
        dir_fastq = path
    else:
        dir_fastq = os.path.join(path, dir_fastq)

    # import batch csv and process samples with count_reads()
    df_batch = pd.read_csv(in_batch)
    for row in df_batch.itertuples():
        t_start = time.perf_counter()
        fastq = os.path.join(dir_fastq, row.fastq_file)
        counts = os.path.join(subdir_counts, row.sample_id + '_counts.csv')
        npcounts = os.path.join(subdir_np, row.sample_id + '_npcounts.csv')
        stats = os.path.join(subdir_stats, row.sample_id + '_stats.txt')
        count_reads(in_fastq=fastq, in_seqinfo=in_ref, out_counts=counts,
                    out_np=npcounts, out_stats=stats, KEY_INTERVAL=KEY_INTERVAL,
                    KEY=KEY, KEY_REV=KEY_REV, DIR=DIR)
        t_end = time.perf_counter()
        print(row.sample_id + ' processed in %.2f sec' % (t_end - t_start))
    batch_end = time.perf_counter()

    return ('all samples processed in %.2f min' % ((batch_end - batch_st) / 60))

#%% Batch Process
    
def batch_process(in_batch, in_ref, counts_path=None, counts_files=None,
                  stats_path=None, stats_files=None, merge_stats=True,
                  pDNA_ID=None, day0_ID=None, out_folder=None, out_prefix=''):

    """
    This function performs all preprocessing on the raw reads
    from count_reads and outputs the aggregated files

    """

    # import ref files and define variables/paths
    path = os.getcwd()
    df_batch = pd.read_csv(in_batch)
    df_ref = pd.read_csv(in_ref)
    if counts_path is None:
        counts_path = path
    if counts_files is None:
        df_batch['counts_files'] = df_batch['sample_id'] + '_counts.csv'
    if stats_path is None:
        stats_path = path
    if stats_files is None:
        df_batch['stats_files'] = df_batch['sample_id'] + '_stats.txt'

    # import csv files and generate df for raw read counts
    df_rawreads = df_ref.copy()
    for row in df_batch.itertuples():
        file = os.path.join(counts_path, row.counts_files)
        df_temp = pd.read_csv(file, names=['sgRNA_seq', row.sample_id])
        # merge on sgRNA_seq to aggregate columns
        df_rawreads = pd.merge(df_rawreads, df_temp, on='sgRNA_seq')
    
    # use guide IDs as the index
    df_rawreads = df_rawreads.set_index(keys='sgRNA_ID', drop=False)
    df_ref = df_ref.set_index(keys='sgRNA_ID', drop=False)
    
    # remove any guides with zero counts in plasmid library or day 0 from
    # further analysis
    if pDNA_ID is not None:
        df_pDNA_zeros = df_rawreads.loc[df_rawreads[pDNA_ID]==0].copy()
        if df_pDNA_zeros.shape[0] > 0:
            print('Following guides have zero count in pDNA and were filtered out:')
            print(df_pDNA_zeros['sgRNA_ID'])
            df_ref = df_ref.loc[df_rawreads[pDNA_ID]!=0]
            df_rawreads = df_rawreads.loc[df_rawreads[pDNA_ID]!=0]
    if day0_ID is not None:
        df_d0_zeros = df_rawreads.loc[df_rawreads[day0_ID]==0].copy()
        if df_d0_zeros.shape[0] > 0:
            print('Following guides have zero count in day0 and were filtered out:')
            print(df_d0_zeros['sgRNA_ID'])
            df_ref = df_ref.loc[df_rawreads[day0_ID]!=0]
            df_rawreads = df_rawreads.loc[df_rawreads[day0_ID]!=0]
        
    # perform log2 normalization (brian/broad method)
    df_log2 = df_ref.copy()
    for row in df_batch.itertuples():
        total_reads = df_rawreads[row.sample_id].sum()
        df_log2[row.sample_id] = df_rawreads[row.sample_id].apply(
            lambda x: np.log2((x * 1000000 / total_reads) + 1))
    # perform t0 normalization
    df_t0 = df_ref.copy()
    t0 = df_batch.loc[df_batch['condition'] == 't0']['sample_id']
    if t0.shape[0] < 1:
        raise Exception('no t0 condition')
    elif t0.shape[0] > 1:
        raise Exception('sorry i am not built for this you broke me')
    t0 = t0[0]
    for row in df_batch.itertuples():
        df_t0[row.sample_id] = df_log2[row.sample_id].sub(df_log2[t0])
    df_t0.drop(columns=t0, inplace=True) # drop the t0 col
    # average replicates by condition
    list_conds = df_batch['condition'].unique().tolist()
    list_conds.remove('t0')
    df_conds = df_ref.copy()
    for cond in list_conds:
        reps = df_batch.loc[df_batch['condition'] == cond]['sample_id'].tolist()
        if len(reps) > 1:
            df_conds[cond] = df_t0[reps].mean(axis=1)
        elif len(reps) == 1:
            df_conds[cond] = df_t0[reps]
        else:
            raise Exception('Error! u broke my script with ur reps')

    if merge_stats:
        df_stats = pd.DataFrame(columns=['parameters'])
        for row in df_batch.itertuples():
            file = os.path.join(stats_path, row.stats_files)
            df_temp = pd.read_csv(file, sep=':', names=['parameters', row.sample_id])
            df_stats = pd.merge(df_stats, df_temp, on='parameters', how='outer')

    # export files
    if out_folder is None:
        destpath = path
    else:
        destpath = os.path.join(path, out_folder)
        if not os.path.exists(destpath):
            pathlib.Path(destpath).mkdir()
    df_rawreads.to_csv(os.path.join(destpath, out_prefix + 'reads.csv'), index=False)
    df_log2.to_csv(os.path.join(destpath, out_prefix + 'log2.csv'), index=False)
    df_t0.to_csv(os.path.join(destpath, out_prefix + 't0norm_reps.csv'), index=False)
    df_conds.to_csv(os.path.join(destpath, out_prefix + 't0norm_conds.csv'), index=False)
    df_stats.to_csv(os.path.join(destpath, out_prefix + 'stats.csv'), index=False)

    return df_ref

#%% Compare conditions
    
def compare_conditions(comparisons, in_t0norm, df_ref,
                       out_comps='comparisons.csv', out_folder=None):
    """
    comparisons: list of tuples (name, sample1, sample2) for pairwise comparisons
        the comparison is always sample1 - sample2 (e.g. treatment - control)
    in_t0norm: csv file with condition-averaged t0 normalized values
    out_comps: output csv file with results

    comparison tuples MUST MATCH the column headers
    """

    # import data
    path = os.getcwd()
    df_import = pd.read_csv(in_t0norm)
    df_import = df_import.set_index(keys='sgRNA_ID', drop=False)
    df_comps = df_ref.copy()
    # perform treatment vs. control comparison
    for name, treatment, control in comparisons:
        df_comps[name] = df_import[treatment].sub(df_import[control])
    # export file
    if out_folder is None:
        destpath = path
    else:
        destpath = os.path.join(path, out_folder)
        if not os.path.exists(destpath):
            pathlib.Path(destpath).mkdir()
    df_comps.to_csv(os.path.join(destpath, out_comps), index=False)

    return 'analyze_comparisons done'

