#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to annotate DNMT3A base editing sgRNA library with key information for
analysis. This version of the script assumes +4 to +8 editing window and C>T
editing.

@author: Nicholas Lue and Ceejay Lee
"""

#%% Import Packages

import os
import math

import numpy as np
import pandas as pd

from sgRNA_functions import (
    processIFile, find_Indexes, generateSubsets, generateMutCom, get_aa
)

#%% Key inputs

# Dictionary with keys being possible values of 'Gene' column (excluding
# controls) and values being the corresponding fasta files
geneDict = {'DNMT3A':'200225_DNMT3A2_Input.fasta'}

# Input annotations. Requires 6 columns: sgRNA_ID, sgRNA_seq, Gene, sgRNA_pos,
# sgRNA_strand (sense or antisense), Gene_strand (plus or minus)
guideInput = 'Library-input-v3.csv'

# Define domains (this info will be used to assign guides to domains).
# First, make dictionaries for each gene, where keys are domains and values are
# tuples defining domain boundaries (inclusive, 1-indexed resi #s). Any exonic
# guides not within any of these boundaries will be labeled 'Undefined'.
# Note: numbering used here corresponds to human DNMT3A isoform 2.
domainDict1 = {'Nterm':(1,102),
               'PWWP':(102.5,161.5),
               'Inter1':(162,292),
               'ADD':(292.5,425.5),
               'Inter2':(426,444),
               'MTase':(444.5,723)}
#domainDict2 = {'domain1':(start,end),
#               'domain2':(start,end)}
# Second, define master domain dictionary with keys being possible genes
# (make sure these match the keys in geneDict!) and values being above dict(s).
domainDict_all = {'DNMT3A':domainDict1}


#%% Analysis set up

# Read each gene's fasta file and parse exon sequences and genomic positions.
# Store within dictionaries where keys are genes and values are corresponding
# parsed information (also dictionaries).
exonDicts = {}
exonPosDicts = {}
for gene, fastaFile in geneDict.items():
    exonDicts[gene], exonPosDicts[gene] = processIFile(fastaFile)

# Prepare reference coding sequences for each gene (stored in dictionary like above)
cdsDict = {}
for gene in geneDict.keys():
    full_CDS = ''
    for i in range(len(exonDicts[gene])):
        seq = str(exonDicts[gene][i])
        firstIntronNuc = len(seq) - 55 # because slicing does not include second index
        full_CDS += seq[55:firstIntronNuc]
    cdsDict[gene] = full_CDS


# Prepare reference lists of exon lengths  (stored in dictionary like above)
lenDict = {}
for gene in geneDict.keys():
    lenDict[gene] = [(len(exonDicts[gene][i]) - 110) for i in range(len(exonDicts[gene]))]

# Read in guides, genomic positions, other required annotations.
# Store in pandas data frame. Each annotation produced by this script will be
# added as a column to annotated_df.
input_df = pd.read_csv(guideInput)
annotated_df = input_df.copy()


#%% Annotate genomic position of editing window for each guide (INCLUSIVE)

def annot_window(guidePos, guideSense, geneStrand):
    
    # Return (-1,-1) for controls
    if guidePos < 0:
        return (-1,-1)
    
    # If gene on + strand, sense = fwd, antisense = rev
    elif geneStrand == 'plus':
        if guideSense == 'sense':
            return (guidePos+3,guidePos+7)
        elif guideSense == 'antisense':
            return (guidePos-7,guidePos-3)
    
    # If gene on - strand, sense = rev, antisense = fwd
    elif geneStrand == 'minus':
        if guideSense == 'sense':
            return (guidePos-7,guidePos-3)
        elif guideSense == 'antisense':
            return (guidePos+3,guidePos+7)

annotated_df = annotated_df.assign(Edit_window=annotated_df.apply(
    lambda x: annot_window(x.sgRNA_pos, x.sgRNA_strand, x.Gene_strand), axis=1))

#%% Annotate the targeted exon (based on first position of protospacer)

def annot_exon(guidePos, geneID, exonPosDicts):
    
    # Return -1 for controls
    if guidePos < 0:
        return -1
    else:
        # Iterate over each exon (+flanking introns) and check if guide is inside
        for exon, exonPos in exonPosDicts[geneID].items():
            if (guidePos >= exonPos[0]) and (guidePos <= exonPos[1]):
                return exon
            elif (exon == len(exonPosDicts[geneID])-1):
                print('Encountered a guide that does not fall within any exon.')
                return 'None'
    
annotated_df = annotated_df.assign(Targeted_exon=annotated_df.apply(
    lambda x: annot_exon(x.sgRNA_pos, x.Gene, exonPosDicts), axis=1))

#%% Annotate window overlap
# Label whether guides (besides controls) target exonic or intronic regions.
# Exon if editing window completely within exon. If completely outside exon,
# label as Intron, 5'-UTR, or 3'-UTR as appropriate. If window spans exon-
# intron junction, label as Exon/Intron, Exon/5'-UTR, or Exon/3'-UTR.

def annot_overlap(exon_num, window, geneStrand, geneID, exonPosDicts):
        
    # Annotate controls
    if exon_num < 0:
        return 'Control'
    
    # Make a list of editing window genomic positions
    window_list = np.arange(window[0], window[1]+1, 1)

    # Get genomic coordinates of exon (inclusive)
    # Note that if gene is on minus strand, exon_i is actually end of exon
    exon_i = exonPosDicts[geneID][exon_num][0] + 55
    exon_f = exonPosDicts[geneID][exon_num][1] - 55
    # Make a list of all genomic positions within exon
    exon_list = np.arange(exon_i, exon_f+1, 1)
    
    # Check if there is complete overlap between editing window and exon
    if all(pos in exon_list for pos in window_list):
        return 'Exon'
    # Check if there is partial overlap between editing window and exon
    elif any(pos in exon_list for pos in window_list):
        if exon_num == 0:
            if (geneStrand == 'plus') and (window[0] < exon_i):
                return "Exon/5'-UTR"
            elif (geneStrand == 'minus') and (window[1] > exon_f):
                return "Exon/5'-UTR"
            else:
                return 'Exon/Intron'
        elif exon_num == len(exonPosDicts[geneID])-1:
            if (geneStrand == 'plus') and (window[1] > exon_f):
                return "Exon/3'-UTR"
            elif (geneStrand == 'minus') and (window[0] < exon_i):
                return "Exon/3'-UTR"
            else:
                return 'Exon/Intron'
        else:
            return 'Exon/Intron'
    # If neither above condition is true, window is completely outside exon
    else:
        if exon_num == 0:
            if (geneStrand == 'plus') and (window[0] < exon_i):
                return "5'-UTR"
            elif (geneStrand == 'minus') and (window[0] > exon_f):
                return "5'-UTR"
            else:
                return 'Intron'
        elif exon_num == len(exonPosDicts[geneID])-1:
            if (geneStrand == 'plus') and (window[0] > exon_f):
                return "3'-UTR"
            elif (geneStrand == 'minus') and (window[0] < exon_i):
                return "3'-UTR"
            else:
                return 'Intron'
        else:
            return 'Intron'

annotated_df = annotated_df.assign(Win_overlap=annotated_df.apply(
    lambda x: annot_overlap(x.Targeted_exon, x.Edit_window,
                            x.Gene_strand, x.Gene, exonPosDicts), axis=1))

#%% Annotate with number of 'C's in editing window and whether C is present

annotated_df = annotated_df.assign(C_count=annotated_df.apply(
    lambda x: x.sgRNA_seq[3:8].count('C'), axis=1))

def annot_isC(C_count):
    
    #Label whether there is at least one C in editing window
    if C_count >= 1:
        return 'C'
    elif C_count == 0:
        return 'No_C'
    else: print('This should not occur')

annotated_df = annotated_df.assign(is_C=annotated_df.apply(
    lambda x: annot_isC(x.C_count), axis=1))

#%% Annotate splice site targeting guides
# Four categories: 'Splice-donor', 'Splice-acceptor', 'None', and 'Control'
# 'Splice-donor' means guide edits the 'GT' at 5' of intron
# 'Splice-acceptor' means guide edits the 'AG' at 3' of intron
# For C editor, only antisense guides can qualify. Strategy here is to check
# (1) that guide is antisense, (2) that its editing window overlaps with
# the 'C' at either splice site within that exon, (3) that there is a 'C' there.

def annot_splice(winOverlap, guideStrand, geneStrand, window, geneID, exon,
                 exonDicts, exonPosDicts):
    
    # Annotate controls
    if winOverlap == 'Control':
        return 'Control'
    else:
        if guideStrand == 'antisense':
            # Define exon boundaries (stored as tuple)
            exon_coords = exonPosDicts[geneID][exon]
            # Identify genomic location of splice site editable bases
            if geneStrand == 'plus':
                donor_loc = exon_coords[1]-54
                acceptor_loc = exon_coords[0]+54
            elif geneStrand == 'minus':
                donor_loc = exon_coords[0]+54
                acceptor_loc = exon_coords[1]-54
            
            # Check that there is an editable base there
            # Define exon+flanking intron sequence (given 5' to 3' along gene)
            exon_seq = str(exonDicts[geneID][exon]).upper()
            donor_present = ('G' == exon_seq[-55])
            acceptor_present = ('G' == exon_seq[54])
            
            # Expand editing window boundaries into a list of all positions
            window_list = np.arange(window[0], window[1]+1, 1)
            
            # Initialize string to hold annotation (doing it this way just in
            # case somehow a guide hits both an acceptor and donor for small exon)
            splice_note = ''
            
            # Check for donor editing
            if all([donor_loc in window_list, donor_present,
                    exon!=(len(exonPosDicts[geneID])-1)]):
                splice_note += 'Splice-donor'
            # Check for acceptor editing
            if all([acceptor_loc in window_list, acceptor_present, exon!=0]):
                if len(splice_note) == 0:
                    splice_note += 'Splice-acceptor'
                else:
                    splice_note += '/Splice-acceptor'
            
            # Return annotation
            if len(splice_note) == 0:
                return 'None'
            else:
                return splice_note
        else:
            return 'None'

annotated_df = annotated_df.assign(Splice_check=annotated_df.apply(
    lambda x: annot_splice(x.Win_overlap, x.sgRNA_strand, x.Gene_strand,
                           x.Edit_window, x.Gene, x.Targeted_exon,
                           exonDicts, exonPosDicts), axis=1))

#%% Annotate residue(s) and CDS position(s) targeted by each guide

# Process each list of exon lengths in lenDict into a list of CDS cumulative
# nucleotide numbers (stored in dictionary with key being gene as before)
shiftDict = {}
for gene, lenList in lenDict.items():
    shiftList = []
    for i in range(len(lenList)):
        shiftList.append(sum(lenList[:i]))
    shiftDict[gene] = shiftList
    
def annot_target_CDS(winOverlap, window, exon_num, geneStrand, geneID,
                     exonPosDicts, shiftDict):
    
    # Assign (-1,-1) to all guides whose windows do not overlap exons
    if not 'Exon' in winOverlap:
        return (-1,-1)
    
    #i and f refer to lower and upper boundaries in genomic coordinates
    #Note if gene is along minus strand, i is downstream of f along the gene
    window_i = window[0]
    window_f = window[1]
    exon_i = exonPosDicts[geneID][exon_num][0] + 55
    exon_f = exonPosDicts[geneID][exon_num][1] - 55
    if window_f > exon_f:
        window_f = exon_f
    if window_i < exon_i:
        window_i = exon_i
    
    #Find targeted interval along CDS (1-indexed, inclusive boundaries)
    if geneStrand == 'plus':
        CDS_i = shiftDict[geneID][exon_num] + (window_i - exon_i + 1)
        CDS_f = shiftDict[geneID][exon_num] + (window_f - exon_i + 1)
    elif geneStrand == 'minus':
        CDS_i = shiftDict[geneID][exon_num] + (exon_f - window_f + 1)
        CDS_f = shiftDict[geneID][exon_num] + (exon_f - window_i + 1)

    return (CDS_i, CDS_f)

annotated_df = annotated_df.assign(Target_CDS=annotated_df.apply(
    lambda x: annot_target_CDS(x.Win_overlap, x.Edit_window, x.Targeted_exon,
                               x.Gene_strand, x.Gene, exonPosDicts,
                               shiftDict), axis=1))

def annot_target_resi(targetCDS):
    
    # Only consider exonic guides
    if targetCDS == (-1,-1):
        return -1
    
    CDS_i = targetCDS[0]
    CDS_f = targetCDS[1]
    
    aa_i = math.ceil(CDS_i / 3)
    aa_f = math.ceil(CDS_f / 3)
    
    return (aa_i, aa_f)

annotated_df = annotated_df.assign(Target_resi=annotated_df.apply(
    lambda x: annot_target_resi(x.Target_CDS), axis=1))

#%% Annotate mutation type for exonic guides with C

def annot_mutations(guideID, winOverlap, num_C, window, guideSense, geneID,
                    spliceCheck, cdsDict):
    # Returns tuples (guideID, Mut_type, Mut_list_all, Mut_list)
    
    # Annotate controls
    if winOverlap == 'Control':
        return (guideID, 'Control', 'None', 'None')
    
    # Annotate guides without a C in the editing window
    if num_C == 0:
        if 'Exon' in winOverlap:
            return (guideID, 'No_C/Exon', 'None', 'None')
        else:
            return (guideID, 'No_C/Non-exon', 'None', 'None')
        
    # Annotate strictly non-exonic guides (not overlapping exons at all)
    if not 'Exon' in winOverlap:
        if 'Splice' in spliceCheck:
            return (guideID, 'Splice', 'None', 'None')
        else:
            return (guideID, 'Non-exon', 'None', 'None')
        
    # Window is the editing window boundaries, inclusive, 1-indexed
    window_seq = cdsDict[geneID][(window[0]-1):(window[1])]
    
    if guideSense == 'sense':
        base = 'C'
        new_base = 'T'
    elif guideSense == 'antisense':
        base = 'G'
        new_base = 'A'
    
    # All remaining guides (1) contain a C in the editing window and
    # (2) are either strictly Exonic or span an Exon/Non-exon junction
    
    # Thus, if there isn't an editable base within the exonic portion of the
    # window, the guide must span a junction, and the editable base must
    # be in the intronic part, so we annotate as Non-Exon (unless splice).
    if base not in window_seq:
        if 'Splice' in spliceCheck:
            return (guideID, 'Splice', 'None', 'None')
        else:
            return (guideID, 'Non-exon', 'None', 'None')
    
    # Generate list of all possible mutated versions of window_seq
    mutList = generateMutCom(window_seq,
                             generateSubsets(find_Indexes(base, window_seq, [])),
                             new_base)
    
    # Identify fully mutated window_seq (all C>T or G>A within window)
    maxMut = generateMutCom(window_seq,
                            [tuple(find_Indexes(base, window_seq, []))],
                            new_base)
    maxMut = maxMut[0] #Unpack from list
    
    #Assemble the codons surrounding the editing site
    window_0 = (window[0]-1, window[1]-1) #Convert editing window to zero-index
    first_codon_i = window_0[0] - window_0[0]%3
    last_codon_i = (window_0[1] + 3) - window_0[1]%3
    first_codon_adj = cdsDict[geneID][first_codon_i:window_0[0]]
    last_codon_adj = cdsDict[geneID][window_0[1]+1:last_codon_i]
    CDS_codonstr = first_codon_adj + window_seq + last_codon_adj #Done
    
    #Assemble analogous codons for potential mutated copies
    mut_codonstr = []
    for mut in mutList:
        mut_codonstr.append(first_codon_adj + mut + last_codon_adj)
    maxMut_str = first_codon_adj + maxMut + last_codon_adj #fully mutated
    
    #Break apart into individual codons
    CDS_codons = [CDS_codonstr[i:i+3] for i in range(0, len(CDS_codonstr), 3)]
    
    #Iterate through list of all possible combinations of base changes and
    #compile list of non-synonymous mutations
    all_mutlist = []
    for mut_seq in mut_codonstr:
        mut_codons = [mut_seq[i:i+3] for i in range(0, len(mut_seq), 3)]
        
        singlecom_mutlist = []
        for i, codon in enumerate(CDS_codons):
            aa_index = int((window_0[0]-window_0[0]%3)/3 + 1 + i) #1-index
            aa_CDS = get_aa(codon)
            aa_mut = get_aa(mut_codons[i])
            if aa_CDS == aa_mut:
                continue
            else:
                mut_note = str(aa_CDS+str(aa_index)+aa_mut)
                singlecom_mutlist.append(mut_note)
        
        if len(singlecom_mutlist) > 0:
            all_mutlist.extend(singlecom_mutlist)
    
    #Make analogous list for fully mutated
    max_all_mutlist = []
    max_mut_codons = [maxMut_str[i:i+3] for i in range(0, len(maxMut_str), 3)]
    for i, codon in enumerate(CDS_codons):
        aa_index = int((window_0[0]-window_0[0]%3)/3 + 1 + i) #1-index
        aa_CDS = get_aa(codon)
        aa_mut = get_aa(max_mut_codons[i])
        if aa_CDS == aa_mut:
            continue
        else:
            mut_note = str(aa_CDS+str(aa_index)+aa_mut)
            max_all_mutlist.append(mut_note)
    
    #Prep annotation for all possible mutations (combinations of base edits)
    if len(all_mutlist) > 0:
        final_all_mutlist = set(all_mutlist)
        all_annot = ''
        for index, mut in enumerate(final_all_mutlist):
            if index == 0:
                temp = mut
            else: temp = ', ' + mut
            all_annot += temp
    else:
        all_annot = 'None'
    
    # Prep fully edited annotation, including mutation type, and produce output
    if len(max_all_mutlist) > 0:
        final_list = ''
        nonsense_check = False
        for index, mut in enumerate(max_all_mutlist):
            if '*' in mut:
                nonsense_check = True
            if index == 0:
                temp = mut
            else: temp = ', ' + mut
            final_list += temp
        if nonsense_check:
            return (guideID, 'Nonsense', all_annot, final_list)
        else:
            if 'Splice' in spliceCheck:
                return (guideID, 'Splice', all_annot, final_list)
            else:
                return (guideID, 'Missense', all_annot, final_list)
    else:
        if 'Splice' in spliceCheck:
            return (guideID, 'Splice', all_annot, 'None')
        else:
            return (guideID, 'Silent', all_annot, 'None')

mut_notes = annotated_df.apply(
    lambda x: annot_mutations(x.sgRNA_ID, x.Win_overlap, x.C_count,
                              x.Target_CDS, x.sgRNA_strand, x.Gene,
                              x.Splice_check, cdsDict), axis=1)
mut_notes_df = pd.DataFrame(mut_notes.tolist(),
                            columns=['sgRNA_ID','Mut_type','Mut_list_all','Mut_list'])
annotated_df = annotated_df.merge(mut_notes_df, on='sgRNA_ID')

#%% Assign each guide a position (residue numbering)
# All guides apart from exonic ones: set to -1
# Exonic guides: assign position based on 1-indexed protein sequence
    # No C: Use the center of the editing window (convert to residue number)
    # Yes C:
        # Silent: Use the center of the editing window
        # Missense: Take mutated residue numbers and average
        # Nonsense: Take nonsense mutation residue numbers (avg if multiple)

def annot_editsite(window, mutType, mutList):
    
    # Set edit site to -1 if non-exonic or control
    if mutType in ['Non-exon', 'No_C/Non-exon', 'Control', 'Splice']:
        return -1
    elif (mutType == 'No_C/Exon') or (mutType == 'Silent'):
        avg_site = sum(window) / len(window)
        # If window contains odd number of positions, avg site is discrete base
        if (window[1]-window[0])%2 == 0:
            return math.ceil(avg_site / 3)
        else:
            # Otherwise avg site is between two bases. Find residues each maps to
            upper_site = math.ceil(avg_site)
            lower_site = math.floor(avg_site)
            upper_aa = math.ceil(upper_site/3)
            lower_aa = math.ceil(lower_site/3)
            # If both bases map to same residue, assign to that residue
            if upper_aa == lower_aa:
                return upper_aa
            # Otherwise take average
            else:
                return (upper_aa+lower_aa)/2
    elif (mutType == 'Nonsense'):
        muts_parsed = mutList.split(', ')
        resi_nums = [int(mut[1:-1]) for mut in muts_parsed if '*' in mut]
        return sum(resi_nums)/len(resi_nums)
    elif (mutType == 'Missense'):
        muts_parsed = mutList.split(', ')
        resi_nums = [int(mut[1:-1]) for mut in muts_parsed]
        resi_nums_unique = set(resi_nums)
        return sum(resi_nums_unique)/len(resi_nums_unique)
    else:
        print('This should not occur!')

annotated_df = annotated_df.assign(Edit_site=annotated_df.apply(
    lambda x: annot_editsite(x.Target_CDS, x.Mut_type, x.Mut_list), axis=1))


#%% Annotate domain

# Label all exonic guides by domain (e.g. any guides with an Edit_site not -1).
# For non-exonic guides, label as Control, Splice, Intron, 5'-UTR, or 3'-UTR.
# Exonic guides that do not fall within a provided domain are labeled 'Undefined'.
def annot_domain(site, winOverlap, mutType, geneID, domainDict_all):
    
    # Annotate controls
    if mutType == 'Control':
        return 'Control'
    # Annotate splice-site targeting guides
    elif mutType == 'Splice':
        return 'Splice'
    # Annotate non-exonic guides
    elif 'Non-exon' in mutType:
        # Assign Intron, 5'-UTR, or 3'-UTR as appropriate
        if "5'-UTR" in winOverlap:
            return "5'-UTR"
        elif "3'-UTR" in winOverlap:
            return "3'-UTR"
        elif 'Intron' in winOverlap:
            return 'Intron'
        else:
            print('This should not occur!')
            return
    # If none of above conditions is true, the guide must be exonic
    else:
        # Use gene annotation to pull out correct domain boundary dictionary
        bounds_Dict = domainDict_all[geneID]
        # Iterate over domains and check if guide lies within each one
        for domain, bounds in bounds_Dict.items():
            if bounds[0] <= site <= bounds[1]:
                return domain
            else:
                continue
        # If guide was not found in any domain, label as 'Undefined'
        return 'Undefined'

annotated_df = annotated_df.assign(Domain=annotated_df.apply(
    lambda x: annot_domain(x.Edit_site, x.Win_overlap, x.Mut_type, x.Gene,
                           domainDict_all), axis=1))


#%% Correct isoform numbering

# Previous annotations have all used short isoform numbering (out of 723).
# Create separate annotations with short vs. long isoform numbering

# First, copy Mut_list column to Mut_list_3A2
annotated_df['Mut_list_3A2'] = annotated_df['Mut_list']

# Function to convert Mut_list to long isoform numbering
def annot_3A1_mutlist(mutList):
    if mutList == 'None':
        return 'None'
    else:
        outlist = []
        for mut in mutList.split(', '):
            resi_3a2 = int(mut[1:-1]) #Residue number in 3A2 numbering
            #If residue number is < 25 then it is not in 3A1
            if resi_3a2 < 25:
                outlist.append('Not_in_3A1')
            else:
                outlist.append(mut[0]+str(resi_3a2+189)+mut[-1])
        #Now make output
        outstr = ', '.join(outlist)
        return outstr

# Annotate long isoform mutations
annotated_df = annotated_df.assign(Mut_list_3A1 = annotated_df.apply(
    lambda x: annot_3A1_mutlist(x.Mut_list), axis=1))

# Second, copy Edit_site to Edit_site_3A2
annotated_df['Edit_site_3A2'] = annotated_df['Edit_site']

# Function to convert Edit_site to long isoform numbering
def annot_corr_numbering(editSite):
    # Return -1 for any non-exon, control, splice guides or for guides that
    # target a residue that's not in the long isoform (short isoform-specific)
    if editSite < 25:
        return -1
    else:
        return editSite + 189

# Annotate long isoform numbering
annotated_df = annotated_df.assign(Edit_site_3A1 = annotated_df.apply(
    lambda x: annot_corr_numbering(x.Edit_site), axis=1))





#%% Export to .csv

# List of key columns to export
key_cols = ['sgRNA_ID', 'sgRNA_seq', 'Gene', 'sgRNA_pos', 'sgRNA_strand',
            'Edit_window', 'Targeted_exon', 'Win_overlap', 'C_count',
            'Splice_check', 'Mut_type', 'Domain',
            'Mut_list_3A1', 'Edit_site_3A1', 'Mut_list_3A2', 'Edit_site_3A2']

annotated_df.to_csv('Annotated_Library_Full.csv', index=False)

# This generates Supplementary Table 1.
annotated_df[key_cols].to_csv('Annotated_Library_Key_Columns.csv', index=False)

