#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Functions for use in base_editor_library_annotation.py script.

@author: Ceejay Lee
"""

#%% Import Packages

from itertools import combinations
from Bio import SeqIO

#%% Functions

# Process .fasta into dictionaries with exons (+flanking introns), number as key
def processIFile(fastaFile):
    #import fasta file for comparison
    #exonDict is a dictionary with key = exon #, value = sequence
    #positionDict is a dictionary with key = exon #, value = (start,end) genomic position, inclusive
    exonDict = {}
    exonPosDict = {}

    for exon in SeqIO.parse(fastaFile,'fasta'): #SeqIO.parse takes a file handle (or filename) and format name, and returns a SeqRecord iterator
        exon_id = parseName(exon.id) #.id returns name of fasta chunk, parseName is a function that splits the name to only retain its numbering
        exon_seq = exon.seq #.seq returns the sequence of the fasta chunk
        exon_pos = parsePos(exon.description)
        exonDict[exon_id] = exon_seq #construct dictionary entry
        exonPosDict[exon_id] = exon_pos #construct dictionary for positions tuples
    
    return exonDict, exonPosDict



# HELPER FUNCTIONS

# Parse seqId so that only exon number remains
def parseName(seqId):
    seqId = int(((seqId.split(' ')[0]).split('_'))[4])
    return seqId

# Parse sequence description to pull out genomic coordinates
def parsePos(seqDesc):
    genPos = ((seqDesc.split(' ')[1]).split(':')[1]).split('-')
    return (int(genPos[0]), int(genPos[1]))

# Find all indexes of a particular base in edit window
def find_Indexes(base, seq_str, indexList):
    last_found = seq_str.rfind(base)
    if last_found == -1:
        return indexList
    indexList.append(last_found)
    new_seqstr = seq_str[:last_found]
    return find_Indexes(base, new_seqstr, indexList)

# Make list of all subsets
def generateSubsets(indexList):
    subsetList = []
    for i in range(1,len(indexList)+1):
        subsetList.extend(combinations(indexList,i))
    return subsetList

# Replace base at specified index to new base (note: original script has an error, fixed here)
def replaceBase(index_tuple, seq_str, new_base):
    for index in index_tuple:
        seq_str = seq_str[:index] + new_base + seq_str[index+1:]
    return seq_str

# Generate all possible combinations of mutations
def generateMutCom(seq_str, subsetList, new_base):
    return [replaceBase(index_tuple, seq_str, new_base) for index_tuple in subsetList]

# Function to translate codon to amino acid
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

