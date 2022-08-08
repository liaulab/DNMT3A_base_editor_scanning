# Base editor scanning charts the DNMT3A activity landscape
Repository containing code used in Lue et al. Base editor scanning charts the DNMT3A activity landscape (2022). Preprint at _bioRxiv_, doi: [10.1101/2022.04.12.487946](https://doi.org/10.1101/2022.04.12.487946). Additional data or code is available from the corresponding author, Dr. Brian B. Liau (<liau@chemistry.harvard.edu>), upon request.

## Abstract
DNA methylation is critical for regulating gene expression, necessitating its accurate placement by enzymes such as the DNA methyltransferase DNMT3A. Dysregulation of this process is known to cause aberrant development and oncogenesis, yet how DNMT3A is regulated holistically by its three domains remains challenging to study. Here we integrate base editing with a DNA methylation reporter to perform in situ mutational scanning of DNMT3A in cells. We identify mutations throughout the protein that perturb function, including ones at an interdomain interface that block allosteric activation. Unexpectedly, we also find mutations in the PWWP domain, a histone reader, that modulate enzyme activity despite preserving histone recognition and protein stability. These effects arise from altered PWWP domain DNA affinity, which we show is a noncanonical function required for full activity in cells. Our findings highlight mechanisms of interdomain crosstalk and demonstrate a generalizable strategy to probe sequence-activity relationships of nonessential chromatin regulators.

## Repository Contents
This repository contains code used for next-generation sequencing data analysis and evolutionary conservation analysis. Specific file or directory names have been replaced where relevant with generic names. Additionally, key input files are provided in the appropriate subdirectories (or otherwise in the main publication), as detailed below. Raw fastq files are available upon request. Contents are organized into the following subdirectories:
1. Base editing sgRNA library annotation code
2. Base editor scanning data analysis code
3. Genotyping analysis code
4. Reporter bisulfite sequencing analysis code
5. PWWP evolutionary conservation analysis code
6. Genomics analysis code

## 1. Base editing sgRNA library annotation code
Code used to annotate sgRNAs with key information for analysis (i.e., to generate Supplementary Table 1). Two inputs are included in this directory:
- 200225_DNMT3A2_Input.fasta, which contains the sequences of DNMT3A isoform 2 (NM_153759.3, hg38) coding sequence exons. Each record corresponds to an exon with flanking intronic sequences (55 bp upstream and 55 bp downstream).
- Library-input-v3.csv, which contains the sequences of library sgRNAs and base information.

## 2. Base editor scanning data analysis code
Code used to process raw reads and analyze DNMT3A base editor scanning results, including PWES analysis. 'Batch_info.csv' input file is included in this directory. Supplementary Table 2 provides the normalized read counts at the replicate level.

## 3. Genotyping analysis code
Code used to analyze genotyping results and construct allele tables. Inputs to scripts are CRISPResso2 output, as well as the .txt CRISPResso2 batch files (included in the subdirectory). Scripts are as follows:
- genotyping_1.py, for analysis of sgW698 (Fig. 1i and Extended Data Fig. 1f) and sgE756 and sgG532 (Fig. 2f and Extended Data Fig. 2i) genotyping.
- genotyping_2.py, for analysis of sgRNA validation genotyping (Fig. 4c and Extended Data Fig. 5b,d) and construction of allele table (Supplementary Table 5).
- genotyping_3.py, for analysis of base editing around E342 codon (Extended Data Fig. 7b).
- genotyping_4.py, for analysis of clonal genotyping for select sgRNAs (Extended Data Fig. 3a,b and Supplementary Table 6).

## 4. Reporter bisulfite sequencing analysis code
Code used to analyze reporter bisulfite amplicon sequencing data (Extended Data Fig. 1c). Input to script is CRISPResso2 output.

## 5. PWWP evolutionary conservation analysis code
Code used to perform PWWP evolutionary conservation analysis (Fig. 4e and Supplementary Table 7).

## 6. Genomics analysis code
Code used to analyze ChIP-seq and RRBS data. Raw and processed ChIP-seq and RRBS data have been deposited to NCBI GEO (GSE199890).
- ChIP-seq_data_processing.sh, for aligning and processing ChIP-seq raw reads to generate input-normalized bigWig files.
- RRBS_data_processing.sh, for trimming, aligning, and processing RRBS raw reads to generate bedGraph files.
- RRBS_and_ChIP-seq_analysis.py, for downstream analysis of genomics data.
- DMR_analysis.sh, for analysis of DMRs.

## Citations
- Vinyard, M. E. et al. CRISPR-suppressor scanning reveals a nonenzymatic role of LSD1 in AML. _Nat. Chem. Biol._ 15, 529–539 (2019).
- Gosavi, P. M. et al. Profiling the Landscape of Drug Resistance Mutations in Neosubstrates to Molecular Glue Degraders. _ACS Cent. Sci._ (2021).
- Hanna, R. E. et al. Massively parallel assessment of human variants with base editor screens. _Cell_ 184, 1064-1080.e20 (2021).
- Clement, K. et al. CRISPResso2 provides accurate and rapid genome editing sequence analysis. _Nat. Biotechnol._ 37, 224–226 (2019).
- Kamburov, A. et al. Comprehensive assessment of cancer missense mutation clustering in protein structures. _Proc. Natl. Acad. Sci._ 112, E5486–E5495 (2015).
- Eddy, S. R. Accelerated Profile HMM Searches. _PLoS Comput. Biol._ 7, e1002195 (2011).
- Suzek, B. E. et al. UniRef clusters: a comprehensive and scalable alternative for improving sequence similarity searches. _Bioinformatics_ 31, 926–932 (2015).
- El-Gebali, S. et al. The Pfam protein families database in 2019. _Nucleic Acids Res._ 47, D427–D432 (2019).
- Sievers, F. et al. Fast, scalable generation of high-quality protein multiple sequence alignments using Clustal Omega. _Mol. Syst. Biol._ 7, 539–539 (2011).
- Li, H. & Durbin, R. Fast and accurate short read alignment with Burrows–Wheeler transform. _Bioinformatics_ 25, 1754–1760 (2009).
- Picard (<https://github.com/broadinstitute/picard>)
- Danecek, P. et al. Twelve years of SAMtools and BCFtools. _Gigascience_ 10, giab008 (2021).
- Ramírez, F. et al. deepTools2: a next generation web server for deep-sequencing data analysis. _Nucleic Acids Res._ 44, W160–W165 (2016).
- Trim Galore (<https://github.com/FelixKrueger/TrimGalore>)
- Krueger, F. & Andrews, S. R. Bismark: a flexible aligner and methylation caller for Bisulfite-Seq applications. _Bioinformatics_ 27, 1571–1572 (2011).
- MethylDackel (<https://github.com/dpryan79/MethylDackel>)
- Quinlan, A. R. & Hall, I. M. BEDTools: a flexible suite of utilities for comparing genomic features. _Bioinformatics_ 26, 841–842 (2010).
- Condon, D. E. et al. Defiant: (DMRs: easy, fast, identification and ANnoTation) identifies differentially Methylated regions from iron-deficient rat hippocampus. _BMC Bioinform._ 19, 31 (2018).




