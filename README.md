# Kaka/kea comparison

This repository contains the analysis code used for:
Martini D, Dussex N, Robertson B, Gemmell NJ, Knapp M. (2021) Evolution of the world’s only alpine parrot - genomic adaptation or phenotypic plasticity, behaviour and ecology. _Molecular Ecology_ [https://doi.org/10.1111/MEC.15978](https://doi.org/10.1111/MEC.15978)

A brief description of the contents of the repository will be provided below.

### 1- Assembly and Annotation
This section describes the steps that we took to go from raw Illumina reads to a clean genome assembly draft and annotation. Briefly, this includes the commands used for preprocessing the sequencing files, assembling with SOAPdenovo2, super-scaffolding with the Satsuma-Chromosemble pipeline and annotating the genome with the MAKER pipeline.

### 2- Gene family evolution
This part of the analysis includes the assignment of protein sequences of different species to gene families using the hmmscan function of HMMER and the TreeFAM database, the processing of these assignments into appropriate input files for phylogenetic analysis and CAFE, and the steps to run the gene family evolution tests in CAFE.

### 3- Phylogenetic analysis
This section goes from the production of a concatenated alignment of orthologous genes (using PRANK and custom scripts) to a description of the commands we used to build a phylogenetic tree with this dataset in RAxML and date the tree with MCMCtree.

### 4- Positive selection tests
This includes the scripts used to run both the free-ratio branch test and the branch-site test for positive selection in PAML4.8, plus all pre-processing and post-processing steps.

### 5- Genome realignment and variants analysis
This section contains scripts for whole-genome realignment and variant calling using BWA and SAMtools/BAMtools pipelines, demographic inference with PSMC, analysis of variant sites using SNPeff and VCFtools/R. 
