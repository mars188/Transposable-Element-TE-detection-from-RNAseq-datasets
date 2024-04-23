# Transposable Element (TE) detection from RNAseq datasets

This page describes the steps to perform the NYUAD Core Bioinformatics hands-on workshop on TE detection from RNA-seq dataset and visulizing the results. 

During this workshop, participants will learn about the concepts and steps involved in analyzing TE data from short read high throughput sequencing data.

Participants are expected to have some basic knowledge of TE and gene expression as well as some biological knowledge on the subject matter.

Although the material and the methods are designed to cater for the NYU Abu Dhabi High Performance Computing environment, it is possible to run this workshop on any other system provided that the neccessary software is installed, and the data is uploaded to that environment. It will also be neccessary to change the input and output directories/files to accomodate such an environment.

We will start off with fastq files sequencing files that are already quality trimmed (Illumina paired end short reads), and throughout the course of the workshop, we will learn how to process and analyze the data so that we end up with readcounts files that can be further used for making graphs and visulizing results. Moving forward, we will use the publically available mouse data downloaded from Benayoun et al. 2019 and can be found here. (https://genome.cshlp.org/content/29/4/697). 

To summarize, we will:

- Work on already quality trimmed reads 
- Neccessary files (Repeatmasker, annotation file, STAR index) needed for the analysis has already been downloaded and available in your folder.  
- Align the data to the reference mouse genome (version mm9).
- Count/Quantify RNAseq reads aligning to TEs
- Performs differential expression analysis on TEs

## Connecting to the HPC using a MAC/Linux machine and copying the data
- Open the "Terminal" app and type ssh NetID@jubail.abudhabi.nyu.edu
- Enter your NYU password when prompted and hit enter
- Once logged in, navigate to your personal "SCRATCH" directory or the any other subdirecotry where you want to run this workshop


## Connecting to the HPC using a Windows machine and copying the data.
- Open the "Putty" app, and fill out the fields as follows Host name jubail.abudhabi.nyu.edu, Port=22, and then click on "Open".
- Enter your NetId, and your password when prompted.
- Once logged in, navigate to your personal "SCRATCH" directory or the any other subdirecotry where you want to run this workshop.
  
## Setting up the environment and copying the data
We will be using the NYUAD High Performance Computing (HPC) cluster for this workshop, however, you can certainly run all of the analysis on any stand alone machine (server, personal laptop/Desktop etc.) provided that you have pre-installed the necessay software packages.

```
mkdir -p /scratch/$USER/TE_workshop
cd /scratch/$USER/TE_workshop
rsync -avP /scratch/ma5877/TE_workshop/data_files/ .
```
## The data
After performing the above action, you will find the following files in your directory
1. squire_fetch: This contains the Repeatmasker, annotation , chromosome fasta files downloaded from UCSC for mouse genome (mm9)
2. squire_clean: Filtered repeatmasker file for repeats of interest after collapsing overlapping repeats. We obtained these files by preforming **squire_clean** step from the workflow. We have already performed this step and you don't have to repeat.   
3. trimmed_reads: This directory contains the dataset that we will be using for this workshop, which are quality trimmed FASTQ sequencing files. They are publicly available to download from the Short Read Archive (SRA [https://www.ncbi.nlm.nih.gov/sra]) using the following accessions ERR674822, SRR24555542, SRR5604284. Our p1, p2, and p3 datasets are reduced from their original size, meaning that not all of the sequencing reads are available. The reason for this is due to time constraints. Running the complete dataset takes much longer than the time we have during this workshop.

As mentioned earlier, each one of you has been preallocated a particular dataset. Your task is to complete the analysis of this data and find out which one of the "case files" relates to the dataset that you have been given.
