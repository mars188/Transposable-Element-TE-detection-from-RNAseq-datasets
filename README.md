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
