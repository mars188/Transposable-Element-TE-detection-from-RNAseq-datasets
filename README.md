# Transposable-Element-TE-detection-from-RNAseq-datasets

The repo is intended to accompany the NYUAD Core Bioinformatics hands-on workshop on Variant Detection and Annotation.

During this 2-day workshop, participants will learn about the steps involved in calling variants (SNPs and Indels) from short read high throughput sequencing data.

Participants are expected to have some basic knowledge of mutation analysis as well as some biological knowledge on the subject matter.

Although the material and the methods are designed to cater for the NYU Abu Dhabi High Performance Computing environment, it is possible to run this workshop on any other system provided that the neccessary software is installed, and the data is uploaded to that environment. It will also be neccessary to change the input and output directories/files to accomodate such an environment.

The workshop will take the format of case studies with the aim of discovering the underlying mutations in 3 separate patients that are exhibiting some disease phenotypes. Our starting point will be raw sequencing files (Illumina paired end short reads), and throughout the course of the workshop, we will learn how to process and analyze the data so that we end up with annotated VCF (Variant Calling Format) files.

More specifically, we will:

Perform Quality Checking and Quality Trimming (QC/QT) on the raw data.
Align the data to the reference human genome (version HG38).
Carry out the necessary alignment post processing steps following established best practice guides.
Call variants (SNPs and Indels) and filter the variants that have been called.
Annotate the Variants and attempt to establish causative mutations.
By the end of this workshop, participants should be able to replicate these analysis steps on any DNA sequencing dataset (WGS/WES/Panels) originating from short read sequencing technologies.

Enjoy!
