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

## Connecting to the HPC using a MAC machine
- Open the "Terminal" app and type ssh NetID@jubail.abudhabi.nyu.edu
- Enter your NYU password when prompted and hit enter
- Once logged in, navigate to your personal "SCRATCH" directory or the any other subdirecotry where you want to run this workshop


## Connecting to the HPC using a Windows machine
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
3. trimmed_reads: This directory contains the dataset that we will be using for this workshop, which are quality trimmed FASTQ sequencing files. They are from the **Aging epigenomics of male mice** and publicly available to download from the Short Read Archive (SRA [https://www.ncbi.nlm.nih.gov/sra]) using the following accessions SRS2248611, SRS2248612, SRS2248613, SRS2248599, SRS2248600, SRS2248601. These data have been downsized from their original dataset. The reason for this is due to time constraints. Running the complete dataset takes much longer than the time we have during this workshop.
4. yml file: This is the workflow file that we will make use of to run the analysis.

**NOTE:** There are a total 6 FASTQ files. Three files belong to the young mice (3-month) and other three came from the old mice (29-month). We will use young mice as control group to detect the differential TE expression in old mice.  

## Step 0: Downloading reference data from UCSC genome broswer / creating your own

You will need to run the following commands in order to download and clean up the the referenece data files.
```
squire_fetch: Downloads input files from RefGene and generates STAR index Only needs to be done once initially to
acquire genomic input files or if a new build is desired.
squire_clean: Filters Repeatmasker file for Repeats of interest, collapses overlapping repeats, and returns as BED file.
Optional: Incorporation of non-reference TE sequence
```
**NOTE:** For this workshop and to save the time, the data have already been downloaded and included in your data file. YOU DO NOT NEED TO RUN THIS STEP. 

## Step 1: Aligns RNAseq reads
```
module load gencore
module load gencore_biosails
module load Miniconda3/
source activate /scratch/gencore/conda3/envs/squire/
squire -h

biox run -w te_squire.yml --select_rules squire_map -o map.sh
hpcrunner.pl submit_jobs -i map.sh
```

Complete script 
```
squire Map \
**put forwrad read flag** {$READ1} \
**put reverse read flag {$READ2} \
-o {$self->root_out_dir}/squire_map/{$sample} \
-f {$self->root_out_dir}/squire_fetch \
**genome build** \
-r 65 \
-p 24 \
-v
```
## Step 2: Quantifies RNAseq reads aligning to TEs
```
biox run -w te_squire.yml --select_rules squire_count -o count.sh
hpcrunner.pl submit_jobs -i count.sh
```
Complete script
```
squire Count \
**flag for map folder** {$squire_map} \
**flag for clean folder** {$squire_clean \
-o {$squire_count \
**flag for fetch folder** {$squire_fetch \
-r 65 \
**genome build** \
-p 24 \
-v
```
## Step 3: Performs differential expression analysis on TEs
```
biox run -w te_squire.yml --select_rules squire_call -o call.sh
hpcrunner.pl submit_jobs -i call.sh
```
Complete script
```
squire Call \
-1 liver_29m* \
-2 liver_3m* \
-A 29m \
-B 3m \
-i {$squire_count \
-o {$squire_call/locus \
-p 24 \
-v && \
squire Call \
-1 liver_29m* \
-2 liver_3m* \
-A 29m -B 3m \
-i {$squire_count \
-o {$squire_call/subfamily \
-s \
-p 24 \
-v
```

# DEseq2 and Visualization of TE counted by SQuIRE Locus specific

## loading data 
##############################################################################
```
SQuIRE_TOTcounts_locus <- read.delim("./Input_Data/SQuIRE_gene_locusTE_counttable.txt")
View(SQuIRE_TOTcounts_locus)
```
## preparation of dataframe 
##########################################################################################
```
library("dplyr")
library("tibble")
SQuIRE_TOTcounts_locus <- SQuIRE_TOTcounts_locus %>%
  remove_rownames() %>% # this is needed if dataframe has been filtered already so the picked row numbers are considered rownames
  column_to_rownames(var = "gene_id") %>% # move specific columns as rownames
  select(liver_3m1_1.fastq, liver_3m2_1.fastq, liver_3m3_1.fastq,
         liver_29m1_1.fastq, liver_29m2_1.fastq, liver_29m3_1.fastq) %>% # reorder columns as called
  rename(liver_3m_1 = liver_3m1_1.fastq, liver_3m_2 = liver_3m2_1.fastq, liver_3m_3 = liver_3m3_1.fastq,
         liver_29m_1 = liver_29m1_1.fastq, liver_29m_2 = liver_29m2_1.fastq, liver_29m_3 = liver_29m3_1.fastq) # rename columns as called
View(SQuIRE_TOTcounts_locus)
cts_SQuIRE_TOTcounts_locus <- as.matrix(SQuIRE_TOTcounts_locus)
```

## preparation of annotation 
##########################################################################################
```
condition <- factor(c(rep("liver_3m", 3), rep("liver_29m", 3)))
coldata<- data.frame(row.names=colnames(cts_SQuIRE_TOTcounts_locus), condition)
View(coldata)

### should return TRUE

all(rownames(coldata) == colnames(cts_SQuIRE_TOTcounts_locus))
```

## DEseq2
##########################################################################################
```
library("DESeq2")
dds_SQuIRE_TOTcounts_locus <- DESeqDataSetFromMatrix(countData = cts_SQuIRE_TOTcounts_locus,
                                                      colData = coldata,
                                                      design = ~ condition)
dds_SQuIRE_TOTcounts_locus

dds_SQuIRE_TOTcounts_locus <- DESeq(dds_SQuIRE_TOTcounts_locus)
View(assay(dds_SQuIRE_TOTcounts_locus))
```
## extract normalized counts
```
library("dplyr")
library("tibble")
SQuIRE_normCounts_locus <- counts(dds_SQuIRE_TOTcounts_locus, normalized=TRUE)
df_SQuIRE_normCounts_locus <- as.data.frame(SQuIRE_normCounts_locus)
df_SQuIRE_normCounts_locus <- df_SQuIRE_normCounts_locus %>%
  rownames_to_column(var = "gene_id") # move rownames as specific column 
View(df_SQuIRE_normCounts_locus)

df_SQuIRE_normCounts_locus_TE_only <- df_SQuIRE_normCounts_locus[26316:31963,]
df_SQuIRE_normCounts_locus_TE_only <- df_SQuIRE_normCounts_locus_TE_only %>%
  rename(TE_ID = gene_id) # rename column using new_name = old_name syntax  
View(df_SQuIRE_normCounts_locus_TE_only)
```
TE_ID contains concatenated TE_chr | TE_start | TE_stop| TE_name | milliDiv | TE_ strand

# Variance stabilizing transformation
##########################################################################################
## Get vsd for PCA
```
vsd_SQuIRE_TOTcounts_locus <- vst(dds_SQuIRE_TOTcounts_locus, blind=TRUE)
vsd_SQuIRE_TOTcounts_locus_Gene_only <- vsd_SQuIRE_TOTcounts_locus[1:26315,]
vsd_SQuIRE_TOTcounts_locus_TE_only <- vsd_SQuIRE_TOTcounts_locus[26316:31963,]
```
## Principal Component Analysis on vst
```
plotPCA(vsd_SQuIRE_TOTcounts_locus, intgroup=c("condition"))
plotPCA(vsd_SQuIRE_TOTcounts_locus_Gene_only, intgroup=c("condition"))
plotPCA(vsd_SQuIRE_TOTcounts_locus_TE_only, intgroup=c("condition"))
```
## export pdf plot on TE
```
dir.create("./Results")
dir.create("./Results/DEseq2")
pdf("./Results/DEseq2/vsd_SQuIRE_TOTcounts_locus_TE_only.pdf", onefile = FALSE, paper = "special", width = 10, height = 7.5)
plotPCA(vsd_SQuIRE_TOTcounts_locus_TE_only, intgroup=c("condition"))
dev.off()
```
# compute Euclidian distance matrix 
##########################################################################################
```
sampleDists_SQuIRE_TOTcounts_locus <- dist(t(assay(vsd_SQuIRE_TOTcounts_locus)))
sampleDists_SQuIRE_TOTcounts_locus_Gene_only <- dist(t(assay(vsd_SQuIRE_TOTcounts_locus[1:26315,])))
sampleDists_SQuIRE_TOTcounts_locus_TE_only <- dist(t(assay(vsd_SQuIRE_TOTcounts_locus[26316:31963,])))

library("RColorBrewer")
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
library("pheatmap")
sampleDistMatrix_SQuIRE_TOTcounts_locus <- as.matrix(sampleDists_SQuIRE_TOTcounts_locus)
pheatmap(sampleDistMatrix_SQuIRE_TOTcounts_locus,
         clustering_distance_rows=sampleDists_SQuIRE_TOTcounts_locus,
         clustering_distance_cols=sampleDists_SQuIRE_TOTcounts_locus,
         col=colors)
sampleDistMatrix_SQuIRE_TOTcounts_locus_Gene_only <- as.matrix(sampleDists_SQuIRE_TOTcounts_locus_Gene_only)
pheatmap(sampleDistMatrix_SQuIRE_TOTcounts_locus_Gene_only,
         clustering_distance_rows=sampleDists_SQuIRE_TOTcounts_locus_Gene_only,
         clustering_distance_cols=sampleDists_SQuIRE_TOTcounts_locus_Gene_only,
         col=colors)
sampleDistMatrix_SQuIRE_TOTcounts_locus_TE_only <- as.matrix(sampleDists_SQuIRE_TOTcounts_locus_TE_only)
pheatmap(sampleDistMatrix_SQuIRE_TOTcounts_locus_TE_only,
         clustering_distance_rows=sampleDists_SQuIRE_TOTcounts_locus_TE_only,
         clustering_distance_cols=sampleDists_SQuIRE_TOTcounts_locus_TE_only,
         col=colors)
```
## export pdf plot on TE
```
pdf("./Results/DEseq2/sampleDistMatrix_SQuIRE_TOTcounts_locus_TE_only.pdf", onefile = FALSE, paper = "special", width = 10, height = 7.5)
sampleDistMatrix_SQuIRE_TOTcounts_locus_TE_only <- as.matrix(sampleDists_SQuIRE_TOTcounts_locus_TE_only)
pheatmap(sampleDistMatrix_SQuIRE_TOTcounts_locus_TE_only,
         clustering_distance_rows=sampleDists_SQuIRE_TOTcounts_locus_TE_only,
         clustering_distance_cols=sampleDists_SQuIRE_TOTcounts_locus_TE_only,
         col=colors)
dev.off()
```
## calculate DEG 
### the second group name represent the ctrl vs. which calculate difference
##########################################################################################
```
res_SQuIRE_TOTcounts_locus <- results(dds_SQuIRE_TOTcounts_locus, contrast = c("condition","liver_29m", "liver_3m"))
res_SQuIRE_TOTcounts_locus
summary(res_SQuIRE_TOTcounts_locus)

res_SQuIRE_TOTcounts_locus_Gene_only <- res_SQuIRE_TOTcounts_locus[1:26315,]
res_SQuIRE_TOTcounts_locus_TE_only <- res_SQuIRE_TOTcounts_locus[26316:31963,]

plotMA(res_SQuIRE_TOTcounts_locus, ylim=c(-7,7), alpha = 0.05)
plotMA(res_SQuIRE_TOTcounts_locus_Gene_only, ylim=c(-7,7), alpha = 0.05)
plotMA(res_SQuIRE_TOTcounts_locus_TE_only, ylim=c(-7,7), alpha = 0.05)

pdf("./Results/DEseq2/res_SQuIRE_TOTcounts_locus_TE_only_MAplot.pdf", onefile = FALSE, paper = "special", width = 10, height = 7.5)
plotMA(res_SQuIRE_TOTcounts_locus_TE_only, ylim=c(-7,7), alpha = 0.05)
dev.off()
```
# export DEseq2 results
##########################################################################################
```
library("dplyr")
library("tibble")
SQuIRE_DESeq2_locus <- as.data.frame(res_SQuIRE_TOTcounts_locus)
SQuIRE_DESeq2_locus <- SQuIRE_DESeq2_locus %>%
  rownames_to_column(var = "gene_id") # move rownames as specific column 
View(SQuIRE_DESeq2_locus)

SQuIRE_DESeq2_locus_TE_only <- SQuIRE_DESeq2_locus[26316:31963,]
SQuIRE_DESeq2_locus_TE_only <- SQuIRE_DESeq2_locus_TE_only %>%
  rename(TE_ID = gene_id) # rename column using new_name = old_name syntax  
View(SQuIRE_DESeq2_locus_TE_only)
```
TE_ID contains concatenated TE_chr | TE_start | TE_stop| TE_name | milliDiv | TE_ strand
##########################################################################################
##########################################################################################

# preparation of merged dataframe 
##########################################################################################
## merging with sort=FALSE keep the order of X list
```
DESeq2_normCounts_locus_TE_only <- merge(SQuIRE_DESeq2_locus_TE_only, df_SQuIRE_normCounts_locus_TE_only, by = "TE_ID", sort = FALSE)
View(DESeq2_normCounts_locus_TE_only)
```
## Separate a column in multiple column base on special characters
```
library("dplyr")
library("tidyr")
DESeq2_normCounts_locus_TE_only <- DESeq2_normCounts_locus_TE_only %>% 
  separate(TE_ID, c("chr", "start", "end", "repName", "repFamily", "repClass", "milliDiv", "strand", "extra"), sep = "([|:,])", extra = "merge", fill = "right", remove = FALSE)
View(DESeq2_normCounts_locus_TE_only)
```
## count by Class of TE 
```
library("dplyr")
DESeq2_normCounts_locus_TE_only_count <- DESeq2_normCounts_locus_TE_only %>%
  group_by(repClass) %>% count()
View(DESeq2_normCounts_locus_TE_only_count)
```
## subset by pvalue and count by Class of TE 
```
DESeq2_normCounts_locus_TE_only_pvalue005 <- subset(DESeq2_normCounts_locus_TE_only, pvalue< 0.05)
View(DESeq2_normCounts_locus_TE_only_pvalue005)
DESeq2_normCounts_locus_TE_only_pvalue005_count <- DESeq2_normCounts_locus_TE_only_pvalue005 %>%
  group_by(repClass) %>% count()
View(DESeq2_normCounts_locus_TE_only_pvalue005_count)
```
## subset by pvalue and by LTR repClass
```
DESeq2_normCounts_locus_TE_only_pvalue005_LTR <- subset(DESeq2_normCounts_locus_TE_only, pvalue< 0.05 & repClass=="LTR")
View(DESeq2_normCounts_locus_TE_only_pvalue005_LTR)
```
## heatmap on normalized counts
##########################################################################################
```
library("pheatmap")
dir.create("./Results/pheatmap")
pdf("./Results/pheatmap/DESeq2_normCounts_locus_TE_only_pvalue005.pdf", onefile = FALSE, paper = "special", width = 7, height = 6)
pheatmap(DESeq2_normCounts_locus_TE_only_pvalue005[,17:22], cluster_rows=TRUE, show_rownames=FALSE,
         cluster_cols=FALSE, scale = "row", labels_row = DESeq2_normCounts_locus_TE_only_pvalue005$repName)
dev.off()

pdf("./Results/pheatmap/DESeq2_normCounts_locus_TE_only_pvalue005_LTR.pdf", onefile = FALSE, paper = "special", width = 7, height = 6)
pheatmap(DESeq2_normCounts_locus_TE_only_pvalue005_LTR[,17:22], cluster_rows=TRUE, show_rownames=FALSE,
         cluster_cols=FALSE, scale = "row", labels_row = DESeq2_normCounts_locus_TE_only_pvalue005_LTR$repName)
dev.off()
```
## Enanched Volcano 
##########################################################################################
### prep color function based on repClass
```
keyvals.colour <- ifelse(
  DESeq2_normCounts_locus_TE_only$pvalue > 0.05 , 'grey80', # is a lighter grey
  ifelse(DESeq2_normCounts_locus_TE_only$repClass == "DNA", 'orange',
         ifelse(DESeq2_normCounts_locus_TE_only$repClass == "LTR", 'purple1',
                ifelse(DESeq2_normCounts_locus_TE_only$repClass == "LINE", 'olivedrab1',
                       ifelse(DESeq2_normCounts_locus_TE_only$repClass == "SINE", 'dodgerblue1',
                              'grey20')))))# anything else is the significant non-TEs
keyvals.colour[is.na(keyvals.colour)] <- 'grey80'
names(keyvals.colour)[keyvals.colour == 'grey80'] <- 'NotSig'
names(keyvals.colour)[keyvals.colour == 'grey20'] <- 'not-TEs'
names(keyvals.colour)[keyvals.colour == 'orange'] <- 'DNA'
names(keyvals.colour)[keyvals.colour == 'purple1'] <- 'LTR'
names(keyvals.colour)[keyvals.colour == 'olivedrab1'] <- 'LINE'
names(keyvals.colour)[keyvals.colour == 'dodgerblue1'] <- 'SINE'
View(keyvals.colour)

library("EnhancedVolcano")
dir.create("./Results/EnhancedVolcano")
pdf("./Results/EnhancedVolcano/DESeq2_normCounts_locus_TE_only_groupCol.pdf", onefile = FALSE, paper = "special", width = 7, height = 6)
EnhancedVolcano(DESeq2_normCounts_locus_TE_only,
                lab = paste(DESeq2_normCounts_locus_TE_only$repName, DESeq2_normCounts_locus_TE_only$repFamily, sep = ":"),
                x = 'log2FoldChange',
                y = 'pvalue',
                xlim = c(-9, 9),
                ylim = c(0, 19),
                title = 'TE_DE_locus',
                subtitle = '29m VS. 3m livers',
                pCutoff = 0.05,
                FCcutoff = 0,
                pointSize = 2.0,
                labSize = 2.0,
                colCustom = keyvals.colour,
                colAlpha = 0.6,
                legendPosition = 'bottom',
                legendLabSize = 14,
                legendIconSize = 4.0,
                labCol = 'black',
                labFace = 'bold',
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                max.overlaps = 30,
                colConnectors = 'black')
dev.off()
```

## Annotation by genomation package
########################################################################################################################################
```
library("genomation")

bed.file = "/Users/fm1442/Sadler Edepli Lab Dropbox/Filippo Macchi/BioInformatics/07_Utilities/Annotation_BED-files/refGene_Mmusulus_mm10.gz"
gene.parts = readTranscriptFeatures(bed.file)

GR_DESeq2_normCounts_locus_TE_only <- as(DESeq2_normCounts_locus_TE_only, "GRanges")
GR_DESeq2_normCounts_locus_TE_only_pvalue005 <- as(DESeq2_normCounts_locus_TE_only_pvalue005, "GRanges")
GR_DESeq2_normCounts_locus_TE_only_pvalue005_LTR <- as(DESeq2_normCounts_locus_TE_only_pvalue005_LTR, "GRanges")

annot_GR_DESeq2_normCounts_locus_TE_only = annotateWithGeneParts(GR_DESeq2_normCounts_locus_TE_only, gene.parts, strand=TRUE, intersect.chr=TRUE)
dir.create("./Results/Genomation")
pdf("./Results/Genomation/annot_GR_DESeq2_normCounts_locus_TE_only.pdf", onefile = FALSE, paper = "special", width = 10, height = 7.5)
plotTargetAnnotation(annot_GR_DESeq2_normCounts_locus_TE_only)
dev.off()
```
## take the numbers (with promoter > exon > intron precedence)
```
annot_GR_DESeq2_normCounts_locus_TE_only
annot_GR_DESeq2_normCounts_locus_TE_only_pvalue005 = annotateWithGeneParts(GR_DESeq2_normCounts_locus_TE_only_pvalue005, gene.parts, strand=TRUE, intersect.chr=TRUE)
pdf("./Results/Genomation/annot_GR_DESeq2_normCounts_locus_TE_only_pvalue005.pdf", onefile = FALSE, paper = "special", width = 10, height = 7.5)
plotTargetAnnotation(annot_GR_DESeq2_normCounts_locus_TE_only_pvalue005)
dev.off()
```
## take the numbers (with promoter > exon > intron precedence)
```
annot_GR_DESeq2_normCounts_locus_TE_only_pvalue005
annot_GR_DESeq2_normCounts_locus_TE_only_pvalue005_LTR = annotateWithGeneParts(GR_DESeq2_normCounts_locus_TE_only_pvalue005_LTR, gene.parts, strand=TRUE, intersect.chr=TRUE)
pdf("./Results/Genomation/annot_GR_DESeq2_normCounts_locus_TE_only_pvalue005_LTR.pdf", onefile = FALSE, paper = "special", width = 10, height = 7.5)
plotTargetAnnotation(annot_GR_DESeq2_normCounts_locus_TE_only_pvalue005_LTR)
dev.off()
```
## take the numbers (with promoter > exon > intron precedence)
annot_GR_DESeq2_normCounts_locus_TE_only_pvalue005_LTR

##########################################################################################
##########################################################################################
