# The chemosynthetic microbiome of deep-sea hydrothermal vents across space and time

## 16S rRNA sequences analysis

################################################################################
## Pipeline for the analysis of sequences using Dada2 modified from https://github.com/giovannellilab/dada2_pipeline. 
## Designed for medium size dadasets, can also be used with large dataset where multiple runs are merged together. 
##The specific steps used for merging of multiple runs are highlighted below. 
## For information contact Matteo Selci at the Rutgers University, selcimatteo@gmail.com 
## 10/08/2024
################################################################################

# The analysis starts from within a working forlder where one or more folders
# containing the raw fastq files are stored. The different raw reads folder in
# this pipeline are called raw_1 raw_2 etc.., and cotain the fastq.gz files.

## Loading the library used for the pipeline
library(dada2); packageVersion("dada2") # Amplicon Sequencing Variants inference
#library(DECIPHER) # Multiple sequence alignment and phylogenetic analysis
#library(ape) # importing and handling phylogenetic trees

## STATING THE DADA2 PIPELINE
## This pipeline is based on the pipeline published in:
## Callahan et al. 2016. Bioconductor Workflow for Microbiome Data Analysis: from raw reads to community analyses. F1000Res 5, 1492.
## Several steps have been modified and adapted to the specific analyses carried
## out in the Giovannelli Lab. Please be sure that waht we do makes sense and
## applies to your specific project. This pipeline is shared as is.


#########################################################################################
#### ANALYSIS BLOCK 1 ###################################################################
########################################################################################

path1 <- "/sequences" # CHANGE it to the directory containing the fastq files after unzipping
list.files(path1) # Verify the file list

# separare forward e reverse reads. Change the pattern if necessary
fnFs1 <- sort(list.files(path1, pattern=".fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
# The Name is selected as ending at the FIRST underscore. This can be changed
# according to the sample name selected
sample.names1 <- sapply(strsplit(basename(fnFs1), ".fastq"), `[`, 1)

# Place filtered files in path1/filtered/ subdirectory
filtFs1 <- file.path(path1, "filtered", paste0(sample.names1, "_F_filt.fastq.gz"))

## TRIMMING
# Checking the sequence quality. Different samples can be inspected by changing the
# numbers within the brackets. Refer to the original publication by Callahan et al. 2016
# for interpretation

png("/home/ms3468/biofilm_agrosche/dada2/run_13dec23/plot/QplotF_1.png", width=1600, height=1200)
plotQualityProfile(fnFs1[4:12]) # reverse sequences
dev.off()


# Select where the sequences should be truncated. usually where the median quality
# (the solid organge line in the plots) starts to drop

out1 <- filterAndTrim(fnFs1, filtFs1, truncLen=172, maxN=0, maxEE=2, truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=18)
head(out1)
mean(out1[,2]/out1[,1]) #between 0.9-0.95 it is fine. Increase the maxEE if you get less than 0.9 

# Verify how many reads are lost. It is important to try different parameters
# (truncLen and maxEE, truncQ) to preserve a large number of sequences, without compormising
# on the quality of the reads. In case of low quality very few sequences will
# pair in the later steps and been retained. This step is the most important
# in the entire Dada2 pipeline. feel free to try several parameter settings and
# choose the best settings.Refer to Callahan et al. 2016 for details.
head(out1)

# Once the filter and trimming has been completed to satisfaction, proceede to
# the error lerning. This is computationally intensive

errF1 <- learnErrors(filtFs1, multithread=20, randomize=TRUE)

png("/home/ms3468/biofilm_agrosche/dada2/run_13dec23/plot/errorFw_plot.png", width=1600, height=1200)
plotErrors(errF1, nominalQ=TRUE)
dev.off()


# Dereplication step
derepFs1 <- derepFastq(filtFs1, verbose=TRUE)


# Name the derep-class objects by the sample names
names(derepFs1) <- sample.names1


# Infer the ASVs using the error prifiles ;earnined in the previous steps
dadaFs1 <- dada(derepFs1, err=errF1, pool="pseudo", multithread=18)


# Create the ASV table from run_1
seqtab1 <- makeSequenceTable(dadaFs1)
dim(seqtab1) # Gives you info on the number of Samples and ASVs identified in run_1

# Inspect distribution of sequence lengths
#table(nchar(getSequences(seqtab1)))

# Drop the tail sequences that are either too short or too long compared to the majority
# Chage the seq numbers within brackets to match the min and max lenght desired
#seqtab1 <- seqtab1[,nchar(colnames(seqtab1)) %in% seq(371,383)]

# Inspect distribution of sequence lengths after dropping the outliers
#table(nchar(getSequences(seqtab1)))
#dim(seqtab1) # Gives you info on the number of Samples and ASVs after tail sequence dropping

# Sanity check for run_1
getN <- function(x) sum(getUniques(x))
track1 <- cbind(out1, sapply(dadaFs1, getN))
colnames(track1) <- c("input", "filtered", "denoisedF")
rownames(track1) <- sample.names1
write.csv(track1,"csv/run1_asv_stats.csv") # Export a csv file as archive
head(track1) # check the number of reads retained at each step

#### END OF ANALYSIS BLOCK 1 ##############################################################
## From here on you can either move to ANALYSIS BLOCK 3 if you have a single run to work with,
## or continue to ANALYSIS BLOCK 2 below if you need to analyze and merge multiple
## runs together
###########################################################################################


###########################################################################################
#### ANALYSIS BLOCK 2 #####################################################################
###########################################################################################

## Repeat analysis block 1 for every run you have, changing all the names from 1 to two
## for example run_1 to run_2 for the path, sample.names1 to sample.names2 or seqtab1 to
## seqtab2. ## For simplicity I have already transformed and copied the commands for run_2
## removing the comments below this line. The REPEAT block is highligthed wiht indents.
## Two comments are present only near the command where it is fundamental to interact with the script.

#### REPEAT FROM HERE <--------------------------------------------------------!

path2 <- "/home/ms3468/biofilm_agrosche/dada2/run_13dec23/sequences2"
list.files(path2)
fnFs2 <- sort(list.files(path2, pattern=".fastq", full.names = TRUE))
sample.names2 <- sapply(strsplit(basename(fnFs2), ".fastq"), `[`, 1)
filtFs2 <- file.path(path2, "filtered", paste0(sample.names2, "_F_filt.fastq.gz"))

png("/home/ms3468/biofilm_agrosche/dada2/run_13dec23/plot/QplotF_2.png", width=1600, height=1200)
    plotQualityProfile(fnFs2[1:4])
dev.off()
    
# Select where the sequences should be truncated iteratively if necessary
out2 <- filterAndTrim(fnFs2, filtFs2, truncLen=194, maxN=0, maxEE=4, truncQ=2, rm.phix=TRUE,trimLeft=27, compress=TRUE, multithread=18)
head(out2)
mean(out2[,2]/out2[,1]) #between 0.9-0.95 it is fine. Increase the maxEE if you get less than 0.9 

errF2 <- learnErrors(filtFs2, multithread=12, randomize=TRUE)

png("/home/ms3468/biofilm_agrosche/dada2/run_13dec23/plot/ERR_2.png", width=1600, height=1200)
plotErrors(errF2, nominalQ=TRUE)
dev.off()

derepFs2 <- derepFastq(filtFs2, verbose=TRUE)
names(derepFs2) <- sample.names2
dadaFs2 <- dada(derepFs2, err=errF2, pool="pseudo", multithread=12)
seqtab2 <- makeSequenceTable(dadaFs2)
dim(seqtab2)
    #table(nchar(getSequences(seqtab2)))
    #seqtab2 <- seqtab2[,nchar(colnames(seqtab2)) %in% seq(371,383)] # Drop the tail sequences
    #table(nchar(getSequences(seqtab2)))
    #dim(seqtab2)
track2 <- cbind(out2, sapply(dadaFs2, getN))
colnames(track2) <- c("input", "filtered", "merged")
rownames(track2) <- sample.names2
write.csv(track2,"csv/run2_asv_stats.csv") # Export a csv file as archive
head(track2) # check the number of reads retained at each step for run_2
#### TO HERE <------------------------------------------------------------------!
## Repeat As many time as necessary changing the number in the names. The follow the lines below

## Merging of multiple runs. Add as many seqtab.nochim with the approapriate number
## separated by a comma as needed.
seqtab <- mergeSequenceTables(seqtab1, seqtab2) # Merging of different runs
dim(seqtab)
#### END OF ANALYSIS BLOCK 2 ##############################################################
## Continue the analysis on ANALYSIS BLOCK 3 below
###########################################################################################


###########################################################################################
#### ANALYSIS BLOCK 3 #####################################################################
###########################################################################################

## From here on we refer to the seqtab.nochim object otained from the merger done in
## ANALYSIS BLOCK 2. Change the name to seqtab.nochim1 if you have a single run or
## feel free to renave you objest to seqtab.nochim to avoid changing the code below

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=15)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab) # Gives you the percentage of sequences recovered after chimera removal
## This is a ggod check points. Even if a lot of ASVs have been removed, the majority of reads sould
## be retained. Usually >0.80 (aka 80%) are retained

## Assign Taxonomy. Point to where the silva database actually is
taxa <- assignTaxonomy(seqtab.nochim,"/home/ms3468/silva_db_16S/silva_nr99_v138.1_train_set.fa.gz", multithread=20)
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
write.csv(taxa.print, "csv/taxa_print.csv") # For inspection in bash or excel

saveRDS(seqtab.nochim,"/home/ms3468/biofilm_agrosche/dada2/run_13dec23/rds/seqtabnochim_run13dec23.rds")
saveRDS(taxa,"/home/ms3468/biofilm_agrosche/dada2/run_13dec23/rds/taxa_run13dec23.rds")
save.image("/home/ms3468/biofilm_agrosche/dada2/run_13dec23/rdata/taxa_assigned_run13dec23.RData")

####################################################################################################################

library(DECIPHER) # Multiple sequence alignment and phylogenetic analysis
library(ape) # importing and handling phylogenetic trees
## Phylogenetic tree building. MSA and Phangorn used in the original publication
## scale quadratically with the number of ASV, and became quickly unisable. I have moved to
## AlignSeqs in the Dechipher package and FastTree in bash for making the phylogenetic tree
## Before continuing be sure that the fasttree execuatble is installed in your system and
## visible in your PATH

seqs <- getSequences(seqtab.nochim) ## Get the sequences
names(seqs) <- seqs
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA) # Generating the multiple sequence alignment
writeXStringSet(alignment, "tree/alignment_13dec23.fasta", format="fasta") # Exporting the alignment to an external file

## Build the tree using FastTree in bash using the GTR model. See more detail about
## FastTree at http://www.microbesonline.org/fasttree/
system('fasttree -gtr -nt tree/alignment_13dec23.fasta > alignment_13dec23.tree', intern=TRUE)

tree <- read.tree(file = "alignment_13dec23.tree") # Reading back the tree into the R environment

## Sanity check. The two numbers should match!
dim(seqtab.nochim)
tree$Nnode+2

## From here you can just keep analyzing your dataset. Take a look at our the phyloseq 
## pipeline or get the analysis file from the repository of one of our paper on GitHub
################################################################################################
