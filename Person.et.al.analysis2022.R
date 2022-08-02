##### DADA2 Pipeline for California Ground Squirrel Microbiome Analysis #####
#### By Erin Person, adapted from DADA2 Pipeline Tutorial (1.16) by Benjamin Callahan ####
#### Also following workflow by Callahan et al. 2017 on Bioconductor ####
#### https://bioconductor.org/help/course-materials/2017/BioC2017/Day1/Workshops/Microbiome/MicrobiomeWorkflowII.html ####

#### SET-UP ####

# start with an empty workspace
rm(list = ls())

# load packages
library(dada2)
library(phyloseq)
library(Biostrings)
library(tidyverse)
library(vegan)
library(lme4)
library(rstatix)
library(DECIPHER)
library(phangorn)
library(car)
library(MASS)
library(optimx)

# set working directory and path
setwd("D:/Microbiome Analysis")
path <- "D:/Microbiome Analysis/Raw.Data"

#### IMPORT DATA ####

# Reads in file names following the format described in the string
# This line creates a list of the forward reads
fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq.gz
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#### PROCESSING IN DADA2 ####

# plot read quality profiles for forward reads for the first two samples
plotQualityProfile(fnFs[1:2])

# plot read quality profiles for reverse reads
plotQualityProfile(fnRs[1:2])

# Create new filtered subdirectory and place filtered files there
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# truncating and filtering reads
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,150),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE
head(out)

# learning error rates (forward and reverse)
errF <- learnErrors(filtFs, multithread=FALSE)
errR <- learnErrors(filtRs, multithread=FALSE)

# dereplication (reduces computation time)
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)

# name the derep-class objects by sample names
names(derepFs) <- sample.names
names(derepFs) <- sample.names

# plot errors
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

# applies core sample inference algorithm to trimmed and filtered data
dadaFs <- dada(derepFs, err=errF, multithread=TRUE, pool = TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE, pool = TRUE)

# merges forward and reverse reads into contigs
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

# Inspect the merger data.frame from the first sample
head(mergers[[1]])

# # construct amplicon sequence variant (ASV) table
# seqtab <- makeSequenceTable(mergers)
seqtabAll <- makeSequenceTable(mergers[!grepl("Mock", names(mergers))])

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

# remove non-target length sequences
seqtab2 <- seqtabAll[,nchar(colnames(seqtabAll)) %in% 250:256]

# remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab2)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab2)

#track reads through pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

# assign taxonomy with Silva db
taxa <- assignTaxonomy(seqtab.nochim, "D:/Microbiome Analysis/Raw.Data/silva_nr99_v138_train_set.fa.gz", multithread=TRUE)

# inspect taxonomic assignments
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

# creating a taxonomy tree in DECIPHER for use in phangorn tree
seqs <- getSequences(seqtab.nochim)
names(seqs) <- seqs # this propagates to the tip labels of the tree
alignment <- DECIPHER::AlignSeqs(Biostrings::DNAStringSet(seqs), anchor=NA,verbose=FALSE)

# use phangorn to construct a tree
phangAlign <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phangAlign)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phangAlign)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
detach("package:phangorn", unload=TRUE)

#### PHYLOSEQ ####

#creating the phyloseq objects
OTU = otu_table(seqtab.nochim, taxa_are_rows = FALSE)
TAX = tax_table(taxa)
PHY = phy_tree(fitGTR$tree)

#reading in metadata
#including rownames to match OTU table
samdf <- data.frame(read_csv('microbiome_metadata.csv'))
rownames(samdf) <- sample_names(OTU)
samdf$id.num <- as.factor(samdf$id.num)

#finishing phyloseq objects
SAM = sample_data(samdf)
ps = phyloseq(OTU, TAX, SAM, PHY)

#creates short names for ASVs instead of long DNA strings
#DNA sequences can be found using "refseq" function (refseq(ps))
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

# Create table, number of features for each phyla
table(tax_table(ps)[, "Phylum"], exclude = NULL)

#remove un-characterized phyla
ps1 <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))

# Compute prevalence of each feature, store as data.frame
prevdf = apply(X = otu_table(ps1),
               MARGIN = ifelse(taxa_are_rows(ps1), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps1),
                    tax_table(ps1))

### same as above but includes unID phyla
# Compute prevalence of each feature, store as data.frame
prevdf_unid = apply(X = otu_table(ps),
               MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf_unid = data.frame(Prevalence = prevdf_unid,
                    TotalAbundance = taxa_sums(ps),
                    tax_table(ps))

#compute total and average prevalences of each phylum
#this determines if there are any phyla that are mostly low-prevalence
prev_table <- plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})

#creating shannon diversity column
Shannon.Diversity <- estimate_richness(ps1, split = TRUE, measures = c("Shannon"))
samdf$shannon <- Shannon.Diversity$Shannon

### creating abundance column
abun <- sample_sums(ps1)
samdf$abundance <- abun

#### STATS ####

### Linear Models ###

## abundance ##

# checking fit of negative binomial distribution
nbinom <- fitdistr(samdf$abundance, densfun = "negative binomial")
qqp(samdf$abundance, "nbinom", size = nbinom$estimate[[1]], mu = nbinom$estimate[[2]])

# generalized linear mixed model with abundance, treatment, sampling week and cort as fixed effects, squirrel ID as random effect
# negative binomial distribution (non-normal discrete data)
glmm.abundance1 <- glmer.nb(abundance ~ Treatment + Sampling.Week + ln.fecal.cort + body.mass + (1 | Squirrel.ID),
                            data = samdf)
summary(glmm.abundance1)

## shannon ##
# checking fit of distribution
gamma <- fitdistr(samdf$shannon, densfun = "gamma")
qqp(samdf$shannon, "gamma", shape = gamma$estimate[[1]], rate = gamma$estimate[[2]])

# GLMM with Gamma distribution for non-normally distributed continuous data
glmm.shannon1 <- glmer(shannon ~ Treatment + Sampling.Week + ln.fecal.cort + body.mass + (1 | Squirrel.ID),
                            family = Gamma(link = "inverse"), data = samdf, glmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
summary(glmm.shannon1)

## Mann-Whitney U post-hoc ##
shannon.mw <- wilcox.test(samdf$shannon ~ samdf$Sampling.Week)
shannon.mw

abun.mw <- wilcox.test(samdf$abundance ~ samdf$Sampling.Week)
abun.mw

### PERMANOVAs

# Calculate weighted UniFrac distances
unifrac_w <- phyloseq::distance(ps1, method = "unifrac", weighted=TRUE)

# Calculate unweighted UniFrac distances
unifrac_uw <- phyloseq::distance(ps1, method = "unifrac", weighted=FALSE)

# setting blocks (random effects or strata)
perm <- how(nperm = 999)
setBlocks(perm)<-with(samdf, Squirrel.ID)

# Adonis test for Treatment + sampling week weighted, squirrel ID as strata 
perm_w_treat <- adonis2(unifrac_w ~ Treatment + Sampling.Week, data = samdf, by = "terms", permutations = perm)
perm_w_treat

# Adonis for squirrel ID alone, weighted, no strata
perm_w_id <- adonis2(unifrac_w ~ Squirrel.ID, data = samdf, permutations = 999)
perm_w_id

# Adonis test for Treatment + sampling week unweighted, squirrel ID as strata
perm_u_treat <- adonis2(unifrac_uw ~ Treatment + Sampling.Week, data = samdf, by = "terms", permutations = perm)
perm_u_treat

# Adonis for squirrel ID alone, unweighted, no strata
perm_u_id <- adonis2(unifrac_uw ~ Squirrel.ID, data = samdf, permutations = 999)
perm_u_id