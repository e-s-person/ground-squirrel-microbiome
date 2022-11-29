##### DADA2 Pipeline for California Ground Squirrel Microbiome Analysis #####
#### By Erin Person, adapted from DADA2 Pipeline Tutorial (1.16) by Benjamin Callahan ####
#### Also following workflow by Callahan et al. 2017 on Bioconductor ####
#### https://bioconductor.org/help/course-materials/2017/BioC2017/Day1/Workshops/Microbiome/MicrobiomeWorkflowII.html ####

#### SET-UP ####

# start with an empty workspace
rm(list = ls())
# set seed
set.seed(42)

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
library(decontam)
library(DESeq2)
library(cowplot)

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
filtFs <- file.path(path, "filtered_negcontrol", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered_negcontrol", paste0(sample.names, "_R_filt.fastq.gz"))
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
seqtab <- makeSequenceTable(mergers)
seqtabAll <- makeSequenceTable(mergers[!grepl("Mock", names(mergers))])

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtabAll)))

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
phangAlign <- phangorn::phyDat(as(alignment, "matrix"), type="DNA")
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
samdf <- data.frame(read_csv('microbiome_metadata_w_controls.csv'))
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

#remove un-characterized phyla
#removes 8 taxa
ps1 <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))

# decontam process
df <- as.data.frame(sample_data(ps1)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(ps1)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=sample_or_control)) + geom_point()

# find contaminants at standard threshold (0.05)
sample_data(ps1)$is.neg <- sample_data(ps1)$sample_or_control == "control"
contamdf.prev <- isContaminant(ps1, method="prevalence", neg="is.neg")
table(contamdf.prev$contaminant)

#  identify and remove contaminants
contaminants <- filter(contamdf.prev, contaminant == TRUE)
ps2 <- prune_taxa(!contamdf.prev$contaminant, ps1)

# identify mitochondrial sequences
mitochondrial.seqs <- subset_taxa(ps2, Family == "Mitochondria")
mitochondria.taxa <- colnames(otu_table(mitochondrial.seqs))

# identify non-bacterial sequences
nonbacterial.seqs <- subset_taxa(ps2, Kingdom == "Archaea")
nonbacterial.taxa <- colnames(otu_table(nonbacterial.seqs))

#Function to remove bad taxa with thanks to Sondra Turjeman for code from
#"Comparing invasive and noninvasive faecal sampling in wildlife microbiome studies: A case study on wild common cranes"
pop_taxa = function(physeq, badTaxa){
  allTaxa = taxa_names(physeq)
  allTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(allTaxa, physeq))
}
#Apply function
#removed 2 mitochondrial ASVs and 4 archaea/non-bacterial ASVs
ps2_filtered = pop_taxa(ps2,mitochondria.taxa)
ps2_filtered = pop_taxa(ps2_filtered,nonbacterial.taxa)

# remove negative controls
ps3 <- prune_samples(sample_data(ps2_filtered)$sample_or_control  == "sample", ps2_filtered)


# Compute prevalence of each feature, store as data.frame
prevdf_filtered = apply(X = otu_table(ps3),
               MARGIN = ifelse(taxa_are_rows(ps3), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf_filtered = data.frame(Prevalence = prevdf_filtered,
                    TotalAbundance = taxa_sums(ps3),
                    tax_table(ps3))

phylum_abun_filtered <- prevdf_filtered %>%
  group_by(Phylum) %>%
  summarise(Abundance = sum(TotalAbundance))

table(tax_table(ps3)[, "Phylum"], exclude = NULL)

### rarefaction ###

# rarefy to minimum sampling depth (~8500)
# 62 ASVs removed during rarefaction
ps3_rare <- rarefy_even_depth(ps3, sample.size = min(sample_sums(ps3)), rngseed = 42)

# Compute prevalence of each feature, store as data.frame
prevdf_rare = apply(X = otu_table(ps3_rare),
               MARGIN = ifelse(taxa_are_rows(ps3_rare), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf_rare = data.frame(Prevalence = prevdf_rare,
                    TotalAbundance = taxa_sums(ps3_rare),
                    tax_table(ps3_rare))

# creating abundance column
abun <- sample_sums(ps2_filtered)
samdf$abundance <- abun

#creating shannon diversity column (pre-rarefaction)
Shannon.Diversity <- estimate_richness(ps2_filtered, split = TRUE, measures = c("Shannon"))
samdf$shannon <- Shannon.Diversity$Shannon

# create sample df with no controls
samdf_rare <- samdf %>%
  filter(sample_or_control == "sample")

#creating shannon diversity column for rarefied samples only
Shannon.Diversity.Rare <- estimate_richness(ps3_rare, split = TRUE, measures = c("Shannon"))
samdf_rare$shannon_rare <- Shannon.Diversity.Rare$Shannon

#### STATS ####

### Linear Models ###

## abundance ##

# repeated measures non-parametric test (Friedman test)
abun.friedman <- samdf_rare %>% friedman_test(abundance ~ Treatment |Squirrel.ID)
abun.friedman

## shannon ##
# checking fit of distribution
gamma <- fitdistr(samdf$shannon, densfun = "gamma")
qqp(samdf$shannon, "gamma", shape = gamma$estimate[[1]], rate = gamma$estimate[[2]])

# GLMM with Gamma distribution for non-normally distributed continuous data - with rarefied data
glmm.shannon2 <- glmer(shannon_rare ~ Treatment + Sampling.Week + ln.fecal.cort + body.mass + (1 | Squirrel.ID),
                       family = Gamma(link = "inverse"), data = samdf_rare, glmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
summary(glmm.shannon2)

## Mann-Whitney U post-hoc with rarefaction ##
shannon.mw2 <- wilcox.test(samdf_rare$shannon_rare ~ samdf_rare$Sampling.Week)
shannon.mw2

### PERMANOVAs

### rarefied PERMANOVAS
# Calculate weighted UniFrac distances
wunifrac_dist_r <- phyloseq::distance(ps3_rare, method = "unifrac", weighted=TRUE)

# Calculate unweighted UniFrac distances
unifrac_dist_r <- phyloseq::distance(ps3_rare, method = "unifrac", weighted=FALSE)

# setting blocks (random effects or strata)
perm <- how(nperm = 999)
setBlocks(perm)<-with(samdf_rare, Squirrel.ID)

# Adonis test for Treatment + sampling week weighted, squirrel ID as strata 
perm_w_treat_r <- adonis2(wunifrac_dist_r ~ Treatment + Sampling.Week, data = samdf_rare, by = "terms", permutations = perm)
perm_w_treat_r

# Adonis for squirrel ID alone, weighted, no strata
perm_w_id_r <- adonis2(wunifrac_dist_r ~ Squirrel.ID, data = samdf_rare, permutations = 999)
perm_w_id_r

# Adonis test for Treatment + sampling week unweighted, squirrel ID as strata
perm_u_treat_r <- adonis2(unifrac_dist_r ~ Treatment + Sampling.Week, data = samdf_rare, by = "terms", permutations = perm)
perm_u_treat_r

# Adonis for squirrel ID alone, unweighted, no strata
perm_u_id_r <- adonis2(unifrac_dist_r ~ Squirrel.ID, data = samdf_rare, permutations = 999)
perm_u_id_r

#### DESeq2 ####
alpha = 0.01

## DESeq2 with agglomeration ##

## by genus
ps3_rare_genus <- tax_glom(ps3_rare, "Genus", NArm = FALSE)

DStreat_genus = phyloseq_to_deseq2(ps3_rare_genus, ~ Group)
DStreat_genus = DESeq(DStreat_genus, test="Wald", fitType="parametric")

# pairwise comparisons - ln2 -  no ln2
res.ln2.noln2_genus = results(DStreat_genus, contrast=c("Group", "LN2", "No LN2"), cooksCutoff = FALSE)
sigtab.ln2.noln2_genus = res.ln2.noln2_genus[which(res.ln2.noln2_genus$padj < alpha), ]
sigtab.ln2.noln2_genus = cbind(as(sigtab.ln2.noln2_genus, "data.frame"), as(tax_table(ps3_rare_genus)[rownames(sigtab.ln2.noln2_genus), ], "matrix"))
head(sigtab.ln2.noln2_genus)

# pairwise comparisons - ln2 -  behav
res.ln2.behav_genus = results(DStreat_genus, contrast=c("Group", "LN2", "behav"), cooksCutoff = FALSE)
sigtab.ln2.behav_genus = res.ln2.behav_genus[which(res.ln2.behav_genus$padj < alpha), ]
sigtab.ln2.behav_genus = cbind(as(sigtab.ln2.behav_genus, "data.frame"), as(tax_table(ps3_rare_genus)[rownames(sigtab.ln2.behav_genus), ], "matrix"))
head(sigtab.ln2.behav_genus)

# pairwise comparisons - no ln2 -  behav
res.noln2.behav_genus = results(DStreat_genus, contrast=c("Group", "No LN2", "behav"), cooksCutoff = FALSE)
sigtab.noln2.behav_genus = res.noln2.behav_genus[which(res.noln2.behav_genus$padj < alpha), ]
sigtab.noln2.behav_genus = cbind(as(sigtab.noln2.behav_genus, "data.frame"), as(tax_table(ps3_rare_genus)[rownames(sigtab.noln2.behav_genus), ], "matrix"))
head(sigtab.noln2.behav_genus)

#### PLOTTING ####

## Figure 1

# Figure 1A. abundance plot for ms - treatment
abunplot_rare <- ggplot (data = samdf_rare, mapping = aes(x = Treatment, y = abundance, fill = Treatment)) +
  geom_text(color = "black", aes(x = Treatment, y=130000, label= "n = 12", family = "serif", size = 12)) +
  geom_boxplot(color = "black") +
  theme_classic() +
  scale_fill_grey() +
  ggtitle("A. Sequence read depth") +
  theme(legend.position="none",
        axis.text = element_text(size=14, color = "black", family = "serif"),
        axis.title = element_text(size=14, family = "serif"),
        plot.title = element_text(size=15, hjust = -0.01, family = "serif", face = "bold")) +
  scale_y_continuous(n.breaks = 6) +
  xlab ("Collection Method") +
  ylab ("Sequence Read Depth")
abunplot_rare

# Figure 1B. diversity box plot for manuscript - Treatment
shanplot_rare <- ggplot (data = samdf_rare, mapping = aes(x = Treatment, y = shannon_rare, fill = Treatment)) +
  geom_boxplot(color = "black") +
  geom_text(color = "black", aes(x = Treatment, y=6, label= "n = 12", family = "serif", size = 12)) +
  theme_classic() +
  scale_fill_grey() +
  ggtitle("B. Shannon diversity index")+
  theme(legend.position="none",
        axis.text = element_text(size=14, color = "black", family = "serif"),
        axis.title = element_text(size=14, family = "serif"),
        plot.title = element_text(size=15, hjust = -0.01, family = "serif", face = "bold")) +
  xlab ("Collection Method") +
  ylab ("Shannon Diversity Index")
shanplot_rare

# creating plot panels
png('treatment_rare.png', width = 10, height = 5, units = 'in', res = 600)
plot_grid(abunplot_rare, shanplot_rare, 
          ncol = 2, nrow = 1)
dev.off()

## Figure 2

# unweighted 
unifrac_r <- ordinate(ps3_rare, method="PCoA", distance = unifrac_dist_r)
# weighted
wunifrac_r <- ordinate(ps3_rare, method="PCoA", distance = wunifrac_dist_r)

# create unweighted axis columns
samdf_rare$Axis.1.u = unifrac_r$vectors[, 1]
samdf_rare$Axis.2.u = unifrac_r$vectors[, 2]

# create unweighted axis columns
samdf_rare$Axis.1.w = wunifrac_r$vectors[, 1]
samdf_rare$Axis.2.w = wunifrac_r$vectors[, 2]

# A. unweighted unifrac by individual ID
pcoaID_u_r <- ggplot(data = samdf_rare, mapping = aes(x=Axis.1.u, y=Axis.2.u, shape = id.num, fill = id.num)) +
  geom_point(size = 4, color = "black") +
  scale_shape_manual(name = "ID", labels = c(1:12), values = rep(21:24, len = 12)) +
  scale_fill_manual(name = "ID", labels = c(1:12), values =rep(c("white","gray","black"),  len = 12)) +
  guides(color=guide_legend(title="ID", ncol = 2), shape = guide_legend(title = "ID", ncol = 2)) +
  theme_classic() +
  theme(axis.text = element_text(size=11, color = "black", family = "serif"),
        axis.title = element_text(size=12, family = "serif"),
        plot.title = element_text(size=12, hjust = -0.01, family = "serif", face = "bold"),
        legend.text = element_text(size=11, family = "serif"),
        legend.title = element_text(size=11, family = "serif"))+
  xlab("Axis 1 [13.4%]") +
  ylab("Axis 2 [7.1%]") +
  annotate(geom="text", x=0.11, y=0.25, label="P < 0.001")+
  ggtitle("A. Individual (unweighted)")
pcoaID_u_r

# B. weighted unifrac by individual ID
pcoaID_w_r <- ggplot(data = samdf_rare, mapping = aes(x=Axis.1.w, y=Axis.2.w, shape = id.num, fill = id.num)) +
  geom_point(size = 4, color = "black") +
  scale_shape_manual(name = "ID", labels = c(1:12), values = rep(21:24, len = 12)) +
  scale_fill_manual(name = "ID", labels = c(1:12), values =rep(c("white","gray","black"),  len = 12)) +
  guides(color=guide_legend(title="ID", ncol = 2), shape = guide_legend(title = "ID", ncol = 2)) +
  theme_classic() +
  theme(axis.text = element_text(size=11, color = "black", family = "serif"),
        axis.title = element_text(size=12, family = "serif"),
        plot.title = element_text(size=12, hjust = -0.01, family = "serif", face = "bold"),
        legend.text = element_text(size=11, family = "serif"),
        legend.title = element_text(size=11, family = "serif"))+
  xlab("Axis 1 [39%]") +
  ylab("Axis 2 [26.1%]") +
  annotate(geom="text", x=0.08, y=0.3, label="P < 0.001")+
  ggtitle("A. Individual (weighted)")
pcoaID_w_r

# C. unweighted unifrac by treatment
pcoatreat_u_r <- ggplot(data = samdf_rare, mapping = aes(x=Axis.1.u, y=Axis.2.u, shape = Treatment, fill = Treatment)) +
  geom_point(size = 4) +
  scale_shape_manual(name = "Collection method", labels = c("Soil LN2", "Soil No LN2", "Tub"), values = 21:23) +
  scale_fill_manual(name = "Collection method", labels = c("Soil LN2", "Soil No LN2", "Tub"), values =c("black","gray","white")) +
  theme_classic() +
  theme(axis.text = element_text(size=11, color = "black", family = "serif"),
        axis.title = element_text(size=12, family = "serif"),
        plot.title = element_text(size=12, hjust = -0.01, family = "serif", face = "bold"),
        legend.text = element_text(size=11, family = "serif"),
        legend.title = element_text(size=11, family = "serif"))+
  xlab("Axis 1 [13.4%]") +
  ylab("Axis 2 [7.1%]") +
  annotate(geom="text", x=0.26, y=0.38, label="P < 0.001")+
  ggtitle("C. Collection method (unweighted)") +
  stat_ellipse(aes(color = Treatment), type = "norm") +
  guides(color="none", shape=guide_legend(title="Collection method")) +  
  scale_color_grey()
pcoatreat_u_r

# D. weighted unifrac by treatment
pcoatreat_w_r <- ggplot(data = samdf_rare, mapping = aes(x=Axis.1.w, y=Axis.2.w, shape = Treatment, fill = Treatment)) +
  geom_point(size = 4) +
  scale_shape_manual(name = "Collection method", labels = c("Soil LN2", "Soil No LN2", "Tub"), values = 21:23) +
  scale_fill_manual(name = "Collection method", labels = c("Soil LN2", "Soil No LN2", "Tub"), values =c("black","gray","white")) +
  theme_classic() +
  theme(axis.text = element_text(size=11, color = "black", family = "serif"),
        axis.title = element_text(size=12, family = "serif"),
        plot.title = element_text(size=12, hjust = -0.01, family = "serif", face = "bold"),
        legend.text = element_text(size=11, family = "serif"),
        legend.title = element_text(size=11, family = "serif"))+
  xlab("Axis 1 [39%]") +
  ylab("Axis 2 [26.1%]") +
  annotate(geom="text", x=0.26, y=0.47, label="P = 0.067")+
  ggtitle("D. Collection method (weighted)")+
  stat_ellipse(aes(color = Treatment), type = "norm") +
  guides(color="none", shape=guide_legend(title="Collection method")) +  
  scale_color_grey()
pcoatreat_w_r

# E. unweighted unifrac by sampling week
pcoaweek_u_r <- ggplot(data = samdf_rare, mapping = aes(x=Axis.1.u, y=Axis.2.u, shape = Sampling.Week, fill = Sampling.Week)) +
  geom_point(size = 4) +
  scale_shape_manual(name = "Sampling week", labels = c("Week 1", "Week 2"), values = 21:22) +
  scale_fill_manual(name = "Sampling week", labels = c("Week 1", "Week 2"), values = c("black","gray")) +
  theme_classic() +
  annotate(geom="text", x=0.25, y=0.35, label="P < 0.001")+
  theme(axis.text = element_text(size=11, color = "black", family = "serif"),
        axis.title = element_text(size=12, family = "serif"),
        plot.title = element_text(size=12, hjust = -0.01, family = "serif", face = "bold"),
        legend.text = element_text(size=11, family = "serif"),
        legend.title = element_text(size=11, family = "serif"))+
  xlab("Axis 1 [13.4%]") +
  ylab("Axis 2 [7.1%]") +
  ggtitle("E. Sampling week (unweighted)") +
  stat_ellipse(aes(color = Sampling.Week), type = "norm") +
  guides(color="none", fill=guide_legend(title="Sampling week")) +  
  scale_color_grey()
pcoaweek_u_r

# F. weighted unifrac by sampling week
pcoaweek_w_r <- ggplot(data = samdf_rare, mapping = aes(x=Axis.1.w, y=Axis.2.w, shape = Sampling.Week, fill = Sampling.Week)) +
  geom_point(size = 4) +
  scale_shape_manual(name = "Sampling week", labels = c("Week 1", "Week 2"), values = 21:22) +
  scale_fill_manual(name = "Sampling week", labels = c("Week 1", "Week 2"), values = c("black","gray")) +
  theme_classic() +
  annotate(geom="text", x=0.21, y=0.38, label="P = 0.067")+
  theme(axis.text = element_text(size=11, color = "black", family = "serif"),
        axis.title = element_text(size=12, family = "serif"),
        plot.title = element_text(size=12, hjust = -0.01, family = "serif", face = "bold"),
        legend.text = element_text(size=11, family = "serif"),
        legend.title = element_text(size=11, family = "serif"))+
  xlab("Axis 1 [39%]") +
  ylab("Axis 2 [26.1%]") +
  ggtitle("F. Sampling week (weighted)") +
  stat_ellipse(aes(color = Sampling.Week), type = "norm") +
  guides(color="none", fill=guide_legend(title="Sampling week")) +  
  scale_color_grey()
pcoaweek_w_r

# plot panels for ms
# all
png('pcoa.rare.png', width = 7.5, height = 8, units = 'in', res = 600)
plot_grid(pcoaID_u_r, pcoaID_w_r, pcoatreat_u_r, pcoatreat_w_r, pcoaweek_u_r, pcoaweek_w_r,
          ncol = 2, nrow = 3)
dev.off()

## Figure 3

# diversity box plot for ms - sampling week
shanplot_week <- ggplot (data = samdf_rare, mapping = aes(x = Sampling.Week, y = shannon_rare, fill = Sampling.Week)) +
  geom_boxplot(color = "black") +
  # geom_text(color = "black", aes(x = Sampling.Week, y=6, label= c("n = 8", "n = 4"), family = "serif", size = 12)) +
  annotate("text", x = "Week 1", y = 6, label = "n = 8", family = "serif", size = 5) +
  annotate("text", x = "Week 2", y = 6, label = "n = 4", family = "serif", size = 5) +
  theme_classic() +
  scale_fill_grey() +
  ggtitle("Shannon diversity index")+
  theme(legend.position="none",
        axis.text = element_text(size=14, color = "black", family = "serif"),
        axis.title = element_text(size=14, family = "serif"),
        plot.title = element_text(size=15, hjust = -0.01, family = "serif", face = "bold")) +
  xlab ("Sampling Week") +
  ylab ("Shannon Diversity Index")
shanplot_week

png('week_rare.png', width = 5, height = 5, units = 'in', res = 600)
shanplot_week
dev.off()
