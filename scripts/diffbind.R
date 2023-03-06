
library(DiffBind)

# this script only needs to use one contrast
# snakemake will run this for all listed contrasts in config file

# Get command line arguments
args <- commandArgs(TRUE)
bam_files <- strsplit(args[1], ",")[[1]]
peak_files <- strsplit(args[2], ",")[[1]]
design_file <- args[3]
analysis_method <- args[4]
contrast <- args[5]
min_overlap <- args[6]
fdr_cutoff <- args[7]

# Get group names
target <- strsplit(contrast, "_vs_")[[1]][1]
baseline <- strsplit(contrast, "_vs_")[[1]][2]
group_names <- c(target, baseline)

# set outdir prefix for exported files
outdir <- "data/diffbind/"

args[1]
args[2]
args[3]
args[4]
args[5]
args[6]
args[7]

#bam_files
#class(bam_files)
#length(bam_files)
#head(bam_files)
#peak_files
#design_file

#-------------------------#

# Load sample information from design file into sample sheet df
# note: diffBind dba function has set colnames for a sample sheet
sample_sheet <- read.csv(design_file,
    header = TRUE,
    stringsAsFactors = FALSE)

#head(sample_sheet)

# add bam files and peak files to sample sheet data frame
sample_sheet$bamReads <- bam_files
sample_sheet$Peaks <- peak_files
# add PeakCaller column
sample_sheet$PeakCaller <- "narrow"

#head(sample_sheet)

#######################
# 1) Read in Datasets #
#######################

# create dba object
mydat <- dba(sampleSheet = sample_sheet,
    config = data.frame(RunParallel = TRUE,
    reportInit = "DBA",
    DataType = DBA_DATA_GRANGES,
    AnalysisMethod = analysis_method,
    minQCth = 30,
    th = fdr_cutoff,
    bUsePval = "FALSE",
    bRemoveRandom = TRUE,
    bCorPlot = FALSE))

# using peak caller data, produce correlation heatmap using cross-correlations
# of each row in binding matrix
png(paste0(outdir, "corplot.png"))
plot(mydat)
dev.off()

###########################################
# 2) Counting Reads in Consensus Peak Set #
###########################################

# The next step is to calculate a binding matrix with scores based on
#   read counts for every sample (affinity scores), rather than confidence
#   scores for only those peaks called in a specific sample (occupancy scores).
# Also, choose to re-center the peaks at the summit and extend 150bp on
#   either side so all peaks will be the same width.
# By default, minOverlap parameter is set to 2, meaning only include peaks
#   in at least 2 peaksets
mydat <- dba.count(mydat,
    summits = 150,
    minOverlap = min_overlap,
    score = DBA_SCORE_TMM_READS_FULL)

# we can make a new correlation heatmap based on the affinity scores
png(paste0(outdir, "corplot2.png"))
plot(mydat)
dev.off()

# make PCA of the consensus peak counts
png(paste0(outdir, "pca.mydat.png"))
dba.plotPCA(mydat,
    attributes = DBA_CONDITION,
    label = DBA_ID)
dev.off()

#################
# 3) Normalize  #
#################

print("attempt normalization")
mydat_norm <- dba.normalize(mydat,
    #method = analysis_method,
    normalize = DBA_NORM_NATIVE,
    background = TRUE)
print("completed normalization")

###########################
# 4) Establish a Contrast #
###########################

mydat_comp <- dba.contrast(mydat_norm,
    group1 = mydat_norm$masks[[group_names[1]]],
    group2 = mydat_norm$masks[[group_names[2]]],
    name1 = group_names[1],
    name2 = group_names[2])
    #categories = DBA_CONDITION)

mydat_comp

####################################
# 5) Perform Differential Analysis #
####################################

print("begin dba.analyze")
mydat4 <- dba.analyze(mydat_comp)
print("end dba.analyze")
# Can make new correlation heatmap based on the diff bound sites
png(paste0(outdir, contrast, ".corrplot.diff.png"))
plot(mydat4,
    contrast = 1)
    #method = analysis_method)
dev.off()

# Can make PCA plots based on the diff bound sites
png(paste0(outdir, contrast, ".diff.pca.png"))
dba.plotPCA(mydat4,
    #method = analysis_method,
    contrast = 1,
    attributes = DBA_CONDITION,
    label = DBA_ID)
dev.off()

#########################################################
# 6) PValue Distributions of Differentially Bound Sites #
#########################################################

# obtain stats for all binding sites in order to produce
#   pvalue distribution histograms
mydat_db <- dba.report(mydat4,
    th = 1,
    #method = analysis_method,
    bCounts = TRUE)

# make pval distributions
png(paste0(outdir, contrast, ".pval.hist.png"))
hist(as.data.frame(mydat_db)$p.value,
    breaks = 40,
    col = "grey")
dev.off()

#########################
# 7) Export ALL results #
#########################

res1 <- as.data.frame(mydat_db)
res1$PeakID <- paste("Peak", rownames(res1), sep = "_")
write.table(res1,
  paste0(outdir, contrast, ".diff.ALL.results.txt"),
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE)

####################################
# 8) Export SIG Results FDR < 0.05 #
####################################

# Obtain and export sig regions (FDR < 0.05, both gains and losses)
mydat_db_sig <- res1[res1$FDR < as.numeric(fdr_cutoff), ]
write.table(mydat_db_sig,
  paste0(outdir, contrast, ".diff.SIG.results.txt"),
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE)

###########
# Boxplot #
###########

png(paste0(outdir, contrast, ".diff.boxplot.png"))
dba.plotBox(mydat4)
    #method = analysis_method)
dev.off()

###########
# MA Plot #
###########

png(paste0(outdir, contrast, ".diff.MA_plot.png"))
dba.plotMA(mydat4,
    th = as.numeric(fdr_cutoff),
    bUsePval = FALSE,
    #method = analysis_method,
    contrast = 1)
dev.off()
