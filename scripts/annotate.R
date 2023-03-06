library(GenomicFeatures)
library(ChIPseeker)
library(writexl)

# Get command line arguments
args <- commandArgs(TRUE)
gtf <- args[1]
data_source <- args[2]
organism <- gsub("_", " ", args[3])
infile <- args[4]
contrast <- args[5]

# set outdir prefix for exported files
outdir <- "data/annot_dmrs/"

#---------------#
# make a TxDB object from the annotation gtf file
txdb <- makeTxDbFromGFF(gtf,
    format = "gtf",
    data_source,
    organism)

# read in the region file to be annotated
#   here, SIG diff peaks txt file from diffbind
mydat <- makeGRangesFromDataFrame(read.delim(infile),
    keep.extra.columns = TRUE)
    #seqnames.field = "seqnames",
    #start.field = "start",
    #end.field = "end",
    #ignore.strand = TRUE)

# annotate via Chipseeker
res <- annotatePeak(mydat,
    tssRegion = c(-3000, 3000),
    TxDb = txdb,
    level = "transcript",
    assignGenomicAnnotation = TRUE,
    overlap = "TSS")

# output annotations
write.table(as.data.frame(res),
    paste0(outdir, contrast, ".annot.txt"),
    sep = "\t",
    col.names = TRUE,
    row.names = FALSE,
    quote = FALSE)

# anno pie
pdf(paste0(outdir, contrast, ".anno_pie_chart.pdf"))
plotAnnoPie(res)
dev.off()

# vennpie
pdf(paste0(outdir, contrast, ".venn_pie_chart.pdf"))
#par(mar=c(0.5,2,1,3))
vennpie(res)
dev.off()
