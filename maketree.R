#!/usr/bin/env Rscript

suppressMessages(require(getopt))
suppressMessages(require(ape))
suppressMessages(require(phangorn))

# Set up options
spec <- matrix(c(
  'help'     , 'h', 0, "logical",
  'alignment', 'a', 1, "character",
  'model'    , 'm', 2, "character"
), byrow = T, ncol = 4)

# Read options
opt <- getopt(spec)

# Do help
if ( !is.null(opt$help) ){
  cat(getopt(spec, usage = T))
  q(status = 1)
}

# Set defaults
if( is.null(opt$alignment) ){
  stop("Error: path to alignment is required")
}

if ( is.null(opt$model) ) opt$model = "F84"

# Testing
# opt$alignment <-"gra/5_coleoptera_fftnsi.fa"
# opt$alignment <- "amm/6_coleoptera_fftnsi.fasta"

# Load in data
alignment <- read.FASTA(opt$alignment)

# Create distance matrix
distmat <- dist.dna(alignment, model = opt$model, pairwise.deletion = T)

# Check for overdistance pairs
# if(any(is.nan(distmat))){
#   nancounts <- rowSums(is.nan(as.matrix(distmat)))
#   pc_overdistance <- round((sum(nancounts > 0.8 * length(alignment)) / length(alignment)) * 100, 0)
#   warning(paste(paste0(pc_overdistance, "%"), "of haplotypes are too dissimilar to 80% of other sequences to compute distances. Consider re-aligning. Uncomputed distances will be set to 1."))
#   distmat[is.nan(distmat)] <- 1
# }
distmat[is.nan(distmat)] <- ceiling(max(distmat, na.rm = T))

# Create tree
tree <- upgma(distmat)

# Write out data
write(write.tree(tree), stdout())
#write.table(as.matrix(distmat), file.path(opt$outdir, paste0(filename, "_distance", ".tsv")))

# Quit
q(status = 1)
