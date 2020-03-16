#!/usr/bin/env Rscript

suppressMessages(require(getopt))
suppressMessages(require(ape))

# Set up options
spec <- matrix(c(
  'help'   , 'h', 0, "logical",
  'tree'   , 't', 1, "character",
  'height' , 'i', 1, "double"
), byrow = T, ncol = 4)

# Read options
opt <- getopt(spec)

# Do help
if ( !is.null(opt$help) ){
  cat(getopt(spec, usage = T))
  q(status = 1)
}

# Set defaults

# Testing
# opt$tree <- "test/TEN_100subset_align_retree2_maxit0_upgma.nwk
# opt$tree <- "amm/numtdumper/6_coleoptera_fftnsi_UPGMA.nwk"

# opt$height = 0.03

# Load in data
if ( is.null(opt$tree) ){
  input_con <- file("stdin")
  open(input_con)
  tree <- read.tree(input_con)
  close(input_con)
} else {
  tree <- read.tree(opt$tree)
}

# Cut tree
cuts <- cutree(as.hclust.phylo(tree), h = opt$height)

# REMEMBER:
# the upgma function is a wrapper for hclust(method = "average"), converting the output to a phylogeny
# In an ultrametric phylogeny, the branching time of a node is equal to the distance from a node to its tips
# The distance between two tips in an ultrametric phylogeny is 2x the distance from either tip to their MRCA
# Therefore if looking for 3% clusters in a phylogeny, need to find nodes with a branching time of <= 0.15
# In a hclust dendrogram, the height value of a merge (=node) is equal to the dissimilarity between its tips.
# Therefore a hclust dendrogram can be considered to be 2x amplified
# So to to look for 3% clusters in a dendrogram, need to cut at a height of 0.3
# as.hclust.phylo() converts from phylo to hclust to run cutree on a hclust object, hclust$heights = 2x branching.times(phylo)

cuts <- cbind(opt$height, gsub('-', ';', names(cuts)), cuts)

# Write out data
write.table(cuts, stdout(), row.names = F, col.names = F, quote = F, sep = "\t")

# Quit
q(status = 1)
