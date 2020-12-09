#!/usr/bin/env Rscript

options(stringsAsFactors = F)

# Load libraries ----------------------------------------------------------

suppressMessages(require(getopt))
suppressMessages(require(ape))
suppressMessages(require(fastcluster))
suppressMessages(require(parallel))
suppressMessages(require(cluster))

# Load functions ----------------------------------------------------------

ncombos <- function(x) (x * (x - 1)) / 2

upgma <- function (D, method = "average", ...){
  # Using hclust from fastcluster 
  suppressMessages(require(fastcluster))
  
  DD <- as.dist(D)
  hc <- hclust(DD, method = method, ...)
  result <- as.phylo(hc)
  result <- reorder(result, "postorder")
  result
}

dist_2d_to_1d <- function (x, y, n) {
  i <- max(c(x, y))
  j <- min(c(x, y))
  valid <- (i >= 1) & (j >= 1) & (i <= n) & (j <= n)
  k <- (2 * n - j) * (j - 1) / 2 + (i - j)
  k[!valid] <- NA_real_
  k
}

getsets <- function(allvalues, maxsize, optimize = TRUE){
  if(length(allvalues) <= maxsize){
    return(list(allvalues))
  } 
  if(maxsize < 2){
    stop("Error: maxsize must be greater than 1")
  }
  nsets <- function(l, x){
    if(x < 2) return(Inf)
    nc <- floor((l - x)/(x - 1))
    r <- (l - x) %% (x - 1)
    return(1 + nc + (nc * (nc + 1) * (x - 1)/2) + ifelse(r > 0, ceiling((l - r) / (x - r)), 0))
  }
  
  # Find the optimal maxsize
  if( optimize ){
    if(maxsize > 100000){
      maxvalues <- seq.int(2, maxsize, length.out = 100000)
    } else {
      maxvalues <- 2:maxsize
    }
    nestcomparisons <- maxvalues ^ 2 * sapply(maxvalues, function(x) nsets(length(allvalues), x))
    maxsize <- rev(maxvalues[nestcomparisons == min(nestcomparisons)])[1]
  }
  
  # Set up the search starting points
  t <- ncombos(length(allvalues))   # Total number of unique combinations needed
  i <- 1                            # Starting indices of values for first combination
  j <- 2:maxsize                    #        [i = 'columns', j = 'rows']
  d <- 1                            # Starting dimensional increment for i
  values <- allvalues[c(i, j)]      # Current set of values
  sets <- list(values)              # List of each set of combinations             
  expn <- nsets(length(allvalues), maxsize) # Number of sets required  
  
  while( length(sets) < expn ){
    if( rev(i)[1] + 1 < j[1] ){
      # STEP RIGHT
      if( rev(i)[1] + d <= length(allvalues) ){
        i <- i + d
      } else {
        # This should only happen when d > 1, i.e. we're on the last row,
        # therefore max of j will be the last row index and i does not need to include the last column index
        i <- (i[1] + d):(length(allvalues) - 1)
      }
    } else {
      # STEP DOWN
      i <- 1:(1 + d - 1) # Reset i
      if( rev(j)[1] + maxsize - 1 <= length(allvalues) ){
        # STAY AT CURRENT DIMENSION
        j <- j + maxsize - 1  
      } else {
        # FLATTEN TO FIT REMAINDER EFFICIENTLY
        j <- (j[1] + maxsize - 1):length(allvalues)
        d <- maxsize - length(j)
        i <- i:(i + d - 1)
      }
    }
    values <- allvalues[unique(c(i, j))]
    sets <- c(sets, list(values))
  }
  return(sets)
}

runnparsesets <- function(sets, alignment, model, cores){
  calci <- function(i, j, l) (l - 0.5) * i - l - (i ^ 2)/2 + j
  getindices <- function(labels, alllabels){
    l <- length(alllabels)
    ilis <- lapply(1:(length(labels) - 1), function(li){
      n1 <- which(alllabels == labels[li])
      n2 <- which(alllabels %in% labels[(li + 1):length(labels)])
      sapply(n2, function(j) calci(n1, j, l))
    }) 
    return(unlist(ilis))
  } 

  allvalues <- names(alignment)
  bf <- base.freq(alignment)

  distout <- mclapply(sets, function(set){
    ds <- dist.dna(alignment[set], model = model, pairwise.deletion = T, base.freq = bf)
    return(cbind(getindices(attr(ds, "Labels"), allvalues), ds))
  }, mc.cores = cores)
  
  distout <- do.call('rbind', distout)
  distout <- distout[!duplicated(distout[,1]), ]
  distout <- distout[order(distout[,1]),2]
  return(makedist(distout, labels = allvalues, method = model))
}

dist_1d_to_2d <- function (k, dist_obj) {
  if (!inherits(dist_obj, "dist")) stop("please provide a 'dist' object")
  n <- attr(dist_obj, "Size")
  valid <- (k >= 1) & (k <= n * (n - 1) / 2)
  k_valid <- k[valid]
  j <- rep.int(NA_real_, length(k))
  j[valid] <- floor(((2 * n + 1) - sqrt((2 * n - 1) ^ 2 - 8 * (k_valid - 1))) / 2)
  i <- j + k - (2 * n - j) * (j - 1) / 2
  c(i, j)
}

makedist <- function(x, size = NULL, labels = NULL, method = NULL){
  if( is.null(size) & is.null(labels) ){
    stop("requires at least one of size or labels")
  } else if ( ! is.null(labels) ){
    attr(x, "Size") <- length(labels)
    attr(x, "Labels") <- labels
  } else if ( ! is.null(size) ){
    attr(x, "Size") <- size
  }
  attr(x, "Diag") <- attr(x, "Upper") <- F
  if( ! is.null(method) ){
    attr(x, "method") <- method
  }
  class(x) <- "dist"
  return(x)
}

subsetdist <- function(dist, subset){
  if( is.character(subset) ){
    mat.idx <- match(subset, attr(dist, "Labels"))
    labels <- subset
    size <- NULL
  } else if ( is.logical(subset) ){
    mat.idx <- which(subset)
    labels <- NULL
    size <- sum(subset)
  } else if ( is.numeric(subset) && all(subset == floor(subset)) ){
    mat.idx <- subset
    labels <- NULL
    size <- length(subset)
  } else {
    stop("subset should be character, logical or interger numeric vector")
  }
  size <- attr(dist, "Size")
  dist.idx <- unlist(sapply(1:(length(mat.idx)-1), function(x){
    sapply( (x+1):length(mat.idx), function(y){
      dist_2d_to_1d(mat.idx[x], mat.idx[y], size)
    })
  }))
  return(makedist(dist[dist.idx], size = size, labels = labels, 
                  method = attr(dist, "method")))
}

upgma_partial_withdistmat <- function(distmat, distmax, cores, startsize){
  # Cluster the sequences
  suppressMessages(require(cluster))
  ngroups <- floor(attr(distmat, "Size") / startsize)
  cl <- pam(distmat, ngroups, diss = T)
  
  # Make sure the largest cluster is not larger than distmax
  while( max(table(cl$clustering)) > distmax ){
    ngroups <- floor(ngroups + 0.25 * ngroups)
    cl <- pam(distmat, ngroups, diss = T)
  }
  
  # Build trees for each cluster
  grouptrees <- mclapply(1:ngroups, function(g){
    contents <- names(cl$clustering)[cl$clustering == g]
    distmat.sub <- subsetdist(distmat, contents)
    upgma(distmat.sub)
  }, mc.cores = cores)
  names(grouptrees) <- cl$medoids
  
  # Build a tree of the centroids
  distmat.centroids <- subsetdist(distmat, cl$medoids)
  tree <- upgma(distmat.centroids)
  
  # Correct the branches of the centroid tree so that bound subtrees will 
  # remain ultrametric with respect to one another
  tip.edges.idx <- match(1:Ntip(tree), tree$edge[,2])
  grouptree.heights <- sapply(tree$tip.label, function(cn){
    node.depth.edgelength(grouptrees[[cn]])[1] })
  corrected.edges <- tree$edge.length[tip.edges.idx] - grouptree.heights
  if( any(corrected.edges < 0) ){
    corrected.edges <- corrected.edges + abs(min(corrected.edges))
  } 
  tree$edge.length[tip.edges.idx] <- corrected.edges
  
  # Bind the subtrees to the centroid tree
  for(cn in tree$tip.label){
    tree <- bind.tree(tree, grouptrees[[cn]], where = which(tree$tip.label == cn))
  }
  
  # This is the maximum height value for which cutting the tree will return 
  # realistic values
  message("Heights on the constructed UPGMA tree are only realistic below ",
          min(grouptree.heights))
  
  return(tree)
}


partial.dist.dna <- function(x, i, j, cores, model){
  j <- j[! j %in% i]
  pdist <- mclapply(data.frame(t(expand.grid(i, j))), function(ij){
    dist.dna(x[ij], model = model, pairwise.deletion = T)
  }, mc.cores = cores)
  return(matrix(unlist(pdist), nrow = length(i), ncol = length(j), dimnames = list(i, j)))
}

greedy_cluster <- function(alignment, clustrad, model, cores, verbose = F){
  start <- Sys.time()
  clusters <- setNames(list(names(alignment)[1]), names(alignment[1]))
  c <- 1
  for(i in 2:length(alignment)){
    n <- names(alignment)[i]
    ndist <- partial.dist.dna(alignment, n, names(clusters), cores, model)
    if( min(ndist) <= clustrad ){
      cent <- colnames(ndist)[which.min(ndist)]
      clusters[[cent]] <- append(clusters[[cent]], n)
    } else {
      clusters[[n]] <- n
      c <- c + 1
    }
    if( verbose ){
      time <- as.numeric(difftime(Sys.time(), start, units = "secs"))
      remain <- round(((length(alignment) - i) * (time / i)) / (60 * 60), 3)
      message(round(time / (60*60), 3), " hours passed, clustered ", i, " sequences into ", c, " clusters, ", remain, " hours left\r", appendLF = F)
    }
  }
  return(clusters)
}

lump_small_clusters <- function(clusters, alignment, minsize, model, cores, verbose){
  sizes <- sapply(clusters, length)
  while( any(sizes < minsize) ){
    toosmall <- names(which.min(sizes))
    close <- partial.dist.dna(alignment, toosmall, names(clusters), cores, model)
    addto <- colnames(close)[which.min(close)]
    clusters[[addto]] <- append(clusters[[addto]], clusters[[toosmall]])
    clusters[[toosmall]] <- NULL
    sizes <- sapply(clusters, length)
  }
  return(clusters)
}

upgma_partial <- function(alignment, distmax, model, cores){
  # Determine the cluster radius
  seq1dist <- partial.dist.dna(alignment, names(alignment)[1], names(alignment)[-1], cores, model)
  seq1dist <- sort(seq1dist[1,])
  clustrad <- unname(round(seq1dist[distmax], 4))
  
  # Do clustering, ensuring largest cluster is not larger than distmax
  clusters <- setNames(list(names(alignment)), names(alignment)[1])
  # Make sure the largest cluster is not larger than distmax
  while( max(sapply(clusters, length)) > distmax ){
    clusters <- greedy_cluster(alignment, clustrad, model, cores, verbose = T)
    clusters <- lump_small_clusters(clusters, alignment, 
                                    minsize =  ceiling(0.1 * max(sapply(clusters, length))),
                                    model, cores)
  }
  
  # Build trees for each cluster
  grouptrees <- mclapply(clusters, function(clus){
    distmat <- dist.dna(alignment[clus],  model = model, pairwise.deletion = T)
    upgma(distmat)
  }, mc.cores = cores)
  
  # Build a tree of the centroids
  distmat.centroids <- dist.dna(alignment[names(clusters)],  model = model, pairwise.deletion = T)
  tree <- upgma(distmat.centroids)
  
  # Correct the branches of the centroid tree so that bound subtrees will 
  # remain ultrametric with respect to one another
  tip.edges.idx <- match(1:Ntip(tree), tree$edge[,2])
  grouptree.heights <- sapply(tree$tip.label, function(cn){
    node.depth.edgelength(grouptrees[[cn]])[1] })
  corrected.edges <- tree$edge.length[tip.edges.idx] - grouptree.heights
  if( any(corrected.edges < 0) ){
    corrected.edges <- corrected.edges + abs(min(corrected.edges))
  } 
  tree$edge.length[tip.edges.idx] <- corrected.edges
  
  # Bind the subtrees to the centroid tree
  for(cn in tree$tip.label){
    tree <- bind.tree(tree, grouptrees[[cn]], where = which(tree$tip.label == cn))
  }
  
  # This is the maximum height value for which cutting the tree will return 
  # realistic values
  message("Heights on the constructed UPGMA tree are only realistic below ",
          min(grouptree.heights))
  
  return(tree)
}


# Set up options ----------------------------------------------------------

spec <- matrix(c(
  'help'     , 'h', 0, "logical",
  'alignment', 'a', 1, "character",
  'model'    , 'm', 2, "character",
  'distmax'  , 'd', 2, "integer",
  'cores'    , 'c', 2, "integer"
), byrow = T, ncol = 4)

# Read options
opt <- getopt(spec)

# Testing
#opt$alignment <-"~/programming/bioinformatics/metamate/tests/data/6_coleoptera_fftnsi.fasta"
#opt$distmax <- 70

# Do help -----------------------------------------------------------------

if ( !is.null(opt$help) ){
  cat(getopt(spec, usage = T))
  q(status = 1)
}

# Set defaults ------------------------------------------------------------

if( is.null(opt$alignment) ){
  stop("Error: path to alignment is required")
}

if ( is.null(opt$model)    ) opt$model <- "F84"
if ( is.null(opt$distmax)  ) opt$distmax <- 65536
if ( is.null(opt$cores)    ) opt$cores <- 1

if( opt$distmax < 2 | opt$distmax > 65536){
  stop("Error: -d/--distmax must be an greater than 1 and less than 65,537")
}

# Load in data ------------------------------------------------------------

alignment <- read.FASTA(opt$alignment)

# Create distance matrix and tree -----------------------------------------

# dist.dna can perform maximum 2^31-1 combinations. The maximum number of 
# sequences generating this number of unique combinations is 65,536 
# ncombos(65536) <= 2^31-1 #TRUE
# ncombos(65537) <= 2^31-1 #FALSE
# The limitation comes from the underlying C implementation of dist.dna
# Instead, we do a hacky workaround whereby we split up the data into a
# set of sets that contain all combinations, compute distances independently
# then rebuild the distance matrix.

if( length(alignment) <= opt$distmax ){
  
  distmat <- dist.dna(alignment, model = opt$model, pairwise.deletion = T)
  distmat[is.nan(distmat)] <- ceiling(max(distmat, na.rm = T))
  tree <- upgma(distmat)
  
} else {
  # if( length(alignment) > 1.8 * opt$distmax ){
  #   message("Warning: the number of sequences in the alignment is very high
  #           relative to the maximum size of individual distance matrix 
  #           computations. Partial distance matrix computation will be 
  #           attempted but this may fail or take a very long time")
  # }
  # sets <- getsets(names(alignment), opt$distmax)
  # distmat <- runnparsesets(sets, alignment, opt$model, opt$cores)
  # tree <- upgma_partial_withdistmat(distmat, opt$distmax, opt$cores, floor(opt$distmax/3))
  
  tree <- upgma_partial(alignment, opt$distmax, opt$model, opt$cores)
  
}

# Write out tree ----------------------------------------------------------

write(write.tree(tree), stdout())
#write.table(as.matrix(distmat), file.path(opt$outdir, paste0(filename, "_distance", ".tsv")))

# Quit --------------------------------------------------------------------

q(status = 0)
