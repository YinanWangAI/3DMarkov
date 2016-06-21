# Compute steady-state probability (SSD) from Hi-C
library(data.table)
library(dplyr)
library(reshape2)
library(HiTC)
library(Matrix)
library(GenomicRanges)
library(igraph)
library(rARPACK)
library(RcppArmadillo)

load('chrom.size')
load('centromere.pos')
options(scipen=999)
Rcpp::sourceCpp("computeSSD.cpp")

# read Hi-C raw matrix
LoadRawMatrix <- function(raw.matrix.path, chr.size, resolution) {
  # Load Hi-C raw matrix stored in sparse format and transform it to full matrix
  # Arguments: 
  #   raw.matrix.path: file path for sparse matrix, contain Hi-C data for one 
  #   chromosome
  #   chr.size: the length of chromosome
  #   resolution: Hi-C resolution
  # Output:
  #   Hi-C raw matrix
  raw.data <- fread(raw.matrix.path)
  if (ncol(raw.data) != 3) {
    stop('The input matrix is not sparse')
  }
  colnames(raw.data) <- c('start_bin', 'end_bin', 'counts')
  matrix.dim <- ceiling(chr.size / resolution)
  dim.name <- seq(0, matrix.dim - 1) * resolution
  dim.complement <- setdiff(dim.name, 
                            union(raw.data$start_bin, raw.data$end_bin))
  raw.data.complement <- data.table(
    start_bin = dim.complement,
    end_bin = dim.complement,
    counts = 0
  )
  raw.data <- rbind(raw.data, raw.data.complement)
  
  raw.data.rev <- data.table(
    start_bin = raw.data$end_bin,
    end_bin = raw.data$start_bin,
    counts = raw.data$counts
  )
  raw.data <- rbind(raw.data, raw.data.rev)
  raw.data <- raw.data[!duplicated(raw.data), ]
  raw.matrix <- dcast(raw.data, start_bin ~ end_bin, fill = 0)
  rownames(raw.matrix) <- raw.matrix[[1]]
  raw.matrix[, 'start_bin'] <- NULL
  return(raw.matrix)
}

MatrixToHTCobj <- function(contact.matrix, chr, chr.size, resolution) {
  # transform Hi-C matrix to HTC object
  # Arguments: 
  #   contact.matrix: contact matrix
  #   chr: the chromosome of the matrix
  #   chr.size: the length of chromosome
  #   resolution: Hi-C resolution
  # Output:
  #   HTC object of the contact matrix
  if (ncol(contact.matrix) != nrow(contact.matrix)) {
    stop('The matrix is not symmetric')
  }
  gr.obj <- GRanges(
    seqnames = chr,
    ranges = IRanges(seq(0, chr.size, resolution), 
                     seq(0, chr.size, resolution) + resolution),
    name = seq(0, chr.size, resolution)
  )
  htc.mat <- Matrix(as.matrix(contact.matrix))
  rownames(htc.mat) <- seq(0, chr.size, resolution)
  colnames(htc.mat) <- seq(0, chr.size, resolution)
  htc.obj <- new("HTCexp", htc.mat, gr.obj, gr.obj)
  return(htc.obj)
}

# ICE normalization
ICENormalization <- function(htc.obj, max.iter = 200) {
  # ICE normalization of Hi-C matrix
  # Arguments: 
  #   htc.obj: HTC object
  #   max.iter: maximum number of iteration
  # Output:
  #   HTC object of the ICE normalized matrix
  return(normICE(htc.obj, max_iter = max.iter))
}

# shortest-path correction
ShortestPathCorrection <- function(htc.obj) {
  # shortest-path correction of contact matrix
  # Arguments: 
  #   htc.obj: HTC object
  # Output:
  #   HTC object of the shortest-path corrected matrix
  hic.matrix <- as.matrix(intdata(htc.obj))
  hic.dist <- 1 / hic.matrix
  hic.graph.weight <- melt(hic.dist, id.vars = 'region')
  hic.edges <- rbind(as.numeric(hic.graph.weight[[1]]), 
                     as.numeric(as.character(hic.graph.weight[[2]])))
  hic.graph <- add_edges(make_empty_graph(nrow(hic.dist)),
                         hic.edges / resolution + 1,
                         weight = hic.graph.weight[[3]])
  hic.shortest.path <- distances(hic.graph, algorithm = 'dijkstra')
  colnames(hic.shortest.path) <- colnames(hic.matrix)
  rownames(hic.shortest.path) <- colnames(hic.matrix)
  return(hic.shortest.path)
}

# remove centromere regions and regions without reads
RemoveRegions <- function(hic.matrix, centromere.region, unmapped.region, 
                          resolution, n = 10) {
  # Remove regions (regions without reads mapped, centromere, regions near centromere)
  # Arguments: 
  #   hic.matrix: contact matrix
  #   centromere.region: centromere region
  #   unmapped.region: regions without mapped reads
  #   resolution: resolution of Hi-C
  #   n: the number of bins to be removed around centromere
  # Output:
  #   Hi-C matrix 
  near.centromere.region <- c(centromere.region[1] - n * resolution,
                              centromere.region[2] + n * resolution)
  rm.region <- seq(floor(near.centromere.region[1] / resolution) * resolution, 
                   floor(near.centromere.region[2] / resolution) * resolution, 
                   by = resolution)
  rm.region <- union(rm.region, unmapped.region)
  remain.region.row <- setdiff(rownames(hic.matrix), rm.region)
  remain.region.col <- setdiff(colnames(hic.matrix), rm.region)
  return(hic.matrix[remain.region.row, remain.region.col])
}

# estimate transition matrix
EstimateTransitionMat <- function(ice.matrix, distance.matrix, 
                                  method = 'binomial', 
                                  t = 1e-3, distance_scale = 1) {
  # estimate transition matrix from distance matrix
  # Arguments: 
  #   ice.matrix: ICE normalized matrix
  #   distance.matrix: shortest-path corrected matrix
  #   method: c('binomial', 'Wiener'), binomial or Wiener Process method to 
  #   estimate transition matrix, in the previous manuscript we use binomial method
  # Output:
  #   transition matrix
  if (method != 'binomial' & method != 'Wiener') {
    stop('Please choose the right method')
  }
  if (method == 'binomial') {
    diag.p <- diag(ice.matrix) / rowSums(ice.matrix)
    diag.p <- diag.p[!is.na(diag.p)]
    p.raw <- (1 / distance.matrix) ^ 3
    diag(p.raw) <- 0
    bin.mapped.reads <- rowSums(p.raw)
    transition.mat <- (p.raw / bin.mapped.reads) * (1 - diag.p)
    diag(transition.mat) <- diag.p
  } else if (method == 'Wiener') {
    # p(d) = {1 / [(4 * pi) ^ (3 / 2)]} * exp(-squre(d) / 4)
#     c <- 1
#     d <- 1000
    # transition.mat <- 1 / (4 * pi * c) ^ 1.5 * exp(-(distance.matrix ^ 2) / (4 * c))
    transition.mat <- exp(-((distance.matrix * distance_scale) ^ 2) / (2 * t))
    transition.mat <- transition.mat / rowSums(transition.mat)
  }
  return(transition.mat)
}

# compute SSD
ComputeSSD <- function(transition.mat, iter.epsilon = 1e-8, iter.max = 1e5) {
  # compute steady-state distribution from transition matrix
  # Arguments: 
  #   transition.mat: transition matrix
  # Output:
  #   steady-state distribution
  ssd <- Re(eigs(t(transition.mat), 1)$vector)
  if (length(ssd) == 0 | any(ssd < 0)) {
    ssd <- as.vector(computeSSD(transition.mat, iter.epsilon, iter.max))
  }  # compute ssd through iteration
  ssd <- ssd / sum(ssd)
  names(ssd) <- colnames(transition.mat)
  return(ssd)
}
