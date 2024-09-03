#Permutation.  This will get us the appropriate GWAS significance cutoff.  We 
# *probably* only have to do this once per condition we try - basically, as long as 
# the number of SNPs and Strains tested don't change and the nature of the data
# doesn't change much either (so ctrl vs ISO, for example), things should be fairly stable.
# but if this is all running in the background on the cluster and its fast enough... might as well
# do it every time :) 


# Define sample outcomes function
generate.sample.outcomes.matrix2=function (scan.object, model.type = c("null", "alt"), method = c("bootstrap", 
                                                                                                  "permutation", "subsample"), subsample.prop = 0.63, subsample.chr = NULL, 
                                           use.REML = TRUE, use.BLUP = FALSE, num.samples, seed = 1) 
{
  model.type <- model.type[1]
  method <- method[1]
  if (method == "subsample") {
    model.type <- "alt"
    fit <- scan.object$fit0
    n <- floor(subsample.prop * length(fit$y))
    sim.y.matrix <- matrix(NA, nrow = length(fit$y), ncol = num.samples)
    for (i in 1:num.samples) {
      this.subsample <- sample.int(n = length(fit$y), size = n, 
                                   replace = FALSE)
      sim.y.matrix[, i] <- fit$y
      sim.y.matrix[, i][1:length(fit$y) %in% this.subsample] <- NA
    }
    rownames(sim.y.matrix) <- names(fit$y)
    return.weights <- fit$weights
    K <- fit$K
    if (is.null(subsample.chr)) {
      locus <- grab.locus.from.scan(scan.object)
    }
    else {
      locus <- grab.locus.from.scan(scan.object, chr = subsample.chr)
    }
  }
  else {
    if (model.type == "null") {
      fit <- scan.object$fit0
      locus <- NULL
    }
    if (model.type == "alt") {
      if (is.null(scan.object$fit1)) {
        stop("Scan object does not include the alternative model. Re-run scan.h2lmm with the just.these.loci option specifying the locus.", 
             call. = FALSE)
      }
      fit <- scan.object$fit1
      locus <- scan.object$loci
    }
    subsample.rate <- NULL
    fit0.REML <- scan.object$fit0.REML
    if (class(fit) != "lmerMod") {
      na.coefficients <- is.na(fit$coefficients)
      Xb <- fit$x[, !na.coefficients, drop = FALSE] %*% 
        fit$coefficients[!na.coefficients, drop = FALSE]
      n <- nrow(fit$x)
      K <- fit$K
      weights <- fit$weights
      return.weights <- weights
      # if (is.null(weights)) {
      #     weights <- rep(1, nrow(K))
      # }
      if (use.REML) {
        if (is.null(K)) {
          tau2 <- 0
          sigma2 <- fit$sigma2.mle * (n/(n - 1))
        }
        else {
          tau2 <- fit0.REML$tau2.mle
          sigma2 <- fit0.REML$sigma2.mle
        }
      }
      else {
        tau2 <- fit$tau2.mle
        sigma2 <- fit$sigma2.mle
      }
      sim.y.matrix <- matrix(NA, nrow = n, ncol = num.samples)
      if (is.null(K)) {
        u <- rep(0, n)
      }
      else {
        original.K <- K
        impute.map <- scan.object$impute.map
        K <- reduce.large.K(large.K = K, impute.map = impute.map)
        if (use.BLUP) {
          X <- fit$x
          Sigma <- original.K * tau2 + diag(1/weights) * 
            sigma2
          inv.Sigma <- solve(Sigma)
          u.BLUP <- (original.K * tau2) %*% inv.Sigma %*% 
            (diag(nrow(original.K)) - X %*% solve(t(X) %*% 
                                                    inv.Sigma %*% X) %*% t(X) %*% inv.Sigma) %*% 
            fit$y
        }
      }
      set.seed(seed)
      for (i in 1:num.samples) {
        if (tau2 != 0) {
          if (use.BLUP) {
            u <- u.BLUP
          }
          else {
            u <- c(mnormt::rmnorm(1, mean = rep(0, nrow(K)), 
                                  varcov = K * tau2))
          }
          names(u) <- unique(impute.map[, 2])
          u <- u[impute.map[, 2]]
        }
        else {
          u <- rep(0, n)
        }
        if (is.null(weights)) {
          e <- rnorm(n = n, mean = 0, sd = sqrt(sigma2))
        }
        else {
          e <- c(mnormt::rmnorm(1, mean = rep(0, n), 
                                varcov = diag(1/weights) * sigma2))
        }
        y.sample <- Xb + u + e
        if (method == "bootstrap") {
          sim.y.matrix[, i] <- y.sample
        }
        if (method == "permutation") {
          perm.y.ranks <- order(y.sample)
          sim.y.matrix[, i] <- fit$y[perm.y.ranks]
        }
      }
      rownames(sim.y.matrix) <- names(fit$y)
    }
    else {
      stop("Need to add lmer-based functionality!!")
    }
  }
  sim.threshold.object <- list(y.matrix = sim.y.matrix, formula = scan.object$formula, 
                               model = scan.object$model.type, weights = return.weights, 
                               K = K, method = method, impute.map = scan.object$impute.map, 
                               locus = locus, subsample.prop = subsample.prop)
  return(sim.threshold.object)
}


run.threshold.scans.brian <- function (sim.threshold.object, outcome.type = c("outcome", "index"), 
          keep.full.scans = TRUE, scan.index = NULL, genomecache, data, 
          use.multi.impute = TRUE, num.imp = 11, chr = "all", just.these.loci = NULL, 
          scan.seed = 1, ...) 
{
  outcome.type <- outcome.type[1]
  y.matrix <- sim.threshold.object$y.matrix
  formula <- sim.threshold.object$formula
  model <- sim.threshold.object$model
  weights <- sim.threshold.object$weights
  K <- sim.threshold.object$K
  pheno.id <- names(sim.threshold.object$impute.map)[1]
  geno.id <- names(sim.threshold.object$impute.map)[2]
  if (is.null(scan.index)) {
    scan.index <- 1:ncol(y.matrix)
  }
  h <- DiploprobReader$new(genomecache)
  loci <- h$getLoci()
  loci.chr <- h$getChromOfLocus(loci)
  if (chr[1] != "all") {
    loci.chr <- h$getChromOfLocus(loci)
    loci <- loci[loci.chr %in% chr]
  }
  if (!is.null(just.these.loci)) {
    loci <- loci[loci %in% just.these.loci]
    loci.chr <- loci.chr[loci %in% just.these.loci]
  }
  full.p <- full.lod <- these.chr <- these.pos <- NULL
  if (keep.full.scans) {
    full.p <- full.lod <- matrix(NA, nrow = length(scan.index), 
                                 ncol = length(loci))
    colnames(full.p) <- colnames(full.lod) <- loci
    these.chr <- h$getChromOfLocus(loci)
    these.pos <- list(Mb = h$getLocusStart(loci, scale = "Mb"), 
                      cM = h$getLocusStart(loci, scale = "cM"))
  }
  min.p <- max.lod <- rep(NA, length(scan.index))
  iteration.formula <- formula(paste0("new_y ~ ", unlist(strsplit(formula, 
                                                                  split = "~"))[-1]))
  for (i in scan.index) {
    new.y <- data.frame(y.matrix[, i], rownames(y.matrix))
    names(new.y) <- c("new_y", pheno.id)
    this.data <- merge(x = new.y, y = data, by = pheno.id, 
                       all.x = TRUE)
    if (outcome.type == "index") {
      this.data[, all.vars(as.formula(formula))[1]] <- this.data[, 
                                                                 all.vars(as.formula(formula))[1]][this.data$new_y]
      iteration.formula <- formula(formula)
    }
    this.scan <- scan.h2lmm(genomecache = genomecache, data = this.data, 
                            formula = iteration.formula, K = K, model = model, 
                            use.multi.impute = use.multi.impute, num.imp = num.imp, 
                            pheno.id = pheno.id, geno.id = geno.id, seed = scan.seed, 
                            weights = weights, chr = chr, use.progress.bar = T, use.fix.par = T)
    if (keep.full.scans) {
      full.p[i, ] <- this.scan$p.value
      full.lod[i, ] <- this.scan$LOD
    }
    min.p[i] <- min(this.scan$p.value)
    max.lod[i] <- max(this.scan$LOD)
    cat("\n", "Threshold scan:", i, "complete", "\n")
  }
  return(list(full.results = list(LOD = full.lod, p.value = full.p, 
                                  chr = these.chr, pos = these.pos), max.statistics = list(LOD = max.lod, 
                                                                                           p.value = min.p)))
}

# Load libraries
library(miqtl)
library(tidyverse)
# Read in the first trailing argument for phenotype
args <- commandArgs(trailingOnly = TRUE)
phenotype_of_interest <- args[1]


# Load the data
scan_file <- file.path("data/processed/scans", paste0(as.character(phenotype_of_interest), "_scan_results.rds"))
scan <- readRDS(scan_file)

#So, I've modified this file such that it has a column (strain_clean) that has
#just the strain names, nothing more.
phenotypes <- read.csv("data/processed/phenotypes/mean_cc_panel_08_06_24.csv")


phenotypes_ctrl <- phenotypes |> 
  filter(!is.na(get(phenotype_of_interest)) & Drug_Clean == 0)

# Create a function to remove outliers for each given phenotype
remove_outliers <- function(data, variables) {
  for (var in variables) {
    # Calculate the lower and upper bounds for outliers
    q1 <- quantile(data[[var]], 0.25)
    q3 <- quantile(data[[var]], 0.75)
    iqr <- q3 - q1
    lower_bound <- q1 - 1.5 * iqr
    upper_bound <- q3 + 1.5 * iqr
    
    # Remove outliers
    data <- data[!(data[[var]] < lower_bound | data[[var]] > upper_bound), ]
  }
  
  return(data)
}

phenotypes_ctrl <- remove_outliers(phenotypes_ctrl, phenotype_of_interest)

# Load genome cache
genomecache <- "data/raw/genomes/haplotype_cache_cc_083024"

#not sure what the right N should be here.  10 seems low.
#okay, the original one breaks at a point trying to generate a variable that it
#never actually needs to use, so load up the alternative function at the end
#of the file...

permuted_phenotype <- generate.sample.outcomes.matrix2(scan.object = scan, 
                                                      method = "permutation", num.samples = 50,
                                                      use.BLUP = T)
start <- Sys.time()
permuted_scans <- run.threshold.scans.brian(sim.threshold.object = permuted_phenotype, 
                                      keep.full.scans=TRUE,
                                      genomecache  = genomecache, 
                                      data = phenotypes_ctrl,
                                      use.multi.impute = FALSE, 
                                      scan.seed = 1)
end <- Sys.time()
permute_threshold <- get.gev.thresholds(threshold.scans = permuted_scans, 
                                        percentile = 0.85) # 10.1186/s40168-023-01552-8 uses this percentile


# Ensure the directory exists
output_dir <- "data/processed/scan_thresholds"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Define the output file path based on the phenotype of interest
output_threshold <- file.path(output_dir, paste0(as.character(phenotype_of_interest), "_threshold.rds"))
output_scan <- file.path(output_dir, paste0(as.character(phenotype_of_interest), "_scan.rds"))

# Save the threshold and scan as RDS files
saveRDS(permute_threshold, output_threshold)
saveRDS(permuted_scans, output_scan)