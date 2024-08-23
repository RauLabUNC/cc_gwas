#this is mostly just a modified version of what is here https://github.com/gkeele/miqtl

#take a look at it for more explanations of what each step is doing

library(sparcc)
library(miqtl)
#Run GWAS
setwd("D:/Work/1- Projects/Collaborative Cross/ISO Project/GWAS/miQTL")

genomecache <- "./CC_Genome_Cache_Clean_w_Founders"
genomecache2 = "./CC_Genome_Cache/CC_Genome_Cache"



#This will simulate phenotypes.  It only works without founder lines.

#phenotypes <- sim.CC.data(genomecache = genomecache2, 
 #                         num.lines = 40, 
  #                        num.sim = 2, 
   #                       num.replicates = 2, 
    #                      qtl.effect.size = 0.8)


#So, I've modified this file such that it has a column (strain_clean) that has
#just the strain names, nothing more.

phenotypes=read.csv("Full CC Panel Data_04_16_24.csv")

phenotypes_ctrl=phenotypes[phenotypes$Drug=="Ctrl",]

# QTL scan using ROP

miqtl.rop.scan = scan.h2lmm(genomecache = genomecache,
                            data = phenotypes_ctrl,
                            pheno.id="Strain_Clean",
                            geno.id="Strain_Clean",
                            formula = BW.day.0 ~ 1, #This is how you can incorporate covariates like sex or... other things.
                            use.multi.impute = F,
                            return.allele.effects = T)

genome.plotter.whole(scan.list=list(ROP = miqtl.rop.scan)) #and plot it.

# QTL scan using multiple imputation   - this is another method.  It takes a LOT longer
# and sometimes crashes.  Not sure if its better.

miqtl.mi11.scan <- scan.h2lmm(genomecache = genomecache, 
                              data = phenotypes_ctrl,
                              pheno.id="Strain_Clean",
                              geno.id="Strain_Clean",
                              formula=BW.day.0 ~ 1, 
                              use.multi.impute = TRUE, 
                              num.imp = 11)

genome.plotter.whole(scan.list = list(miqtl.mi11.scan),i)

#Permutation.  This will get us the appropriate GWAS significance cutoff.  We 
# *probably* only have to do this once per condition we try - basically, as long as 
# the number of SNPs and Strains tested don't change and the nature of the data
# doesn't change much either (so ctrl vs ISO, for example), things should be fairly stable.
# but if this is all running in the background on the cluster and its fast enough... might as well
# do it every time :) 

#not sure what the right N should be here.  10 seems low.
#okay, the original one breaks at a point trying to generate a variable that it
#never actually needs to use, so load up the alternative function at the end
#of the file...
permuted_phenotype <- generate.sample.outcomes.matrix2(scan.object = miqtl.rop.scan, 
                                                      method = "permutation", num.samples = 10)

permuted_scans <- run.threshold.scans(sim.threshold.object = permuted_phenotype, 
                                      keep.full.scans=TRUE,
                                      genomecache  = genomecache, 
                                      data = phenotypes_ctrl,
                                      use.multi.impute = FALSE, 
                                      scan.seed = 1)

permute_threshold <- get.gev.thresholds(threshold.scans = permuted_scans, 
                                        percentile = 0.95)

genome.plotter.whole(scan.list=list(ROP = miqtl.rop.scan), 
                     hard.thresholds = permute_threshold)





# READ IN THIS FUNCTION TO DO PERMUTATION
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



