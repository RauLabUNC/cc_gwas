scan.h2lmm.test <- function (genomecache, data, formula, K = NULL, model = c("additive", 
                                                          "full"), locus.as.fixed = TRUE, return.allele.effects = FALSE, 
          p.value.method = c("LRT", "ANOVA"), use.par = "h2", use.multi.impute = TRUE, 
          num.imp = 11, chr = "all", brute = TRUE, use.fix.par = TRUE, 
          seed = 1, pheno.id = "SUBJECT.NAME", geno.id = "SUBJECT.NAME", 
          weights = NULL, do.augment = FALSE, use.full.null = FALSE, 
          added.data.points = 1, just.these.loci = NULL, print.locus.fit = FALSE, 
          use.progress.bar = TRUE, debug.single.fit = FALSE, ...) 
{
  model <- model[1]
  p.value.method <- p.value.method[1]
  h <- DiploprobReader$new(genomecache)
  founders <- h$getFounders()
  num.founders <- length(founders)
  loci <- h$getLoci()
  if (pheno.id != geno.id & is.null(K)) {
    K <- diag(length(unique(data[, geno.id])))
    rownames(K) <- colnames(K) <- unique(data[, geno.id])
  }
  cache.subjects <- rownames(h$getLocusMatrix(loci[1], model = "additive"))
  data.and.K <- make.processed.data(formula = formula, data = data, 
                                    cache.subjects = cache.subjects, K = K, pheno.id = pheno.id, 
                                    geno.id = geno.id)
  data <- data.and.K$data
  K <- data.and.K$K
  if (!is.null(K)) {
    if (all(dim(K) == 0)) {
      stop("The colnames and rownames of K do not match the specified ID geno.id in the data", 
           call. = FALSE)
    }
  }
  cache.subjects <- unique(as.character(data[, geno.id]))
  if (!is.null(weights)) {
    weights <- weights[as.character(data[, pheno.id])]
  }
  loci.chr <- h$getChromOfLocus(loci)
  if (chr[1] != "all") {
    loci.chr <- h$getChromOfLocus(loci)
    loci <- loci[loci.chr %in% chr]
  }
  if (!is.null(just.these.loci)) {
    loci <- loci[loci %in% just.these.loci]
    loci.chr <- loci.chr[loci %in% just.these.loci]
  }
  augment.indicator <- NULL
  formula.string <- Reduce(paste, deparse(formula))
  null.formula <- make.null.formula(formula = formula, do.augment = do.augment)
  original.n <- nrow(data)
  old.data <- data
  if (do.augment) {
    augment.n <- ifelse(model == "additive", num.founders, 
                        num.founders + choose(num.founders, 2))
    augment.indicator <- c(rep(0, original.n), rep(1, augment.n))
    if (!use.full.null) {
      data <- make.simple.augment.data(data = data, formula = formula, 
                                       augment.n = augment.n)
      data <- data.frame(data, augment.indicator = augment.indicator)
      K <- make.simple.augment.K(K = K, augment.n = augment.n)
    }
    if (use.full.null) {
      no.augment.K <- K
      K <- make.full.null.augment.K(K = no.augment.K, original.n = original.n, 
                                    augment.n = augment.n)
      data <- make.full.null.augment.data(formula = formula, 
                                          data = data, no.augment.K = no.augment.K, use.par = use.par, 
                                          brute = brute, original.n = original.n, augment.n = augment.n, 
                                          weights = weights)
    }
    weights <- make.augment.weights(data = data, weights = weights, 
                                    augment.n = augment.n, added.data.points = added.data.points)
  }
  use.lmer <- check.for.lmer.formula(null.formula)
  if (model == "full" | use.multi.impute) {
    happy.locus.old <- paste(loci[1], "RData", sep = ".")
    happy.locus.new <- paste(gsub("([[:upper:]])", "@\\1", 
                                  loci[1]), "RData", sep = ".")
    locus_path <- paste0(genomecache, "/full/chr", h$getChromOfLocus(loci[1]), 
                         "/")
    if (file.exists(paste0(locus_path, "data"))) {
      locus_path <- paste0(locus_path, "data/")
    }
    if (!file.exists(paste0(locus_path, happy.locus.old)) & 
        !file.exists(paste0(locus_path, happy.locus.new))) {
      stop("Full model probabilities not available in genome cache, only additive ROP can be fit", 
           call. = FALSE)
    }
  }
  if (model == "full" & return.allele.effects) {
    return.allele.effects <- FALSE
    cat("Allele effects from regression models currently only available with additive model\n", 
        "Setting return.allele.effects to FALSE\n")
  }
  if (use.lmer & p.value.method == "ANOVA") {
    cat("ANOVA not currently supported with our implementation of LMER, switching to LRT\n")
    p.value.method <- "LRT"
  }
  else if (p.value.method == "ANOVA" & (!is.null(K) | !locus.as.fixed)) {
    cat("Standard ANOVA F-test not valid with mixed effect model, switching to LRT\n")
    p.value.method <- "LRT"
  }
  if (use.lmer & !is.null(K)) {
    stop("Cannot use LMER sparse random effects AND a non-sparse random effect", 
         call. = FALSE)
  }
  if (use.lmer & !locus.as.fixed) {
    stop("Cannot use LMER sparse random effects AND fit locus effect as random", 
         call. = FALSE)
  }
  if (!use.fix.par & !locus.as.fixed) {
    cat("standard ANOVA F-test not valid with mixed effect model, switching to LRT\n")
    use.fix.par <- TRUE
  }
  if (!locus.as.fixed & p.value.method == "LRT") {
    p.value.method <- "LRT.random.locus"
  }
  if (use.lmer) {
    fit0 <- lmmbylmer(formula = null.formula, data = data, 
                      REML = FALSE, weights = weights)
    fit0.REML <- lmmbylmer(formula = null.formula, data = data, 
                           REML = TRUE, weights = weights)
    fix.par <- NULL
  }
  else {
    if (is.null(K)) {
      fit0 <- lmmbygls(formula = null.formula, data = data, 
                       eigen.K = NULL, K = NULL, pheno.id = pheno.id, 
                       use.par = "h2", fix.par = 0, weights = weights, 
                       brute = brute)
      fit0.REML <- lmmbygls(formula = null.formula, data = data, 
                            eigen.K = NULL, K = NULL, pheno.id = pheno.id, 
                            use.par = "h2.REML", fix.par = 0, weights = weights, 
                            brute = brute)
    }
    else {
      if (pheno.id != geno.id) {
        Z <- model.matrix(process.random.formula(geno.id = geno.id), 
                          data = data)
        K <- crossprod(t(Z), tcrossprod(K, Z))
        rownames(K) <- colnames(K) <- as.character(data[, 
                                                        pheno.id])
      }
      if (!is.null(weights)) {
        J <- weights^(1/2) * t(weights^(1/2) * K)
        eigen.J <- process_eigen_decomposition(eigen.decomp = eigen(J, 
                                                                    symmetric = TRUE))
        fit0 <- lmmbygls(null.formula, data = data, pheno.id = pheno.id, 
                         eigen.K = eigen.J, K = K, use.par = use.par, 
                         weights = weights, brute = brute)
        fit0.REML <- lmmbygls(null.formula, data = data, 
                              pheno.id = pheno.id, eigen.K = eigen.J, K = K, 
                              use.par = "h2.REML", weights = weights, brute = brute)
      }
      else {
        eigen.K <- process_eigen_decomposition(eigen.decomp = eigen(K, 
                                                                    symmetric = TRUE))
        fit0 <- lmmbygls(null.formula, data = data, pheno.id = pheno.id, 
                         eigen.K = eigen.K, K = K, use.par = use.par, 
                         weights = weights, brute = brute)
        fit0.REML <- lmmbygls(null.formula, data = data, 
                              pheno.id = pheno.id, eigen.K = eigen.K, K = K, 
                              use.par = "h2.REML", weights = weights, brute = brute)
      }
    }
    if (use.fix.par) {
      fix.par <- ifelse(locus.as.fixed, fit0$h2, fit0.REML$h2)
    }
    if (!use.fix.par) {
      fix.par <- NULL
    }
  }
  MI.LOD <- MI.p.value <- allele.effects <- NULL
  LOD.vec <- p.vec <- df <- rep(NA, length(loci))
  null.data <- data
  if (return.allele.effects) {
    if (use.multi.impute) {
      allele.effects <- array(NA, dim = c(length(founders), 
                                          length(loci), num.imp), dimnames = list(founders, 
                                                                                  loci, paste0("imp", 1:num.imp)))
    }
    else {
      allele.effects <- matrix(NA, nrow = length(founders), 
                               ncol = length(loci), dimnames = list(founders, 
                                                                    loci))
    }
  }
  impute.map <- data.frame(data[, pheno.id], data[, geno.id])
  names(impute.map) <- c(pheno.id, geno.id)
  non.augment.subjects <- as.character(data[, geno.id])[grep(pattern = "augment", 
                                                             x = as.character(data[, geno.id]), invert = TRUE)]
  y <- data$y
  if (!print.locus.fit) {
    if (use.progress.bar) {
      pb <- txtProgressBar(min = 0, max = length(loci), 
                           style = 3)
    }
  }
  for (i in 1:length(loci)) {
    if (use.multi.impute) {
      if (i == 1) {
        MI.LOD <- MI.p.value <- matrix(NA, nrow = num.imp, 
                                       ncol = length(loci))
      }
      diplotype.prob.matrix <- h$getLocusMatrix(loci[i], 
                                                model = "full", subjects = non.augment.subjects)
      if (do.augment) {
        if (model == "additive") {
          augment.matrix <- matrix(0, nrow = augment.n, 
                                   ncol = choose(augment.n, 2) + augment.n)
          for (k in 1:augment.n) {
            augment.matrix[k, k] <- 1
          }
        }
        if (model == "full") {
          augment.matrix <- diag(augment.n)
        }
        sample.names <- rownames(diplotype.prob.matrix)
        diplotype.prob.matrix <- rbind(diplotype.prob.matrix, 
                                       augment.matrix)
        rownames(diplotype.prob.matrix) <- c(sample.names, 
                                             paste0("augment.obs", 1:augment.n))
      }
      if (locus.as.fixed) {
        fit0.for.mi <- fit0
      }
      else {
        fit0.for.mi <- fit0.REML
      }
      fit1 <- multi.imput.lmmbygls(formula = formula, y = y, 
                                   X.probs = diplotype.prob.matrix, weights = weights, 
                                   locus.as.fixed = locus.as.fixed, return.allele.effects = return.allele.effects, 
                                   model = model, p.value.method = p.value.method, 
                                   founders = founders, pheno.id = pheno.id, num.imp = num.imp, 
                                   use.lmer = use.lmer, impute.map = impute.map, 
                                   use.par = use.par, fix.par = fix.par, fit0 = fit0.for.mi, 
                                   do.augment = do.augment, brute = brute, seed = seed)
      MI.LOD[, i] <- fit1$LOD
      MI.p.value[, i] <- fit1$p.value
      LOD.vec[i] <- median(fit1$LOD)
      p.vec[i] <- median(fit1$p.value)
      if (return.allele.effects) {
        allele.effects[, i, ] <- fit1$allele.effects
      }
    }
    else {
      X <- h$getLocusMatrix(loci[i], model = model, subjects = non.augment.subjects)
      keep.col <- 1:ncol(X)
      if (locus.as.fixed) {
        max.column <- which.max(colSums(X, na.rm = TRUE))[1]
        keep.col <- keep.col[keep.col != max.column]
        X <- X[, keep.col]
      }
      colnames(X) <- gsub(pattern = "/", replacement = ".", 
                          x = colnames(X), fixed = TRUE)
      locus.formula <- make.alt.formula(formula = formula, 
                                        X = X, do.augment = do.augment)
      if (do.augment) {
        X.names <- rownames(X)
        if (model == "additive") {
          X <- rbind(X, 2 * diag(augment.n)[, keep.col])
        }
        if (model == "full") {
          X <- rbind(X, diag(augment.n)[, keep.col])
        }
        rownames(X) <- c(X.names, paste0("augment.obs", 
                                         1:augment.n))
      }
      if (use.lmer) {
        data <- cbind(null.data, X)
        fit1 <- lmmbylmer(formula = locus.formula, data = data, 
                          REML = FALSE, weights = weights)
        LOD.vec[i] <- log10(exp(as.numeric(logLik(fit1)) - 
                                  as.numeric(logLik(fit0))))
        p.vec[i] <- pchisq(q = -2 * (as.numeric(logLik(fit0)) - 
                                       as.numeric(logLik(fit1))), df = length(fixef(fit1)) - 
                             length(fixef(fit0)), lower.tail = FALSE)
      }
      else {
        if (locus.as.fixed) {
          X <- cbind(fit0$x, X)
          fit1 <- lmmbygls(formula = locus.formula, y = y, 
                           X = X, eigen.K = fit0$eigen.K, K = fit0$K, 
                           weights = weights, use.par = "h2", fix.par = fix.par, 
                           M = fit0$M, logDetV = fit0$logDetV, brute = brute)
          LOD.vec[i] <- log10(exp(fit1$logLik - fit0$logLik))
          p.vec[i] <- get.p.value(fit0 = fit0, fit1 = fit1, 
                                  method = p.value.method)
          df[i] <- fit1$rank
          if (return.allele.effects) {
            
            #allele.effects[, i] <- get.allele.effects.from.fixef(fit = fit1, 
            #                                                     founders = founders, allele.in.intercept = founders[max.column])
            effects <- fit1$coefficients[founders]
            names(effects) <- founders
            effects <- effects + 0
            effects[founders[max.column]] <- 0
            allele.effects[, i] <- as.vector(scale(effects, center = TRUE, scale = FALSE))
          }
        }
        else {
          fit1 <- lmmbygls.random(formula = null.formula, 
                                  pheno.id = pheno.id, y = y, X = fit0$x, eigen.K = fit0$eigen.K, 
                                  K = fit0$K, Z = X, weights = weights, use.par = "h2", 
                                  null.h2 = fix.par, brute = brute)
          LOD.vec[i] <- log10(exp(fit1$REML.logLik - 
                                    fit0.REML$REML.logLik))
          p.vec[i] <- get.p.value(fit0 = fit0.REML, fit1 = fit1, 
                                  method = p.value.method)
          df[i] <- 1
          if (return.allele.effects) {
            allele.effects[, i] <- get.allele.effects.from.ranef(fit = fit1, 
                                                                 founders = founders)
          }
        }
      }
    }
    if (debug.single.fit) {
      browser()
    }
    if (print.locus.fit) {
      cat(paste("locus", i, "out of", length(loci)), "\n")
    }
    else {
      if (use.progress.bar) {
        setTxtProgressBar(pb, i)
      }
    }
  }
  names(LOD.vec) <- names(p.vec) <- names(df) <- loci
  output <- list(LOD = LOD.vec, p.value = p.vec, MI.LOD = MI.LOD, 
                 MI.p.value = MI.p.value, df = df, pos = list(Mb = h$getMarkerLocation(loci, 
                                                                                       scale = "Mb"), cM = h$getMarkerLocation(loci, scale = "cM")), 
                 loci = loci, chr = h$getChromOfLocus(loci), fit0 = fit0, 
                 fit0.REML = fit0.REML, allele.effects = allele.effects, 
                 y = fit0$y, formula = formula.string, model.type = model, 
                 p.value.method = p.value.method, impute.map = impute.map, 
                 locus.effect.type = fit1$locus.effect.type)
  if (length(just.these.loci) == 1) {
    output$fit1 <- fit1
  }
  if (pheno.id != geno.id & !is.null(K)) {
    rownames(Z) <- as.character(data[, pheno.id])
    output$Z <- Z
  }
  return(output)
}

make.processed.data <- function(formula, data, cache.subjects, K, pheno.id, geno.id){
  all.variables <- all.vars(formula)
  covariates <- all.variables[-1]
  lh.formula.string <- unlist(strsplit(Reduce(paste, deparse(formula)), split="~"))[1]
  lh.formula.string <- gsub("[[:space:]]", "", lh.formula.string)
  covariates <- unique(c(covariates, pheno.id, geno.id))
  formula.string <- paste(lh.formula.string,
                          paste(covariates, collapse="+"),
                          sep="~")
  data <- model.frame(formula(formula.string), data=data)
  names(data) <- c("y", covariates)
  # Selecting those in both data and cache
  include.subjects <- intersect(unique(as.character(data[,geno.id])), cache.subjects)
  data <- data[as.character(data[,geno.id]) %in% include.subjects,]
  #matching <- match(x=as.character(data[,geno.id]), table=include.subjects)
  #data <- data[matching,]
  if(!is.null(K)){
    ### TODO: further selection based on whether K does not have an individual in cach or data
    K <- K[include.subjects, include.subjects]
  }
  if(length(covariates) > 0){
    covariate.matrix <- matrix(NA, nrow=nrow(data), ncol=length(covariates))
    for(i in 1:length(covariates)){
      if(is.factor(data[,covariates[i]])){
        factor.counts <- table(data[,covariates[i]])
        data[,covariates[i]] <- gdata::reorder.factor(x=data[,covariates[i]], new.order=names(sort(factor.counts[factor.counts != 0], decreasing=TRUE)))
      }
    }
  }
  return(list(data=data, K=K))
}

make.simple.augment.K <- function(K, augment.n){
  if(!is.null(K)){
    original.K.names <- colnames(K)
    K <- as.matrix(Matrix::bdiag(K, diag(augment.n)))
    rownames(K) <- colnames(K) <- c(original.K.names, paste0("augment.obs", 1:augment.n))
  }
  return(K)
}

make.simple.augment.data <- function(data, formula, augment.n){
  real.y.names <- data$SUBJECT.NAME
  all.variables <- all.vars(formula)
  covariates <- all.variables[-1]
  augment.y <- rep(mean(data$y), augment.n)
  augment.y.names <- paste0("augment.obs", 1:augment.n)
  covariate.matrix <- NULL
  if(length(all.variables) > 1){
    covariate.matrix <- matrix(NA, nrow=augment.n, ncol=length(covariates))
    for(i in 1:length(covariates)){
      if(is.factor(data[,covariates[i]])){
        covariate.matrix[,i] <- rep(levels(data[,covariates[i]])[1], augment.n)
      }
      if(is.numeric(data[,covariates[i]])){
        covariate.matrix[,i] <- rep(mean(data[,covariates[i]]), augment.n)
      }
    }
  }
  if(is.null(covariate.matrix)){
    partial.augment.data <- data.frame(augment.y, augment.y.names)
    names(partial.augment.data) <- c("y", "SUBJECT.NAME")
  }
  if(!is.null(covariate.matrix)){
    partial.augment.data <- data.frame(augment.y, covariate.matrix, augment.y.names)
    names(partial.augment.data) <- c("y", covariates, "SUBJECT.NAME")
  }
  data <- rbind(data, partial.augment.data)
  return(data)
}

make.augment.weights <- function(data, weights, augment.n, added.data.points){
  if(added.data.points == augment.n & is.null(weights)){
    weights <- NULL
  }
  else if(added.data.points != augment.n & is.null(weights)){
    weights <- c(rep(1, nrow(data) - augment.n), rep(added.data.points/augment.n, augment.n))
    names(weights) <- as.character(data$SUBJECT.NAME)
  }
  else if(!is.null(weights)){
    weights <- c(weights, rep(added.data.points/augment.n, augment.n))
    names(weights) <- as.character(data$SUBJECT.NAME)
  }
  return(weights)
}

make.full.null.augment.K <- function(K, augment.n, original.n){
  if(!is.null(K)){
    original.K.names <- colnames(K)
    K <- as.matrix(Matrix::bdiag(K, diag(augment.n)))
    K[-(1:original.n), 1:original.n] <- K[1:original.n, -(1:original.n)] <- 0.5
    K[-(1:original.n), -(1:original.n)][K[-(1:original.n), -(1:original.n)] == 0] <- 0.5
    rownames(K) <- colnames(K) <- c(original.K.names, paste0("augment.obs", 1:augment.n))
  }
  return(K)
}

make.full.null.augment.data <- function(formula, data, no.augment.K, use.par, brute,
                                        original.n, augment.n, weights){
  all.variables <- all.vars(formula)
  covariates <- all.variables[-1]
  
  null.formula.no.augment <- make.null.formula(formula=formula, is.augmented=FALSE)
  fit0.no.augment <- lmmbygls(formula=null.formula.no.augment, data=data, covariates=covariates, K=no.augment.K,
                              use.par=use.par, brute=brute, null.test=TRUE)
  set.seed(seed)
  y.null.hat <- predict.lmmbygls(fit0.no.augment=fit0.no.augment, original.n=original.n, augment.n=augment.n, 
                                 covariates=covariates, weights=weights)
  real.y.names <- data$SUBJECT.NAME
  augment.y.names <- paste0("augment.obs", 1:augment.n)
  
  covariate.matrix <- NULL
  if(!is.null(covariates)){
    covariate.matrix <- matrix(NA, nrow=augment.n, ncol=length(covariates))
    for(i in 1:length(covariates)){
      if(is.factor(data[,covariates[i]])){
        covariate.matrix[,i] <- rep(levels(data[,covariates[i]])[1], augment.n)
      }
      if(is.numeric(data[,covariates[i]])){
        covariate.matrix[,i] <- rep(mean(data[,covariates[i]]), augment.n)
      }
    }
  }
  partial.augment.data <- data.frame(y.null.hat, covariate.matrix, augment.y.names)
  names(partial.augment.data) <- c(outcome, covariates, "SUBJECT.NAME")
  data <- rbind(data, partial.augment.data)
  data <- cbind(data, data.frame(augment.indicator=c(rep(0, original.n), rep(1, augment.n))))
  return(data)
}

## Formula manipulation functions
make.null.formula <- function(formula, do.augment){
  this.formula.string <- Reduce(paste, deparse(formula))
  this.formula.string <- paste0("y ~ ", unlist(strsplit(this.formula.string, split="~"))[-1])
  this.formula <- as.formula(ifelse(do.augment, paste0(this.formula.string, " + augment.indicator"), this.formula.string))
  return(this.formula)
}
make.alt.formula <- function(formula, X, do.augment){
  this.formula.string <- Reduce(paste, deparse(formula))
  this.formula.string <- paste0("y ~ ", unlist(strsplit(this.formula.string, split="~"))[-1])
  this.formula.string <- ifelse(do.augment, paste0(this.formula.string, " + augment.indicator"), this.formula.string)
  this.formula <- as.formula(paste(this.formula.string, paste(gsub(pattern="/", replacement=".", x=colnames(X), fixed=TRUE), collapse=" + "), sep=" + "))
  return(this.formula)
}
make.snp.null.formula <- function(formula, condition.loci, X.list, model){
  this.formula.string <- Reduce(paste, deparse(formula))
  this.formula.string <- paste0("y ~ ", unlist(strsplit(this.formula.string, split="~"))[-1])
  if(!is.null(condition.loci)){
    for(i in 1:length(condition.loci)){
      if(model == "additive"){
        this.formula.string <- paste(this.formula.string, paste("cond_SNP", i, sep="_"), sep=" + ")
      }
      else if(model == "full"){
        this.formula.string <- paste(this.formula.string, c(paste("cond_SNP", i, "aa", sep="_"), paste("cond_SNP", i, "Aa", sep="_")), sep=" + ")
      }
    }
  }
  return(as.formula(this.formula.string))
}
make.snp.alt.formula <- function(formula, model){
  this.formula.string <- Reduce(paste, deparse(formula))
  this.formula.string <- paste0("y ~ ", unlist(strsplit(this.formula.string, split="~"))[-1])
  if(model == "additive"){
    this.formula <- as.formula(paste(this.formula.string, "SNP", sep=" + "))
  }
  else if(model == "full"){
    this.formula <- as.formula(paste(this.formula.string, "SNP_aa", "SNP_Aa", sep=" + "))
  }
  return(this.formula)
}
remove.whitespace.formula <- function(formula){
  formula.string <- paste0(Reduce(paste, deparse(formula)))
  formula.string <- gsub("[[:space:]]", "", formula.string)
  return(as.formula(formula.string))
}
check.for.lmer.formula <- function(formula){
  formula.string <- paste0(Reduce(paste, deparse(formula)))
  formula.string <- gsub("[[:space:]]", "", formula.string)
  use.lmer <- grepl(pattern="\\([a-zA-Z0-9\\.]+\\|[a-zA-Z0-9\\.]+\\)", x=formula.string, perl=TRUE)
  return(use.lmer)
}
process.random.formula <- function(geno.id){
  random.formula <- as.formula(paste("~", paste(geno.id, "1", sep=" - ")))
  return(random.formula)
}

process_eigen_decomposition <- function(eigen.decomp, 
                                        tol=1e-6){
  # from Robert, who took it from MASS::mvrnorm()
  if(!all(eigen.decomp$values >= -tol * abs(eigen.decomp$values[1L]))){
    stop("K is not positive definite")
  }
  if(any(eigen.decomp$values < 0)){
    if(any(eigen.decomp$values < -tol)){
      message("Zeroing negative eigenvalues: smallest eigenvalue was ", min(eigen.decomp$values), "\n")
    }
    eigen.decomp$values <- pmax(eigen.decomp$values, 0)
  }
  return(eigen.decomp)
}

replicates.eigen <- function(Z, K) {
  eigen <- eigen(K %*% crossprod(Z,Z ), symmetric=FALSE)
  return(list(values=eigen$values,
              vectors=qr.Q(qr(Z %*% eigen$vectors))))
}
