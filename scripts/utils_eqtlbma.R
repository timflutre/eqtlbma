##  `utils_eqtlbma.R' contains utility functions for the eQtlBma package.
##  Copyright (C) 2013 Timothee Flutre
##
##  This program is free software: you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation, either version 3 of the License, or
##  (at your option) any later version.
##
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License
##  along with this program.  If not, see <http://www.gnu.org/licenses/>.

TransformGeneExpInStdNormal <- function(mat, break.ties.rand=TRUE,
                                        seed=1859){
  ## Transform gene expression levels into a N(0,1) via quantile normalization
  ## (gene by gene, across samples).
  ##
  ## Args:
  ##  mat: matrix with genes in rows and samples in columns
  ##  break.ties.rand: should the ties be broken randomly? yes by default
  ##  seed: to make it reproducible
  ##
  ## Returns:
  ##  a matrix with transformed genes in rows and samples in columns
  stopifnot(is.matrix(mat), is.logical(break.ties.rand))
  if(nrow(mat) < ncol(mat))
    warning("input matrix doesn't seem to have genes in rows and samples in columns")
  if(break.ties.rand && ! is.null(seed))
    set.seed(seed)
  
  mat.qn <- t(apply(mat, 1, function(exp.per.gene){
    if(break.ties.rand){
      idx <- sample(length(exp.per.gene))
      tmp <- qqnorm(exp.per.gene[idx], plot.it=FALSE)$x
      tmp[sort(idx, index.return=TRUE)$ix]
    } else
      qqnorm(exp.per.gene, plot.it=FALSE)$x
  }))
  colnames(mat.qn) <- colnames(mat)
  
  return(mat.qn)
}

RemoveConfoundersFromGeneExp <- function(X, confounders){
  ## Remove a set of confounders (e.g. PCs or PEER factors) from a matrix 
  ## of gene expression levels (linear regression per gene).
  ##
  ## Args:
  ##  X: matrix with samples in rows and genes in columns
  ##     (it will be centered and scaled before PCs are removed)
  ##  confounders: matrix with samples in rows and confounders in columns
  ##
  ## Returns:
  ##  matrix of residuals
  stopifnot(is.matrix(X),
            is.matrix(confounders),
            nrow(X) == nrow(confounders))
  if(nrow(X) > ncol(X))
    warning("input matrix doesn't seem to have samples in rows and genes in columns")
  
  res <- lm.fit(x=confounders, y=scale(X, center=TRUE, scale=TRUE))
  return(t(res$residuals))
}

MakeGrid <- function(grid.type="general", no.het=FALSE){
  ## Make the grid used to compute Bayes Factors.
  ##
  ## Args:
  ##  grid.type: "general" indicates the meta-analysis grid (large),
  ##             otherwise the configuration grid (small) is returned
  ##  no.het: if TRUE, the grid is built without heterogeneity
  ##
  ## Returns:
  ##  matrix with the grid
  stopifnot(is.character(grid.type), is.logical(no.het))
  oma2.plus.phi2 <- c(0.1^2, 0.2^2, 0.4^2, 0.8^2, 1.6^2) # avg eff size
  oma2.over.oma2.plus.phi2 <- c(0, 1/4, 1/2, 3/4, 1) # homogeneity
  if(grid.type != "general"){
    if(no.het){
      oma2.over.oma2.plus.phi2 <- c(1)
    } else
      oma2.over.oma2.plus.phi2 <- c(3/4, 1)
  }
  grid <- matrix(NA, nrow=length(oma2.plus.phi2) *
                 length(oma2.over.oma2.plus.phi2), ncol=2)
  colnames(grid) <- c("phi2", "oma2")
  i <- 1
  for(aes in oma2.plus.phi2){
    for(hom in oma2.over.oma2.plus.phi2){
      grid[i,"phi2"] <- aes * (1 - hom)
      grid[i,"oma2"] <- aes * hom
      i <- i + 1
    }
  }
  return(grid)
}

CalcActivityProbasPerSubgroup <- function(configs, plot.it=FALSE,
                                          main=NULL){
  ## Calculate the probability to be active in s subgroups
  ## by summing configuration probabilities.
  ##
  ## Args:
  ##  configs: data.frame with column "id" and column "proba"
  ##           ("id" should be "1" or "1-3-4", etc)
  ##
  ## Returns:
  ##  list with probas per subgroup
  stopifnot(is.data.frame(configs),
            "id" %in% colnames(configs),
            "proba" %in% colnames(configs),
            is.logical(plot.it))
  if(plot.it && is.null(main))
    main <- "Estimates of configuration probabilities"
  
  ## nb of subgroups
  if("0" %in% configs$id){
    S <- log2(nrow(configs))
  } else
    S <- log2(nrow(configs) + 1)
  
  probas1 <- rep(0, 9)
  probas2 <- list()
  for(s in 1:S){
    probas2[[s]] <- c(NA)
    for(i in 1:nrow(configs)){
      if(configs$id[i] == "0")
        next
      if(length(strsplit(configs$id[i], "-")[[1]]) == s){
        probas1[s] <- probas1[s] + configs$proba[i]
        probas2[[s]] <- c(probas2[[s]], configs$proba[i])
      }
    }
    probas2[[s]] <- probas2[[s]][! is.na(probas2[[s]])]
  }
  
  if(plot.it){
    if("0" %in% configs$id){
      plot(x=0:S, y=log10(c(configs$proba[configs$id == "0"], probas1)),
           xlim=c(0,S), ylim=c(-15.1,1),
           pch=19, yaxt="n", col="blue",
           xlab="nb of subgroups in which an eQTL is active",
           ylab="log10(sum over configuration probabilities)",
           main=main)
    } else{
      plot(x=1:S, y=log10(probas1),
           xlim=c(1,S), ylim=c(-15.1,1),
           pch=19, yaxt="n", col="blue",
           xlab="nb of subgroups in which an eQTL is active",
           ylab="log10(sum over configuration probabilities)",
           main=main)
    }
    axis(side=2, at=log10(c(1e-15, 1e-10, 1e-5, 1e-2, 1)),
         labels=c("1e-15", "1e-10", "1e-5", "1e-2", "1"))
    for(s in 1:(S-1))
      points(x=jitter(rep(s, length(probas2[[s]]))), y=log10(probas2[[s]]), pch=1)
  }
  
  return(list(subgroups=probas1, configs=probas2))
}

Log10WeightedSum <- function(x, weights=NULL){
  ## Compute log_{10}(\sum_i w_i 10^x_i) stably.
  ##
  ## Args:
  ##  x: vector of log10(values)
  ##  weights: optional vector of weights (equal weights if not specified)
  ##
  ## Returns:
  ##  log10(\sum_i weights[i] 10^x[i])
  stopifnot(is.vector(x))
  if(! is.null(weights)){
    stopifnot(is.vector(weights))
  } else
    weights <- rep(1/length(x), length(x))
  max <- max(x)
  max + log10(sum(weights * 10^(x - max)))
}

CalcAsymptoticBayesFactorWakefield <- function(theta.hat, V, W, log10=TRUE){
  ## Calculate the asymptotic Bayes Factor proposed by Wakefield
  ## in Genetic Epidemiology 33:79-86 (2009)
  ##
  ## Notes:
  ##  ABF = (\int_theta p(theta.hat|theta) p(theta) dtheta) / p(theta.hat|theta=0)
  ##
  ## Args:
  ##  theta.hat: MLE of the additive genetic effect
  ##  V: variance of theta.hat
  ##  W: variance of the prior on theta
  ##  log10: return the log10 of the ABF is TRUE
  ##
  ## Returns:
  ##  numeric
  z2 <- theta.hat^2 / V # Wald statistic
  log10.ABF <- 0.5 * log10(V) - 0.5 * log10(V + W) +
    (0.5 * z2 * W / (V + W)) / log(10)
  if(log10)
    return(log10.ABF)
  else
    return(10^log10.ABF)
}

CalcExactBayesFactorServinStephens <- function(G, Y, sigma.a, sigma.d, log10=TRUE){
  ## Calculate the exact Bayes Factor proposed by Servin and Stephens
  ## in PLoS Genetics 3 7 (2007)
  ##
  ## Args:
  ##  G: vector of genotypes
  ##  Y: vector of phenotypes
  ##  sigma.a: variance of the prior on the additive genetic effect
  ##  sigma.d: variance of the prior on the dominance genetic effect
  ##  log10: return the log10 of the BF is TRUE
  ##
  ## Returns:
  ##  numeric
  stopifnot(is.vector(G), is.vector(Y))
  subset <- complete.cases(Y) & complete.cases(G)
  Y <- Y[subset]
  G <- G[subset]
  stopifnot(length(Y) == length(G))
  N <- length(G)
  X <- cbind(rep(1,N), G, G == 1)
  inv.Sigma.B <- diag(c(0, 1/sigma.a^2, 1/sigma.d^2))
  inv.Omega <- inv.Sigma.B + t(X) %*% X
  inv.Omega0 <- N
  tY.Y <- t(Y) %*% Y
  log10.BF <- as.numeric(0.5 * log10(inv.Omega0) -
                         0.5 * log10(det(inv.Omega)) -
                         log10(sigma.a) - log10(sigma.d) -
                         (N/2) * (log10(tY.Y - t(Y) %*% X %*% solve(inv.Omega)
                                        %*% t(X) %*% cbind(Y)) -
                                  log10(tY.Y - N*mean(Y)^2)))
  if(log10)
    return(log10.BF)
  else
    return(10^log10.BF)
}
