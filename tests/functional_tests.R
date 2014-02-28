#!/usr/bin/env Rscript

## `functional_tests.R' simulates a typical eQTL data set in multiple
## tissues and analyzes it with R in order to test `bf'.
## Copyright (C) 2012-2013 Timothee Flutre
## License: GPLv3+

## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

#options(warn=2, error=recover)

rm(list=ls())
sink(file=stdout(), type="message")
prog.name <- "functional_tests.R"
prog.version <- "1.0"

R.v.maj <- as.numeric(R.version$major)
R.v.min.1 <- as.numeric(strsplit(R.version$minor, "\\.")[[1]][1])
if(R.v.maj < 2 || (R.v.maj == 2 && R.v.min.1 < 15))
    stop("require R >= 2.15 (for paste0)", call.=FALSE)

##' Display the help on stdout
##'
##' The format complies with help2man (http://www.gnu.org/s/help2man)
##' @title Help
help <- function(){
  txt <- paste0("`", prog.name, "' simulates a typical eQTL data set in multiple tissues and analyzes it with R in order to test `bf'.\n")
  txt <- paste0(txt, "\nUsage: ", prog.name, " [OPTIONS] ...\n")
  txt <- paste0(txt, "\nOptions:\n")
  txt <- paste0(txt, "  -h, --help\tdisplay the help and exit\n")
  txt <- paste0(txt, "  -V, --version\toutput version information and exit\n")
  txt <- paste0(txt, "  -v, --verbose\tverbosity level (0/default=1/2/3)\n")
  txt <- paste0(txt, "      --dir\tpath to the working directory\n")
  txt <- paste0(txt, "      --rgs\tremove some genes from some subgroups\n")
  txt <- paste0(txt, "      --cvrt\tuse covariates\n")
  txt <- paste0(txt, "      --mvlr\tuse multivariate model\n")
  txt <- paste0(txt, "      --ris\tremove some individuals from some subgroups\n")
  txt <- paste0(txt, "      --rgsi\tremove some exp levels in some subgroups for some individuals\n")
  txt <- paste0(txt, "      --lik\tlikelihood (default=norm/pois/qpois)\n")
  txt <- paste0(txt, "      --gfmt\tformat of the genotypes (default=custom/impute)\n")
  message(txt)
}

##' Display version and license information on stdout
##'
##' To comply with help2man (http://www.gnu.org/s/help2man)
##' @title Version
version <- function(){
  txt <- paste0(prog.name, " ", prog.version, "\n")
  txt <- paste0(txt, "\n")
  txt <- paste0(txt, "Copyright (C) 2011-2014 Timothee Flutre.\n")
  txt <- paste0(txt, "License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>\n")
  txt <- paste0(txt, "This is free software; see the source for copying conditions.  There is NO\n")
  txt <- paste0(txt, "warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n")
  txt <- paste0(txt, "\n")
  txt <- paste0(txt, "Written by Timothee Flutre.\n")
  message(txt)
}

##' Parse the command-line arguments
##'
##' Allow short and long options
##' @title Command-line
##' @param params list of parameters initialized with default values
##' @return List of parameters
parseCmdLine <- function(params){
  args <- commandArgs(trailingOnly=TRUE)
  ## print(args)
  
  i <- 0
  while(i < length(args)){ # use "while" loop for options with no argument
    i <- i + 1
    if(args[i] == "-h" || args[i] == "--help"){
      help()
      quit("no", status=0)
    }
    else if(args[i] == "-V" || args[i] == "--version"){
      version()
      quit("no", status=0)
    }
    if(args[i] == "-v" || args[i] == "--verbose"){
      params$verbose <- as.numeric(args[i+1])
      i <- i + 1
    }
    else if(args[i] == "--dir"){
      params$dir.name <- args[i+1]
      i <- i + 1
    }
    else if(args[i] == "--rgs"){
      params$rmvGenesFromSbgrp <- TRUE
      i <- i + 1
    }
    else if(args[i] == "--cvrt"){
      params$withCovars <- TRUE
      i <- i + 1
    }
    else if(args[i] == "--mvlr"){
      params$mvlr <- TRUE
      i <- i + 1
    }
    else if(args[i] == "--ris"){
      params$rmvIndsFromSbgrp <- TRUE
      i <- i + 1
    }
    else if(args[i] == "--rgsi"){
      params$rmvExpFromSbgrpsInds <- TRUE
      i <- i + 1
    }
    else if(args[i] == "--lik"){
      params$lik <- args[i+1]
      i <- i + 1
    }
    else if(args[i] == "--gfmt"){
      params$geno.format <- args[i+1]
      i <- i + 1
    }
  }
  
  if(params$verbose > 1){
    message("parameters:")
    print(params)
  }
  
  if(params$mvlr || params$rmvIndsFromSbgrp)
    suppressPackageStartupMessages(require(MASS)) # for ginv()
  
  return(params)
}

##' Check the values of the command-line parameters
##'
##' @param params list of parameters
checkParams <- function(params){
  stopifnot(! is.null(params$dir.name),
            file.exists(params$dir.name),
            params$lik %in% c("norm","pois","qpois"),
            params$geno.format %in% c("custom","impute"))
}

print.mat2 <- function(mat, a, b){
  for(i in 1:min(a, nrow(mat))){
    for(j in 1:min(b, ncol(mat))){
      cat(sprintf("%.6e  ", mat[i,j]))
    }
    cat("\n")
  }
}

getGrid <- function(grid.type="general", no.het=FALSE){
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

getGrids <- function(){
  grids <- list()
  grids$gridL <- getGrid("general")
  grids$gridS <- getGrid("small")
  return(grids)
}

writeGrids <- function(grids=NULL, verbose=0){
  stopifnot(! is.null(grids))
  if(verbose > 0)
    message("write grids ...")
  write.table(x=grids$gridL, quote=FALSE, sep="\t",
              file=gzfile("grid_phi2_oma2_general.txt.gz"),
              row.names=FALSE, col.names=FALSE)
  write.table(x=grids$gridS, quote=FALSE, sep="\t",
              file=gzfile("grid_phi2_oma2_with-configs.txt.gz"),
              row.names=FALSE, col.names=FALSE)
}

getBinaryConfigs <- function(nb.subgroups=1, verbose=0){
  if(verbose > 0)
    cat("list all configurations for the effect of a SNP in each subgroup\n")
  if(nb.subgroups >= 20)
    warning("nb of subgroups may be too high, better to keep it below 20!",
            call.=FALSE, immediate.=TRUE)
  states <- c(0,1)  # must be sorted!
  if(nb.subgroups == 1)
    matrix(data=states, nrow=2, ncol=1)
  else{
    inner <- Recall(nb.subgroups-1)
    cbind(rep(states, rep(nrow(inner), 2)),
          matrix(t(inner), ncol=ncol(inner), nrow=nrow(inner)*2,
                 byrow=TRUE))
  }
}

getSimulatedData <- function(rmvGenesFromSbgrp=FALSE, rmvIndsFromSbgrp=FALSE,
                             rmvExpFromSbgrpsInds=FALSE, lik="norm",
                             verbose=0){
  if(verbose > 0)
    message("simulate data with R ...")
  
  params <- list()
  params$seed <- 1859
  set.seed(params$seed)
  
  nb.inds <- 200 # don't change, see rmvIndsFromSbgrp
  inds <- data.frame(id=paste0("ind", 1:nb.inds),
                     name=paste0("individual ", 1:nb.inds),
                     stringsAsFactors=FALSE)
  inds$sex <- sample(c(0, 1), nb.inds, replace=TRUE)
  
  nb.pairs <- 14
  params$nb.pairs <- nb.pairs
  null.pairs <- c(c(TRUE,TRUE), c(FALSE,FALSE), c(TRUE,FALSE),
                  c(TRUE,FALSE), c(FALSE,TRUE),
                  TRUE, FALSE, TRUE, FALSE)
  nb.chrs <- 2
  nb.genes <- 10 # 5 with 2 SNP, 4 with 1 SNP, 1 with no SNP
  nb.snps <- 15  # 14 with 1 gene, 1 with no gene
  params$len.cis <- 5 # in bp
  gene.coords <- data.frame(chr=c(rep("chr1", 6),
                              rep("chr2", 4)),
                            start=c(seq(11, 261, length.out=6),
                              seq(11, 161, length.out=4)),
                            end=c(seq(30, 270, length.out=6),
                              seq(30, 180, length.out=4)),
                            id=paste0("gene", 1:nb.genes),
                            score=rep(1000, nb.genes),
                            strand=rep("+", nb.genes),
                            stringsAsFactors=FALSE)
  snp.coords <- data.frame(chr=c(rep("chr1", 11),
                             rep("chr2", 4)),
                           start=c(9,12, 59,62, 109,112, 159,162, 209,212,
                             259, 13, 63, 113, 666),
                           end=NA,
                           id=paste0("snp", 1:nb.snps),
                           stringsAsFactors=FALSE)
  snp.coords$end <- snp.coords$start + 1
  
  maf <- 0.3
  freq.hwe <- function(maf){ # Hardy-Weinberg equilibrium
    f2 <- maf^2 # proba of having 2 copies of minor allele
    f1 <- 2 * maf * (1-maf) # one copy
    f0 <- 1 - f2 - f1 # zero copy
    return(c(f0, f1, f2))
  }
  geno.counts <- matrix(data=replicate(nb.snps, sample(x=0:2, size=nb.inds,
                          replace=TRUE, prob=freq.hwe(maf))),
                        nrow=nb.snps, ncol=nb.inds, byrow=TRUE,
                        dimnames=list(snp=snp.coords$id, ind=inds$id))
  
  nb.subgroups <- 3
  truth <- data.frame(gene=rep("", nb.pairs), snp="", config="", pve.g=0.0,
                      het=0.0, oma=0.0, phi=0.0, bbar=0.0,
                      stringsAsFactors=FALSE)
  for(s in 1:nb.subgroups)
    truth <- cbind(truth, 0.0, 0.0, 0.0)
  colnames(truth)[seq(9,9+3*nb.subgroups-1, 3)] <- paste0("mu",
                                                          1:nb.subgroups)
  colnames(truth)[seq(10,10+3*nb.subgroups-1, 3)] <- paste0("b",
                                                            1:nb.subgroups)
  colnames(truth)[seq(11,11+3*nb.subgroups-1, 3)] <- paste0("sigma",
                                                            1:nb.subgroups)
  subgroups <- data.frame(id=paste0("s", 1:nb.subgroups),
                          name=paste0("tissue ", 1:nb.subgroups),
                          stringsAsFactors=FALSE)
  phenos <- list()
  for(s in 1:nb.subgroups){
    sbgrp.id <- subgroups$id[s]
    phenos[[sbgrp.id]] <- matrix(data=NA, nrow=nb.genes, ncol=nb.inds,
                                 dimnames=list(gene=gene.coords$id,
                                   ind=inds$id))
  }
  configs <- getBinaryConfigs(nb.subgroups)
  pves.g <- runif(n=nb.pairs, min=0.1, max=0.4) # explained by genotype
  hets <- runif(n=nb.pairs, min=0, max=0.2)# heterogeneity
  mus <- rnorm(n=nb.subgroups, mean=4, sd=2)
  sigmas <- rep(1, nb.subgroups) # st dev of errors
  gs.pair <- 0 # index of the gene-SNP pair
  for(g in 1:nb.genes){
    for(p in 1:nb.snps){
      if(snp.coords$chr[p] != gene.coords$chr[g] ||
         snp.coords$start[p] < gene.coords$start[g] - params$len.cis ||
         snp.coords$start[p] > gene.coords$start[g] + params$len.cis)
        next # SNP not in cis
      gs.pair <- gs.pair + 1
      truth[gs.pair,"gene"] <- gene.coords$id[g]
      truth[gs.pair,"snp"] <- snp.coords$id[p]
      if(null.pairs[gs.pair]){ # pair is not an eQTL
        truth[gs.pair,"config"] <- paste(rep("0", nb.subgroups),
                                         collapse="")
        truth[gs.pair,"pve.g"] <- 0
        truth[gs.pair,"het"] <- NA
        truth[gs.pair,"oma"] <- NA
        truth[gs.pair,"phi"] <- NA
        truth[gs.pair,"bbar"] <- NA
        for(s in 1:nb.subgroups){
          truth[gs.pair,paste0("mu", s)] <- mus[s]
          truth[gs.pair,paste0("b", s)] <- 0
          truth[gs.pair,paste0("sigma", s)] <- sigmas[s]
        }
      } else{ # pair is an eQTL
        config <- configs[sample(x=2:nrow(configs), size=1),]
        truth[gs.pair,"config"] <- paste(config, collapse="")
        pve.g <- pves.g[gs.pair]
        truth[gs.pair,"pve.g"] <- pve.g
        het <- hets[gs.pair]
        truth[gs.pair,"het"] <- het
        oma2.plus.phi2 <- pve.g / ((1 - pve.g) * 2 * maf * (1 - maf))
        phi2 <- het * oma2.plus.phi2
        oma2 <- oma2.plus.phi2 - phi2
        truth[gs.pair,"oma"] <- sqrt(oma2)
        truth[gs.pair,"phi"] <- sqrt(phi2)
        bbar <- rnorm(n=1, mean=0, sd=sqrt(oma2))
        truth[gs.pair,"bbar"] <- bbar
        for(s in 1:nb.subgroups){
          truth[gs.pair,paste0("mu", s)] <- mus[s]
          if(config[s] == 0){
            truth[gs.pair,paste0("b", s)] <- 0
          } else
            truth[gs.pair,paste0("b", s)] <- rnorm(n=1, mean=bbar,
                                                   sd=sqrt(phi2))
          truth[gs.pair,paste0("sigma", s)] <- sigmas[s]
        }
      }
    }
  }
  
  for(gs.pair in 1:nb.pairs){
    for(s in 1:nb.subgroups){
      if(length(grep("pois", lik)) == 0){ # Normal likelihood
        phenos[[s]][truth$gene[gs.pair],] <- truth[gs.pair,paste0("mu",s)] +
          truth[gs.pair,paste0("b",s)] * truth[gs.pair,paste0("sigma",s)] *
            geno.counts[truth[gs.pair,"snp"],] +
              rnorm(n=nb.inds, mean=0, sd=truth[gs.pair,paste0("sigma",s)])
      } else{ # Poisson likelihood
        for(i in 1:nb.inds){
          lambda <- exp(truth[gs.pair,paste0("mu",s)] +
                        truth[gs.pair,paste0("b",s)] * truth[gs.pair,paste0("sigma",s)] *
                        geno.counts[truth[gs.pair,"snp"],i])
          phenos[[s]][truth$gene[gs.pair],i] <- rpois(n=1, lambda=lambda)
        }
      }
    }
  }
  
  ## rmv gene9 in s3
  if(rmvGenesFromSbgrp){
    phenos[["s3"]] <- phenos[["s3"]][-9,]
    tmp <- truth
    tmp[tmp$gene == "gene9", "config"] <-
      paste(do.call(c, strsplit(tmp[tmp$gene == "gene9", "config"], ""))[-3],
            collapse="")
    tmp[tmp$gene == "gene9", c("mu3", "b3", "sigma3")] <- NA
    truth <- tmp
  }
  
  ## rmv some individuals from s3
  if(rmvIndsFromSbgrp)
    phenos[["s3"]] <- phenos[["s3"]][,1:100]
  
  ## rmv some exp levels for gene9, s3, ind25
  if(rmvExpFromSbgrpsInds)
    phenos[["s3"]][9,25] <- NA
  
  return(list(inds=inds,
              gene.coords=gene.coords,
              snp.coords=snp.coords,
              geno.counts=geno.counts,
              phenos=phenos,
              truth=truth,
              params=params,
              configs=configs))
}

writeSimulatedData <- function(data=NULL, geno.format="custom",
                               verbose=0){
  stopifnot(! is.null(data))
  
  if(verbose > 0)
    message("write simulated data ...")
  write.table(x=data$gene.coords, file=gzfile("gene_coords.bed.gz"),
              quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
  if(geno.format == "custom"){
    write.table(x=data$snp.coords, file=gzfile("snp_coords.bed.gz"),
                quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
    tmp <- rbind(colnames(data$geno.counts), data$geno.counts)
    tmp <- cbind(c("id", rownames(data$geno.counts)), tmp)
    write.table(x=tmp, file=gzfile("genotypes.txt.gz"),
                quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
    tmp <- data.frame(subgroup=names(data$phenos),
                      file=rep("genotypes.txt.gz",
                        length(data$phenos)))
    write.table(x=tmp, file="list_genotypes.txt", quote=FALSE,
                row.names=FALSE, col.names=FALSE)
  } else if(geno.format == "vcf"){
    message("ERROR: saving genotypes in VCF is not yet implemented")
  } else if(geno.format == "impute"){
    message("ERROR: saving genotypes in IMPUTE is not yet implemented")
    tmp <- data.frame(chr=data$snp.coords$chr,
                      name=data$snp.coords$id,
                      coord=data$snp.coords$end,
                      a1=rep("A", nrow(data$geno.counts)),
                      a2=rep("C", nrow(data$geno.counts)))
    for(i in 1:ncol(data$geno.counts)){
      tmp[[paste0(colnames(data$geno.counts)[i], "_a1a1")]] <-
        ifelse(data$geno.counts[,i] == 0, 1, 0)
      tmp[[paste0(colnames(data$geno.counts)[i], "_a1a2")]] <- 
        ifelse(data$geno.counts[,i] == 1, 1, 0)
      tmp[[paste0(colnames(data$geno.counts)[i], "_a2a2")]] <- 
        ifelse(data$geno.counts[,i] == 2, 1, 0)
    }
    write.table(x=tmp, file=gzfile("genotypes_imp.txt.gz"),
                quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
    tmp <- data.frame(subgroup=names(data$phenos),
                      file=rep("genotypes_imp.txt.gz",
                        length(data$phenos)))
    write.table(x=tmp, file="list_genotypes.txt", quote=FALSE,
                row.names=FALSE, col.names=FALSE)
  }
  
  for(s in 1:length(data$phenos)){
    tmp <- rbind(colnames(data$phenos[[s]]), data$phenos[[s]])
    tmp <- cbind(c("id", rownames(data$phenos[[s]])), tmp)
    write.table(x=tmp, file=gzfile(paste0("phenotypes_",names(data$phenos)[s],
                         ".txt.gz")), quote=FALSE, sep="\t", row.names=FALSE,
                col.names=FALSE)
  }
  tmp <- data.frame(subgroup=names(data$phenos),
                    file=paste0("phenotypes_",names(data$phenos),".txt.gz"))
  write.table(x=tmp, file="list_phenotypes.txt", quote=FALSE,
              row.names=FALSE, col.names=FALSE)
  
  tmp <- rbind(c("id", data$inds$id),
               c("sex", data$inds$sex))
  write.table(x=tmp, file="covariates.txt", quote=FALSE, sep="\t",
              row.names=FALSE, col.names=FALSE)
  tmp <- data.frame(subgroup=names(data$phenos),
                    file=rep("covariates.txt", length(data$phenos)))
  write.table(x=tmp, file="list_covariates.txt", quote=FALSE,
              row.names=FALSE, col.names=FALSE)
  
  write.table(x=data$truth, file=gzfile("truth.txt.gz"),
              quote=FALSE, sep="\t", row.names=FALSE)
}

## Freely available online explanations at:
## http://sfb649.wiwi.hu-berlin.de/fedc_homepage/xplore/ebooks/html/spm/spmhtmlnode27.html
irls <- function(y, X, threshold, family, verbose=0){
  stopifnot(is.matrix(X), is.vector(y), nrow(X) == length(y),
            length(grep("poisson", family)) != 0)
  
  ## internal functions ---------------
  init.mu <- function(y, family){
    if(family != "binomial"){
      mu <- y
      mu[mu == 0] <- 0.01
    } else{
      mu <- (y + 0.5) / (k + 1) # k: nb of trials
    }
    return(mu)
  }
  compute.z <- function(y, mu, family){
    if(length(grep("poisson", family)) != 0){
      return(log(mu) + (y - mu) * (1/mu))
    }
  }
  compute.weights <- function(mu, family){
    if(length(grep("poisson", family)) != 0){
      return(as.vector(mu))
    }
  }
  weighted.leastsquare.fit <- function(X, w, z){
    fit <- lm(z ~ X[,2], weights=w)
    chisq <- sum(w * fit$residuals^2)
    return(list(rank=fit$rank, chisq=chisq, beta.hat=as.numeric(coefficients(fit))))
  }
  compute.mu <- function(beta, X){
    return(exp(X %*% beta))
  }
  compute.variance <- function(X, w, family, y=NULL, beta=NULL, rank=NULL){
    if(length(grep("poisson", family)) != 0){
      if(family == "poisson"){
        return(list(cov=solve(t(X) %*% diag(w) %*% X), sigma2=1))
      } else if(family == "quasipoisson"){
        mu <- as.vector(compute.mu(beta, X))
        sigma2 <- (1/(nrow(X) - rank)) * sum((y - mu)^2 / mu)
        return(list(cov=sigma2 * solve(t(X) %*% diag(w) %*% X), sigma2=sigma2))
      }
    }
  }
  compute.pvalue <- function(family, beta.hat, se.beta.hat, N, rank){
    if(family == "poisson"){
      return(2 * pnorm(-abs(beta.hat / se.beta.hat)))
    } else if(family == "quasipoisson"){
      return(2 * pt(-abs(beta.hat / se.beta.hat), df=N-rank))
    }
  }
  ##-----------------------------------
  
  mu <- init.mu(y, family)
  old.chisq <- -1
  nb.iters <- 0
  
  while(TRUE){
    if(verbose > 0)
      message(nb.iters)
    z <- compute.z(y, mu, family)
    w <- compute.weights(mu, family)
    wls.fit <- weighted.leastsquare.fit(X, w, z)
    if(verbose > 0){
      tmp <- paste0("chisq=", format(wls.fit$chisq, digits=10))
      tmp <- paste0(tmp, " beta.hat=(", format(wls.fit$beta.hat[1], digits=6))
      for(beta in wls.fit$beta.hat[-1])
        tmp <- paste0(tmp, ",", format(beta, digits=6))
      tmp <- paste0(tmp, ")")
      message(tmp)
    }
    if(abs(wls.fit$chisq - old.chisq) < threshold){
      var.beta.hat <- compute.variance(X, w, family, y, wls.fit$beta.hat,
                                       wls.fit$rank)
      pval.beta.hat <- compute.pvalue(family, wls.fit$beta.hat,
                                      sqrt(diag(var.beta.hat$cov)), length(y),
                                      wls.fit$rank)
      break
    }
    old.chisq <- wls.fit$chisq
    mu <- compute.mu(wls.fit$beta.hat, X)
    nb.iters <- nb.iters + 1
  }
  
  return(list(beta.hat=wls.fit$beta.hat,
              se.beta.hat=sqrt(diag(var.beta.hat$cov)),
              pval.beta.hat=pval.beta.hat,
              sigma2=var.beta.hat$sigma2,
              chisq=wls.fit$chisq,
              nb.iters=nb.iters))
}

## Calculate the summary statistics in each tissue
calcSstatsOnSimulatedData <- function(data=NULL, withCovars=FALSE,
                                      lik="norm", verbose=0){
  stopifnot(! is.null(data))
  if(verbose > 0)
    message("calculate the summary statistics in each tissue ...")
  
  sstats <- lapply(data$phenos, function(x){
    tmp <- data.frame(gene=rep(NA, data$params$nb.pairs), snp=NA, maf=NA, n=NA,
                      pve=NA, sigmahat=NA, betahat.geno=NA, sebetahat.geno=NA,
                      betapval.geno=NA, stringsAsFactors=FALSE)
    if(withCovars && "sex" %in% names(data$inds)){
      tmp$betahat.sex <- NA
      tmp$sebetahat.sex <- NA
      tmp$betapval.sex <- NA
    }
    tmp
  })
  
  i <- 0
  for(g in 1:nrow(data$gene.coords)){ # loop over genes
    for(p in 1:nrow(data$snp.coords)){ # loop over snps
      if(data$snp.coords$chr[p] != data$gene.coords$chr[g] ||
         data$snp.coords$start[p] < data$gene.coords$start[g] - data$params$len.cis ||
         data$snp.coords$start[p] > data$gene.coords$start[g] + data$params$len.cis)
        next # SNP not in cis
      if(verbose > 0)
        message(paste(data$gene.coords$id[g], data$snp.coords$id[p]))
      i <- i + 1
      gene <- data$gene.coords$id[g]
      
      for(s in 1:length(data$phenos)){ # loop over subgroups
        snp <- data$snp.coords$id[p]
        sstats[[s]]$gene[i] <- gene
        sstats[[s]]$snp[i] <- snp
        sstats[[s]]$maf[i] <- sum(as.numeric(data$geno.counts[p,])) /
          (2 * length(as.numeric(data$geno.counts[p,])))
        if(gene %in% rownames(data$phenos[[s]])){
          common.inds.idx <- which(data$inds$id %in% colnames(data$phenos[[s]]))
          inds.with.exp.idx <- which(! is.na(data$phenos[[s]][g,])) # handle rmvExpFromSbgrpsInds
          inds.tokeep.idx <- c()
          for(idx in 1:length(data$inds$id))
            if(idx %in% common.inds.idx & idx %in% inds.with.exp.idx)
              inds.tokeep.idx <- append(inds.tokeep.idx, idx)
          sstats[[s]]$n[i] <- length(inds.tokeep.idx)
          
          if(length(grep("pois", lik)) == 0){ # Normal likelihood
            if(! (withCovars && "sex" %in% names(data$inds))){
              tmp <- summary(lm(as.numeric(data$phenos[[s]][g,inds.tokeep.idx]) ~
                                as.numeric(data$geno.counts[p,inds.tokeep.idx])))
            } else{
              tmp <- summary(lm(as.numeric(data$phenos[[s]][g,inds.tokeep.idx]) ~
                                as.numeric(data$geno.counts[p,inds.tokeep.idx]) +
                                as.numeric(data$inds$sex[inds.tokeep.idx])))
            }
            sstats[[s]]$pve[i] <- tmp$r.squared
            sstats[[s]]$sigmahat[i] <- tmp$sigma
            sstats[[s]]$betapval.geno[i] <- tmp$coefficients[2,"Pr(>|t|)"]
            sstats[[s]]$betahat.geno[i] <- tmp$coefficients[2,"Estimate"]
            sstats[[s]]$sebetahat.geno[i] <- tmp$coefficients[2,"Std. Error"]
          } else{
            
            if(lik == "pois"){
              fit <- glm(as.numeric(data$phenos[[s]][g,inds.tokeep.idx]) ~
                         as.numeric(data$geno.counts[p,inds.tokeep.idx]),
                         family=poisson(link="log"))
            } else # lik == "qpois"
              fit <- glm(as.numeric(data$phenos[[s]][g,inds.tokeep.idx]) ~
                         as.numeric(data$geno.counts[p,inds.tokeep.idx]),
                         family=quasipoisson(link="log"))
            tmp <- summary(fit)
            sstats[[s]]$sigmahat[i] <- sqrt(tmp$dispersion) # diff than 1 if qpois
            sstats[[s]]$betahat.geno[i] <- tmp$coefficients[2,"Estimate"]
            sstats[[s]]$sebetahat.geno[i] <- tmp$coefficients[2,"Std. Error"]
            if(lik == "pois"){
              sstats[[s]]$betapval.geno[i] <- tmp$coefficients[2,"Pr(>|z|)"]
            } else
              sstats[[s]]$betapval.geno[i] <- tmp$coefficients[2,"Pr(>|t|)"]
            
            ## if(lik == "pois"){
            ##   tmp <- irls(y=as.numeric(data$phenos[[s]][g,inds.tokeep.idx]),
            ##               X=cbind(rep(1,length(inds.tokeep.idx)), as.numeric(data$geno.counts[p,inds.tokeep.idx])),
            ##               threshold=1e-6, family="poisson")
            ## } else
            ##   tmp <- irls(y=as.numeric(data$phenos[[s]][g,inds.tokeep.idx]),
            ##               X=cbind(rep(1,length(inds.tokeep.idx)), as.numeric(data$geno.counts[p,inds.tokeep.idx])),
            ##               threshold=1e-6, family="quasipoisson",
            ##               verbose=ifelse(gene == "gene1" & snp == "snp1" & s == 1, 1, 0))
            ## sstats[[s]]$sigmahat[i] <- sqrt(tmp$sigma2)
            ## sstats[[s]]$betahat.geno[i] <- tmp$beta.hat[2]
            ## sstats[[s]]$sebetahat.geno[i] <- tmp$se.beta.hat[2]
            ## sstats[[s]]$betapval.geno[i] <- tmp$pval.beta.hat[2]
            
          } # end if pois or qpois
          
        } else # case where gene is absent in subgroup
          sstats[[s]]$n[i] <- 0
        
      } # end of for loop over subgroups
    } # end of for loop over snps
  } # end of for loop over genes
  
  return(sstats)
}

getStdSstatsAndCorrSmallSampleSize <- function(data, sstats, g, p,
                                               nbCovars,
                                               correct=TRUE){
  std.sstats.corr <- do.call(rbind, lapply(sstats, function(x){
    N <- x[x$gene == data$gene.coords$id[g] &
           x$snp == data$snp.coords$id[p],
           "n"]
    if(N == 0){
      c(NA, NA, NA)
    } else{
      sigmahat <- x[x$gene == data$gene.coords$id[g] &
                    x$snp == data$snp.coords$id[p],
                    6]
      betahat <- x[x$gene == data$gene.coords$id[g] &
                   x$snp == data$snp.coords$id[p],
                   7]
      sebetahat <- x[x$gene == data$gene.coords$id[g] &
                     x$snp == data$snp.coords$id[p],
                     8]
      bhat <- betahat / sigmahat
      sebhat <- sebetahat / sigmahat
      t <- qnorm(pt(-abs(bhat/sebhat), N-2-nbCovars, log=TRUE), log=TRUE)
      if(correct){
        if(abs(t) > 10^(-8)){
          sigmahat <- abs(betahat) / (abs(t) * sebhat)
          bhat <- betahat / sigmahat
          sebhat <- bhat / t
        } else{ # abs(t) <= 10^(-8)
          bhat <- 0
          sebhat <- Inf
        }
      }
      c(bhat, sebhat, t)
    }
  }))
  colnames(std.sstats.corr) <- c("bhat", "sebhat", "t")
  return(std.sstats.corr)
}

## Calculate the log10(ABF) of Wen & Stephens (AOAS, 2013))
calcL10Abf <- function(sstats, phi2, oma2){
  l10abf <- NA
  
  bbarhat.num <- 0
  bbarhat.denom <- 0
  varbbarhat <- 0
  l10abfs.single <- c()
  for(i in 1:nrow(sstats)){ # for each subgroup
    if(sum(is.na(sstats[i,])) == length(sstats[i,]))
      next
    bhat <- sstats[i,"bhat"]
    varbhat <- sstats[i,"sebhat"]^2
    t <- sstats[i,"t"]
    if(abs(t) > 10^(-8)){
      bbarhat.num <- bbarhat.num + bhat / (varbhat + phi2)
      bbarhat.denom <- bbarhat.denom + 1 / (varbhat + phi2)
      varbbarhat <- varbbarhat + 1 / (varbhat + phi2)
      tmp <- 0.5 * log10(varbhat) -
        0.5 * log10(varbhat + phi2) +
          (0.5 * t^2 * phi2 / (varbhat + phi2)) / log(10)
    } else
      tmp <- 0
    l10abfs.single <- c(l10abfs.single, tmp)
  }

  if(length(l10abfs.single) != 0){
    bbarhat <- ifelse(bbarhat.denom != 0, bbarhat.num / bbarhat.denom, 0)
    varbbarhat <- ifelse(varbbarhat != 0, 1 / varbbarhat, Inf)
    if(bbarhat != 0 & ! is.infinite(varbbarhat)){
      T2 <- bbarhat^2 / varbbarhat
      l10abf.bbar <- ifelse(T2 != 0,
                            0.5 * log10(varbbarhat) - 0.5 * log10(varbbarhat + oma2) +
                            (0.5 * T2 * oma2 / (varbbarhat + oma2)) / log(10),
                            0)
      l10abf <- l10abf.bbar
      for(i in 1:length(l10abfs.single))
        l10abf <- l10abf + l10abfs.single[i]
    } else
      l10abf <- 0
  }    
  
  return(as.numeric(l10abf))
}

calcL10AbfsRawConstGridL <- function(sstats, gridL){
  abfs <- matrix(data=NA, nrow=3, ncol=nrow(gridL))
  rownames(abfs) <- c("gen", "gen-fix", "gen-maxh")
  for(i in 1:nrow(gridL)){
    abfs[1,i] <- calcL10Abf(sstats, gridL[i,"phi2"], gridL[i,"oma2"])
    abfs[2,i] <- calcL10Abf(sstats, 0, gridL[i,"phi2"] + gridL[i,"oma2"])
    abfs[3,i] <- calcL10Abf(sstats, gridL[i,"phi2"] + gridL[i,"oma2"], 0)
  }
  return(abfs)
}

getConfigNames <- function(configs){
  config.names <- rep(NA, nrow(configs)-1)
  for(i in 2:nrow(configs)) # skip null config
    config.names[i-1] <- paste(which(configs[i,] == 1), collapse="-")
  return(config.names)
}

## std.sstats.corr: matrix with rows s1,s2,s3 and cols bhat,sebhat,t
calcL10AbfsRawAllGridS <- function(std.sstats.corr, gridS, configs){
  l10abfs <- matrix(data=NA, nrow=nrow(configs)-1, ncol=nrow(gridS))
  ## rownames(l10abfs) <- getConfigNames(configs) # -> not the good order...
  rownames(l10abfs) <- c("1", "2", "3", "1-2", "1-3", "2-3", "1-2-3")
  isAbsentInSbgrp <- rep(FALSE, 3)
  for(config in rownames(l10abfs)){
    tmp <- std.sstats.corr[as.numeric(do.call(c, strsplit(config, "-"))),]
    if(length(grep("-", config)) == 0 && sum(is.na(tmp)) == length(tmp))
      isAbsentInSbgrp[as.numeric(config)] <- TRUE
    tmp <- matrix(tmp, nrow=ifelse(is.vector(tmp), 1, nrow(tmp)),
                  ncol=ncol(std.sstats.corr))
    colnames(tmp) <- colnames(std.sstats.corr)
    for(i in 1:nrow(gridS))
      l10abfs[config,i] <- calcL10Abf(tmp, gridS[i,"phi2"], gridS[i,"oma2"])
  }
  
  ## handle genes absent in some subgroups
  if(sum(isAbsentInSbgrp) != 0){
    for(config in rownames(l10abfs)){
      if(length(grep("-", config)) != 0){
        sbgrps <- strsplit(config, "-")[[1]]
        absent.sbgrps <- sapply(sbgrps, function(s){isAbsentInSbgrp[as.numeric(s)]})
        if(sum(absent.sbgrps) != 0){
          config.present <- paste(sbgrps[! absent.sbgrps], collapse="-")
          l10abfs[config,] <- l10abfs[config.present,]
        }
      } else if(isAbsentInSbgrp[as.numeric(config)]) # absent singleton
        l10abfs[config,] <- rep(0, ncol(gridS))
    }
  }
  
  return(l10abfs)
}

## Return the log10(ABF) from Wen (Biometrics, accepted)
## using raw data or summary statistics
##
## betag.hat is a vector [SNP 1 in all tissues, SNP 2 in all tissues, etc]
## Vg is a matrix (covariance of betag.hat)
calcL10AbfMvlr <- function(Y=NULL, Xg=NULL, Xc=NULL,
                           betag.hat=NULL, Sigma.hat=NULL, Vg=NULL,
                           gamma, phi2, oma2, alpha, model="ES", debug=FALSE){
  require(MASS)
  if((is.null(Y) || is.null(Xg) || is.null(Xc)) &&
     (is.null(betag.hat) || is.null(Sigma.hat) || is.null(Vg)))
    stop("need to provide Y, Xg and Xc, or betag.hat, Sigma.hat and Vg")
  if(alpha != 0.0)
    warning("the small sample correction is not performed")
  if(debug){message(paste0("phi2=",phi2," oma2=",oma2))}
  
  S <- ifelse(!is.null(Y), ncol(Y), ncol(Vg)) # nb of subgroups
  W <- matrix(rep(oma2,S*S),S,S) + diag(rep(phi2,S),S,S)
  ## if(debug){message("W="); print.mat2(W, 4, 4)}
  P <- ifelse(!is.null(Xg), ncol(Xg), length(betag.hat)/S) # nb of snps
  Wg <- diag(rep(1,P)) %x% W
  Wg <- Wg * (gamma %*% t(gamma))
  ## if(debug){message("W="); print.mat2(W, 4, 4)}
  
  # inputs are raw data
  if(! is.null(Y) && ! is.null(Xg) && ! is.null(Xc)){
    ## add intercept
    if(ncol(Xc) == 0){
      Xc <- matrix(1, nrow=nrow(Xg), ncol=1)
    } else
      Xc <- matrix(cbind(1,Xc), nrow=nrow(Xg), ncol=ncol(Xc)+1)
    X <- cbind(Xc, Xg)
    
    N <- dim(Y)[1] # nb of individuals
    Q <- dim(Xc)[2] # nb of other covariates + intercept
    m <- Q+S+1 # degrees of freedom for the Wishart
    ## if(debug){message(paste0("N=",N," P=",P," S=",S," Q=",Q," m=",m))}
    
    Sigma.hat.full <- (t(Y)%*%(diag(rep(1,N)) - X%*%ginv(t(X)%*%X)%*%t(X))%*%Y + diag(rep(1e-8,S),S,S)) / (N+m-Q-S-1)
    Sigma.hat.null <- (t(Y)%*%(diag(rep(1,N)) - Xc%*%solve(t(Xc)%*%Xc)%*%t(Xc))%*%Y + diag(rep(1e-8,S),S,S)) / (N+m-Q-S-1)
    ## if(debug){message("Sigma.hat.full="); print.mat2(Sigma.hat.full)}
    ## if(debug){message("Sigma.hat.null="); print.mat2(Sigma.hat.null)}
    
    Sigma.hat <- alpha*Sigma.hat.full + (1-alpha)*Sigma.hat.null
    Sigma.hat.inv <- solve(Sigma.hat)
    ## if(debug){message("Sigma.hat="); print.mat2(Sigma.hat)}
    
    Vg.inv <- (t(Xg)%*%Xg - t(Xg)%*%Xc%*%solve(t(Xc)%*%Xc)%*%t(Xc)%*%Xg) %x% Sigma.hat.inv
    
    vec <- matrix(as.vector(t(Y-Xc%*%solve(t(Xc)%*%Xc)%*%t(Xc)%*%Y)), ncol=1)
    
    bVi <- t(vec) %*% (Xg %x% Sigma.hat.inv)
  } else{ # inputs are summary stats
    Vg.inv <- solve(Vg)
    ## if(debug){message("Vg.inv="); print.mat2(Vg.inv)}
    bVi <- t(betag.hat) %*% Vg.inv
    ## if(debug){message("bVi="); print.mat2(bVi)}
  }
  
  if(model == "ES"){
    if(nrow(Sigma.hat) == 1 && ncol(Sigma.hat)){
      Wg <- matrix(Sigma.hat,1,1)^.5 %*% Wg %*% matrix(Sigma.hat,1,1)^.5
    } else
      Wg <- diag(diag(Sigma.hat))^.5 %*% Wg %*% diag(diag(Sigma.hat))^.5
  }
  ## if(debug){message("Wg="); print.mat2(Wg)}
  
  ivw <- diag(rep(1,P*S)) + Vg.inv %*% Wg
  ## if(debug){message("ivw="); print.mat2(ivw)}
  
  ## if(debug){message(paste0("det=",determinant(ivw)$modulus[[1]]))}
  
  l10abf <- (- 0.5 * determinant(ivw)$modulus[[1]] +
             0.5 * bVi %*% Wg %*% solve(ivw) %*% t(bVi)) /
               log(10)
  ## if(debug){message(paste0("l10abf=",l10abf))}
  
  return(l10abf)
}

calcL10AbfsRawConstGridLMvlr <- function(data, g, p, gridL, prop.fit.sigma){
  l10abfs <- matrix(data=NA, nrow=3, ncol=nrow(gridL))
  rownames(l10abfs) <- c("gen", "gen-fix", "gen-maxh")
  
  Y <- matrix(do.call(cbind, lapply(data$phenos, function(X){X[g,]})),
              nrow=nrow(data$inds), ncol=length(data$phenos))
  Xg <- matrix(data$geno.counts[p,], nrow=nrow(data$inds), ncol=1)
  Xc <- matrix(nrow=0, ncol=0)
  gamma <- matrix(rep(1, ncol(Y)), ncol=1)
  
  for(i in 1:nrow(gridL)){
    l10abfs[1,i] <- calcL10AbfMvlr(Y=Y, Xg=Xg, Xc=Xc, gamma=gamma,
                                   phi2=gridL[i,"phi2"], oma2=gridL[i,"oma2"],
                                   alpha=prop.fit.sigma)
    l10abfs[2,i] <- calcL10AbfMvlr(Y=Y, Xg=Xg, Xc=Xc, gamma=gamma,
                                   phi2=0, oma2=gridL[i,"phi2"] + gridL[i,"oma2"],
                                   alpha=prop.fit.sigma)
    l10abfs[3,i] <- calcL10AbfMvlr(Y=Y, Xg=Xg, Xc=Xc, gamma=gamma,
                                   phi2=gridL[i,"phi2"] + gridL[i,"oma2"], oma2=0,
                                   alpha=prop.fit.sigma)
  }
  
  return(l10abfs)
}

calcL10AbfsRawAllGridSMvlr <- function(data, g, p, gridS, configs,
                                       prop.fit.sigma){
  l10abfs <- matrix(data=NA, nrow=nrow(configs)-1, ncol=nrow(gridS))
  rownames(l10abfs) <- c("1", "2", "3", "1-2", "1-3", "2-3", "1-2-3")
  
  Y <- matrix(do.call(cbind, lapply(data$phenos, function(X){X[g,]})),
              nrow=nrow(data$inds), ncol=length(data$phenos))
  Xg <- matrix(data$geno.counts[p,], nrow=nrow(data$inds), ncol=1)
  Xc <- matrix(nrow=0, ncol=0)
  
  for(config in rownames(l10abfs)){
    gamma <- rep(0, ncol(Y))
    gamma[as.numeric(strsplit(config, "-")[[1]])] <- 1
    gamma <- matrix(gamma, ncol=1)
    for(i in 1:nrow(gridS))
      l10abfs[config,i] <- calcL10AbfMvlr(Y=Y, Xg=Xg, Xc=Xc, gamma=gamma,
                                          phi2=gridS[i,"phi2"],
                                          oma2=gridS[i,"oma2"],
                                          alpha=prop.fit.sigma)
  }
  
  return(l10abfs)
}

calcBetahatsAndDiagsPerSubgroup <- function(y, Xg, Xc, alpha=0.0, verbose=0){
  stopifnot(is.list(y), is.list(Xg), is.list(Xc))
  
  S <- length(y) # nb of subgroups
  
  tmp <- lapply(1:S, function(s){
    Y <- y[[s]] # N x 1
    N <- nrow(Y)
    
    ## add intercept
    if(ncol(Xc[[s]]) == 0){
      Xc[[s]] <- matrix(1, nrow=N, ncol=1)
    } else
      Xc[[s]] <- matrix(cbind(1,Xc[[s]]), nrow=N, ncol=ncol(Xc[[s]])+1)
    Q <- ncol(Xc[[s]])
    X <- cbind(Xc[[s]], Xg[[s]]) # N x (Q+P) where P=1 (single SNP)
    idx.geno <- ncol(X) # index of column corresponding to genotype
    
    ## if(verbose > 0){message("MLE of full model")}
    X.svd <- svd(X) # X = U D V'
    if(length(X.svd$d) > 1){
      D.inv <- diag(1/X.svd$d)
    } else
      D.inv <- diag(1/X.svd$d,1,1)
    X.ps <- X.svd$v %*% D.inv %*% t(X.svd$u) # pseudo-inverse
    B.hat.full <- X.ps %*% Y
    E.hat <- Y - X %*% B.hat.full # residuals
    sigma2.hat.full <- ((1/N) * t(E.hat) %*% E.hat)[1,1]
    ## if(verbose > 0){message(paste0("sigma2.hat.full=",sigma2.hat.full))}
    
    ## if(verbose > 0){message("MLE of null model")}
    Xc.svd <- svd(Xc[[s]])
    if(length(Xc.svd$d) > 1){
      Dc.inv <- diag(1/Xc.svd$d)
    } else
      Dc.inv <- diag(1/Xc.svd$d,1,1)
    Xc.ps <- Xc.svd$v %*% Dc.inv %*% t(Xc.svd$u)
    B.hat.null <- Xc.ps %*% Y
    E.hat <- Y - Xc[[s]] %*% B.hat.null
    sigma2.hat.null <- ((1/N) * t(E.hat) %*% E.hat)[1,1]
    ## if(verbose > 0){message(paste0("sigma2.hat.null=",sigma2.hat.null))}
    
    ## final estimates
    betag.hat <- B.hat.full[idx.geno,1]
    sigma2.hat <- alpha*sigma2.hat.full + (1-alpha)*sigma2.hat.null
    varbetag.hat <- (sigma2.hat * X.svd$v %*% D.inv^2 %*% t(X.svd$v))[idx.geno,idx.geno]
    
    return(list(betag.hat=betag.hat,
                sigma2.hat=sigma2.hat,
                varbetag.hat=varbetag.hat))
  })
  
  return(list(betag.hat=do.call(c, lapply(tmp, function(x){x$betag.hat})),
              sigma2.hat=do.call(c, lapply(tmp, function(x){x$sigma2.hat})),
              varbetag.hat=do.call(c, lapply(tmp, function(x){x$varbetag.hat}))))
}

## in subgroup s: A = (Xc'Xc + Xu'Xu)^(-1) %*% Xc'
getMatricesA <- function(inds, X, s1, s2, verbose=0){
  A <- list()
  
  inds.com.s1s2 <- inds[[s1]][which(inds[[s1]] %in% inds[[s2]])]
  Xs1s2 <- X[[s1]][inds.com.s1s2,]
  tXs1s2.Xs1s2 <- t(Xs1s2) %*% Xs1s2
  
  inds.uniq.s1 <- inds[[s1]][which(! inds[[s1]] %in% inds[[s2]])]
  if(length(inds.uniq.s1) == 0){
    A[["s1"]] <- solve(tXs1s2.Xs1s2) %*% t(Xs1s2)
  } else{
    Xu.s1 <- X[[s1]][inds.uniq.s1,]
    A[["s1"]] <- solve(tXs1s2.Xs1s2 + t(Xu.s1) %*% Xu.s1) %*% t(Xs1s2)
  }
  
  inds.uniq.s2 <- inds[[s2]][which(! inds[[s2]] %in% inds[[s1]])]
  if(length(inds.uniq.s2) == 0){
    A[["s2"]] <- solve(tXs1s2.Xs1s2) %*% t(Xs1s2)
  } else{
    Xu.s2 <- X[[s2]][inds.uniq.s2,]
    A[["s2"]] <- solve(tXs1s2.Xs1s2 + t(Xu.s2) %*% Xu.s2) %*% t(Xs1s2)
  }
  
  return(A)
}

## warning, no small sample correction
getErrCovSigmaBtwPairSubgroups <- function(inds, X, y, s1, s2){
  Sigma <- list()
  
  inds.com.s1s2 <- inds[[s1]][which(inds[[s1]] %in% inds[[s2]])]
  N.s1s2 <- length(inds.com.s1s2)
  Y.s1s2 <- cbind(y[[s1]][inds.com.s1s2,], y[[s2]][inds.com.s1s2,])
  X.s1s2 <- X[[s1]][inds.com.s1s2,]
  Xc.s1s2 <- X.s1s2[,-2]
  
  X.s1s2.svd <- svd(X.s1s2)
  Bhat.full <- X.s1s2.svd$v %*% diag(1/X.s1s2.svd$d) %*% t(X.s1s2.svd$u) %*% Y.s1s2
  Sigma[["full"]] <- (1/N.s1s2) * t(Y.s1s2 - X.s1s2 %*% Bhat.full) %*% (Y.s1s2 - X.s1s2 %*% Bhat.full)
  
  Xc.s1s2.svd <- svd(Xc.s1s2)
  if(length(Xc.s1s2.svd$d) > 1){
    Dc.s1s2.inv <- diag(1/Xc.s1s2.svd$d)
  } else
    Dc.s1s2.inv <- diag(1/Xc.s1s2.svd$d,1,1)
  Bhat.null <- Xc.s1s2.svd$v %*% Dc.s1s2.inv %*% t(Xc.s1s2.svd$u) %*% Y.s1s2
  Sigma[["null"]] <- (1/N.s1s2) * t(Y.s1s2 - Xc.s1s2 %*% Bhat.null) %*% (Y.s1s2 - Xc.s1s2 %*% Bhat.null)
  
  return(Sigma)
}

calcOffDiagCovarsFromPairsOfSubgroups <- function(y, Xg, Xc, sigma2.hat,
                                                  varbetag.hat, alpha=0.0,
                                                  verbose=0){
  stopifnot(is.list(y), is.list(Xg), is.list(Xc), is.vector(sigma2.hat),
            is.vector(varbetag.hat))
  
  S <- length(y) # nb of subgroups
  Sigma.hat <- matrix(0, nrow=S, ncol=S)
  diag(Sigma.hat) <- sigma2.hat
  Vg <- matrix(0, nrow=S, ncol=S)
  diag(Vg) <- varbetag.hat
  
  inds <- lapply(y, rownames)
  
  X <- lapply(1:S, function(s){
    if(ncol(Xc[[s]]) == 0){
      tmp <- matrix(cbind(1, Xg[[s]]), nrow=nrow(Xg[[s]]), ncol=2)
    } else
      tmp <- matrix(cbind(1,Xg[[s]],Xc[[s]]), nrow=nrow(Xg[[s]]), ncol=ncol(Xc[[s]])+2)
    rownames(tmp) <- rownames(Xg[[s]])
    tmp
  })
  
  ## for each pair of subgroups, get the covariance between residuals
  ## and between genotype betahats
  for(s1 in 1:(S-1)){
    for(s2 in (s1+1):S){
      A <- getMatricesA(inds, X, s1, s2)
      
      Sigmas.s1s2.hat <- getErrCovSigmaBtwPairSubgroups(inds, X, y, s1, s2)
      
      Sigma.s1s2.hat <- alpha * Sigmas.s1s2.hat[["full"]] +
        (1-alpha) * Sigmas.s1s2.hat[["null"]]
      Sigma.hat[s1,s2] <- Sigma.s1s2.hat[1,2]
      
      cov.betahat1.betahat2 <- Sigma.hat[s1,s2] * A[["s1"]] %*% t(A[["s2"]])
      Vg[s1,s2] <- cov.betahat1.betahat2[2,2]
    }
  }
  
  ## fill the lower part of each symmetric matrix
  Sigma.hat <- Sigma.hat + t(Sigma.hat) - diag(diag(Sigma.hat))
  Vg <- Vg + t(Vg) - diag(diag(Vg))
  
  return(list(Sigma.hat=Sigma.hat, Vg=Vg))
}

calcSstatsBoth <- function(data, g, p, prop.fit.sigma){
  S <- length(data$phenos)
  
  ## prepare raw data for each subgroup
  y <- lapply(1:S, function(s){
    tmp <- matrix(data$phenos[[s]][g,], ncol=1)
    rownames(tmp) <- colnames(data$phenos[[s]])
    tmp
  })
  Xg <- lapply(1:S, function(s){
    tmp <- matrix(data$geno.counts[p,colnames(data$phenos[[s]])], ncol=1)
    rownames(tmp) <- colnames(data$phenos[[s]])
    tmp
  })
  Xc <- lapply(1:S, function(s){
    matrix(nrow=0, ncol=0)
  })
  
  ## for each subgroup separately, calc betag.hat, diag of Vg
  ## and diags of Sigma.hat.full and Sigma.hat.null
  sstats.uvlr <- calcBetahatsAndDiagsPerSubgroup(
    y=y, Xg=Xg, Xc=Xc, alpha=prop.fit.sigma, verbose=0)
  
  ## for common individuals between each pair of subgroups,
  ## calc each off-diagonal element of Vg as well as both Sigma's
  sstats.mvlr <- calcOffDiagCovarsFromPairsOfSubgroups(
    y=y, Xg=Xg, Xc=Xc, sigma2.hat=sstats.uvlr$sigma2.hat,
    varbetag.hat=sstats.uvlr$varbetag.hat, alpha=prop.fit.sigma, verbose=0)
  
  return(list(betag.hat=sstats.uvlr$betag.hat,
              Sigma.hat=sstats.mvlr$Sigma.hat,
              Vg=sstats.mvlr$Vg))
}

calcL10AbfsRawConstGridLFromSstats <- function(sstats.both, gridL,
                                               prop.fit.sigma,
                                               debug=FALSE){
  l10abfs <- matrix(data=NA, nrow=3, ncol=nrow(gridL))
  rownames(l10abfs) <- c("gen", "gen-fix", "gen-maxh")
  
  S <- length(sstats.both$betag.hat)
  gamma <- matrix(rep(1, S), ncol=1)
  
  for(i in 1:nrow(gridL)){
    l10abfs[1,i] <- calcL10AbfMvlr(betag.hat=sstats.both$betag.hat,
                                   Sigma.hat=sstats.both$Sigma.hat,
                                   Vg=sstats.both$Vg,
                                   gamma=gamma,
                                   phi2=gridL[i,"phi2"],
                                   oma2=gridL[i,"oma2"],
                                   alpha=prop.fit.sigma)
    l10abfs[2,i] <- calcL10AbfMvlr(betag.hat=sstats.both$betag.hat,
                                   Sigma.hat=sstats.both$Sigma.hat,
                                   Vg=sstats.both$Vg,
                                   gamma=gamma,
                                   phi2=0,
                                   oma2=gridL[i,"phi2"] + gridL[i,"oma2"],
                                   alpha=prop.fit.sigma)
    l10abfs[3,i] <- calcL10AbfMvlr(betag.hat=sstats.both$betag.hat,
                                   Sigma.hat=sstats.both$Sigma.hat,
                                   Vg=sstats.both$Vg,
                                   gamma=gamma,
                                   phi2=gridL[i,"phi2"] + gridL[i,"oma2"],
                                   oma2=0,
                                   alpha=prop.fit.sigma)
  }

  if(debug)
    tmp <- calcL10AbfMvlr(betag.hat=sstats.both$betag.hat,
                          Sigma.hat=sstats.both$Sigma.hat,
                          Vg=sstats.both$Vg,
                          gamma=gamma,
                          phi2=0.0075,
                          oma2=0.0025,
                          alpha=prop.fit.sigma,
                          debug=TRUE)
  
  return(l10abfs)
}

calcL10AbfsRawAllGridSFromSstats <- function(sstats.both, gridS, configs,
                                             prop.fit.sigma){
  l10abfs <- matrix(data=NA, nrow=nrow(configs)-1, ncol=nrow(gridS))
  rownames(l10abfs) <- c("1", "2", "3", "1-2", "1-3", "2-3", "1-2-3")
  
  S <- length(sstats.both$betag.hat)
  
  for(config in rownames(l10abfs)){
    gamma <- rep(0, S)
    gamma[as.numeric(strsplit(config, "-")[[1]])] <- 1
    gamma <- matrix(gamma, ncol=1)
    for(i in 1:nrow(gridS))
      l10abfs[config,i] <- calcL10AbfMvlr(betag.hat=sstats.both$betag.hat,
                                          Sigma.hat=sstats.both$Sigma.hat,
                                          Vg=sstats.both$Vg,
                                          gamma=gamma,
                                          phi2=gridS[i,"phi2"],
                                          oma2=gridS[i,"oma2"],
                                          alpha=prop.fit.sigma)
  }
  
  return(l10abfs)
}

calcRawAbfsOnSimulatedData <- function(data=NULL, sstats=NULL, grids=NULL,
                                       withCovars=FALSE, mvlr=FALSE,
                                       rmvIndsFromSbgrp=FALSE,
                                       lik="norm", verbose=0){
  stopifnot(! is.null(data), ! is.null(sstats), ! is.null(grids))
  if(verbose > 0)
    message("calculate raw ABFs ...")
  l10abfs.raw <- data.frame(gene=rep(NA, 10*data$params$nb.pairs),
                            snp=NA, config=NA)
  for(i in 1:nrow(grids$gridL)){
    l10abfs.raw <- cbind(l10abfs.raw, NA)
    names(l10abfs.raw)[3+i] <- paste0("l10abf.grid",i)
  }
  i <- 1
  for(g in 1:nrow(data$gene.coords)){
    for(p in 1:nrow(data$snp.coords)){
      if(data$snp.coords$chr[p] != data$gene.coords$chr[g] ||
         data$snp.coords$start[p] < data$gene.coords$start[g] - data$params$len.cis ||
         data$snp.coords$start[p] > data$gene.coords$start[g] + data$params$len.cis)
        next # SNP not in cis
      isAbsentInSbgrp <- do.call(c, lapply(data$phenos, function(X){
        ! data$gene.coords$id[g] %in% rownames(X)}))
      if(mvlr && any(isAbsentInSbgrp))
        next # gene unexpressed in some subgroups
      nbCovars <- 0
      if(withCovars && "sex" %in% names(data$inds))
        nbCovars <-  1
      if(! mvlr && ! rmvIndsFromSbgrp){
        std.sstats.corr <- getStdSstatsAndCorrSmallSampleSize(data, sstats,
                                                              g, p,
                                                              nbCovars)
        l10abfs.raw.const.gridL <- calcL10AbfsRawConstGridL(std.sstats.corr,
                                                            grids$gridL)
        l10abfs.raw.all.gridS <- calcL10AbfsRawAllGridS(std.sstats.corr,
                                                        grids$gridS,
                                                        data$configs)
      } else if(mvlr && ! rmvIndsFromSbgrp){
        l10abfs.raw.const.gridL <- calcL10AbfsRawConstGridLMvlr(data, g, p,
                                                                grids$gridL,
                                                                0.0)
        l10abfs.raw.all.gridS <- calcL10AbfsRawAllGridSMvlr(data, g, p,
                                                            grids$gridS,
                                                            data$configs, 0.0)
      } else{
        sstats.both <- calcSstatsBoth(data, g, p, 0.0)
        l10abfs.raw.const.gridL <- calcL10AbfsRawConstGridLFromSstats(
          sstats.both, grids$gridL, 0.0, debug=ifelse(data$gene.coords$id[g] == "gene1" && data$snp.coords$id[p] == "snp1", TRUE, FALSE))
        l10abfs.raw.all.gridS <- calcL10AbfsRawAllGridSFromSstats(
          sstats.both, grids$gridS, data$configs, 0.0)
      }
      for(m in 1:nrow(l10abfs.raw.const.gridL)){
        l10abfs.raw$gene[i] <- data$gene.coords$id[g]
        l10abfs.raw$snp[i] <- data$snp.coords$id[p]
        l10abfs.raw$config[i] <- rownames(l10abfs.raw.const.gridL)[m]
        l10abfs.raw[i,-c(1:3)] <- l10abfs.raw.const.gridL[m,]
        i <- i + 1
      }
      for(m in 1:nrow(l10abfs.raw.all.gridS)){
        l10abfs.raw$gene[i] <- data$gene.coords$id[g]
        l10abfs.raw$snp[i] <- data$snp.coords$id[p]
        l10abfs.raw$config[i] <- rownames(l10abfs.raw.all.gridS)[m]
        l10abfs.raw[i,4:(4+ncol(l10abfs.raw.all.gridS)-1)] <- l10abfs.raw.all.gridS[m,]
        i <- i + 1
      }
    }
  }
  return(l10abfs.raw)
}

log10WeightedSum <- function(values=NULL, weights=NULL){
  stopifnot(! is.null(values))
  if(is.null(weights))
    weights <- rep(1/length(values), length(values))
  max <- max(values)
  return(max + log10(sum(weights * 10^(values - max))))
}

calcAvgAbfsSinAndGenSin <- function(l10abfs.avg, i){
  tmp <- c()
  weights <- c()
  
  for(config in c("l10abf.1", "l10abf.2", "l10abf.3")){
    if(! is.na(l10abfs.avg[i, config])){
      tmp <- c(tmp, l10abfs.avg[i, config])
      weights <- c(weights, (1/2)*(1/choose(3,1)))
    }
  }
  
  tmp <- c(tmp, l10abfs.avg$l10abf.gen[i])
  weights <- c(weights, (1/2))
  l10abfs.avg$l10abf.gen.sin[i] <- log10WeightedSum(tmp, weights)
  
  return(l10abfs.avg)
}

calcAvgAbfAll <- function(l10abfs.avg, i, configs){
  tmp <- c()
  weights <- c()
  for(config in paste0("l10abf.", configs)){
    k <- length(strsplit(config, "-")[[1]])
    if(! is.na(l10abfs.avg[i, config])){
      tmp <- c(tmp, l10abfs.avg[i, config])
      weights <- c(weights, (1/3) * (1/choose(3,k)))
    }
  }
  l10abfs.avg$l10abf.all[i] <- log10WeightedSum(tmp, weights)
  return(l10abfs.avg)
}

calcAvgAbfsOnSimulatedData <- function(data=NULL, l10abfs.raw=NULL, mvlr=FALSE,
                                       rmvIndsFromSbgrp=FALSE, verbose=0){
  stopifnot(! is.null(data), ! is.null(l10abfs.raw))
  if(verbose > 0)
    message("calculate average ABFs ...")
  l10abfs.avg <- data.frame(gene=rep(NA, data$params$nb.pairs), snp=NA,
                            nb.subgroups=NA)
  for(i in c("l10abf.gen", "l10abf.gen.fix", "l10abf.gen.maxh",
             "l10abf.gen.sin", "l10abf.all")){
    l10abfs.avg <- cbind(l10abfs.avg, NA)
    colnames(l10abfs.avg)[ncol(l10abfs.avg)] <- i
  }
  configs <- unique(l10abfs.raw$config)[-c(1:3)]
  for(i in paste0("l10abf.",configs)){
    l10abfs.avg <- cbind(l10abfs.avg, NA)
    colnames(l10abfs.avg)[ncol(l10abfs.avg)] <- i
  }
  
  i <- 1
  for(g in 1:nrow(data$gene.coords)){
    for(p in 1:nrow(data$snp.coords)){
      if(data$snp.coords$chr[p] != data$gene.coords$chr[g] ||
         data$snp.coords$start[p] < data$gene.coords$start[g] - data$params$len.cis ||
         data$snp.coords$start[p] > data$gene.coords$start[g] + data$params$len.cis)
        next # SNP not in cis
      isAbsentInSbgrp <- do.call(c, lapply(data$phenos, function(X){
        ! data$gene.coords$id[g] %in% rownames(X)}))
      if(mvlr && any(isAbsentInSbgrp))
        next # gene unexpressed in some subgroups
      if(verbose > 0)
        message(paste(data$gene.coords$id[g], data$snp.coords$id[p]))
      
      l10abfs.avg$gene[i] <- data$gene.coords$id[g]
      l10abfs.avg$snp[i] <- data$snp.coords$id[p]
      l10abfs.avg$nb.subgroups[i] <- 0
      
      tmp <- l10abfs.raw[l10abfs.raw$gene == l10abfs.avg$gene[i] &
                         l10abfs.raw$snp == l10abfs.avg$snp[i] &
                         l10abfs.raw$config == "gen",
                         c(4:ncol(l10abfs.raw))]
      l10abfs.avg$l10abf.gen[i] <- log10WeightedSum(tmp)
      tmp <- l10abfs.raw[l10abfs.raw$gene == l10abfs.avg$gene[i] &
                         l10abfs.raw$snp == l10abfs.avg$snp[i] &
                         l10abfs.raw$config == "gen-fix",
                         c(4:ncol(l10abfs.raw))]
      l10abfs.avg$l10abf.gen.fix[i] <- log10WeightedSum(tmp)
      tmp <- l10abfs.raw[l10abfs.raw$gene == l10abfs.avg$gene[i] &
                         l10abfs.raw$snp == l10abfs.avg$snp[i] &
                         l10abfs.raw$config == "gen-maxh",
                         c(4:ncol(l10abfs.raw))]
      l10abfs.avg$l10abf.gen.maxh[i] <- log10WeightedSum(tmp)
      
      for(config in configs){
        tmp <- l10abfs.raw[l10abfs.raw$gene == l10abfs.avg$gene[i] &
                           l10abfs.raw$snp == l10abfs.avg$snp[i] &
                           l10abfs.raw$config == config,
                           c(4:13)]
        if(length(grep("-", config)) == 0 && sum(tmp == 0.0) != length(tmp))
          l10abfs.avg$nb.subgroups[i] <- l10abfs.avg$nb.subgroups[i] + 1
        l10abfs.avg[i, paste0("l10abf.",config)] <- log10WeightedSum(tmp)
      }
      
      l10abfs.avg <- calcAvgAbfsSinAndGenSin(l10abfs.avg, i)
      
      l10abfs.avg <- calcAvgAbfAll(l10abfs.avg, i, configs)
      
      i <- i + 1
    }
  }
  return(l10abfs.avg)
}

getResultsOnSimulatedData <- function(data=NULL, grids=NULL, withCovars=FALSE,
                                      mvlr=FALSE, rmvIndsFromSbgrp=FALSE,
                                      rmvExpFromSbgrpsInds=FALSE, lik="norm",
                                      verbose=0){
  stopifnot(! is.null(data), ! is.null(grids))
  if(verbose > 0)
    message("get results on simulated data with R ...")
  
  sstats <- calcSstatsOnSimulatedData(data, withCovars, lik, verbose-1)
  
  if(length(grep("pois", lik)) == 0){
    l10abfs.raw <- calcRawAbfsOnSimulatedData(data, sstats, grids, withCovars,
                                              mvlr, rmvIndsFromSbgrp, lik,
                                              verbose-1)
    l10abfs.avg <- calcAvgAbfsOnSimulatedData(data, l10abfs.raw, mvlr,
                                              rmvIndsFromSbgrp, verbose-1)
  } else{
    l10abfs.raw <- NA
    l10abfs.avg <- NA
  }
  
  return(list(sstats=sstats,
              l10abfs.raw=l10abfs.raw,
              l10abfs.avg=l10abfs.avg))
}

writeResultsOnSimulatedData <- function(res=NULL, verbose=0){
  stopifnot(! is.null(res))
  
  if(verbose > 0)
    message("write results (expected) ...")
  
  for(s in names(res$sstats)){
    row.ids <- which(res$sstats[[s]]$n != 0) # don't save absent genes
    tmp <- cbind(res$sstats[[s]][row.ids,1:2],
                 sprintf(fmt="%.06e", res$sstats[[s]][row.ids,3]),
                 res$sstats[[s]][row.ids,4])
    for(j in 5:9) #ncol(res$sstats[[s]])) <- don't save covariates results
      tmp <- cbind(tmp, sprintf(fmt="%.06e", res$sstats[[s]][row.ids,j]))
    colnames(tmp) <- colnames(res$sstats[[s]])[1:9]
    if(tmp[1,"pve"] == "NA") # this means "lik=norm" is false
      tmp[,"pve"] <- NA
    write.table(x=tmp, file=gzfile(paste0("exp_bf_sumstats_",s,".txt.gz")),
                quote=FALSE, row.names=FALSE, sep="\t", na="nan")
  }

  if(is.data.frame(res$l10abfs.raw)){
    tmp <- res$l10abfs.raw[,1:3]
    for(j in 4:ncol(res$l10abfs.raw))
      tmp <- cbind(tmp, gsub("NA", "nan", sprintf(fmt="%.06e", res$l10abfs.raw[,j])))
    colnames(tmp) <- colnames(res$l10abfs.raw)
    write.table(x=tmp, file=gzfile("exp_bf_l10abfs_raw.txt.gz"),
                quote=FALSE, row.names=FALSE, sep="\t", na="nan")
    
    tmp <- res$l10abfs.avg[,1:3]
    for(j in 4:ncol(res$l10abfs.avg))
      tmp <- cbind(tmp, gsub("NA", "nan", sprintf(fmt="%.06e", res$l10abfs.avg[,j])))
    colnames(tmp) <- colnames(res$l10abfs.avg)
    write.table(x=tmp, file=gzfile("exp_bf_l10abfs_avg-grids.txt.gz"),
                quote=FALSE, row.names=FALSE, sep="\t", na="nan")
  }
}

run <- function(params){
  setwd(params$dir.name)
  data <- getSimulatedData(params$rmvGenesFromSbgrp, params$rmvIndsFromSbgrp,
                           params$rmvExpFromSbgrpsInds, params$lik,
                           params$verbose)
  writeSimulatedData(data, geno.format=params$geno.format,
                     verbose=params$verbose)
  grids <- getGrids()
  writeGrids(grids, params$verbose)
  res <- getResultsOnSimulatedData(data, grids, params$withCovars,
                                   params$mvlr, params$rmvIndsFromSbgrp,
                                   params$rmvExpFromSbgrpsInds, params$lik,
                                   params$verbose)
  writeResultsOnSimulatedData(res, params$verbose)
}

main <- function(){
  params <- list(verbose=1,
                 dir.name=NULL,
                 rmvGenesFromSbgrp=FALSE,
                 withCovars=FALSE,
                 mvlr=FALSE,
                 rmvIndsFromSbgrp=FALSE,
                 rmvExpFromSbgrpsInds=FALSE,
                 lik="norm",
                 geno.format="custom")
  params <- parseCmdLine(params)
  checkParams(params)
  if(params$verbose > 0)
    message(paste0("START ", prog.name, " (", date(), ")"))
  
  run(params)
  
  if(params$verbose > 0)
    message(paste0("END ", prog.name, " (", date(), ")"))
}

main()
