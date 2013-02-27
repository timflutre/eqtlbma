## `functional_tests.R' simulate a typical eQTL data set in multiple
## tissues and analyze it with R in order to test `eqtlbma'.
## Copyright (C) 2012-2013 Timothee Flutre

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

parseArgs <- function(){
  params <- list(verbose=1,
                 dir.name=NULL,
                 rmvGenesFromSbgrps=FALSE,
                 withCovars=FALSE)
  
  args <- commandArgs(trailingOnly=TRUE)
  ## print(args)
  stopifnot(length(args) != 0)
  
  for(i in seq(1, length(args), 2)){
    if(args[i] == "-v" || args[i] == "--verbose")
      params$verbose <- as.numeric(args[i+1])
    else if(args[i] == "--dir")
      params$dir.name <- args[i+1]
    else if(args[i] == "--rgs")
      params$rmvGenesFromSbgrps <- TRUE
    else if(args[i] == "--cvrt")
      params$withCovars <- TRUE
  }
  if(params$verbose > 1)
    print(params)
  
  return(params)
}

##-----------------------------------------------------------------------------
## list of functions

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

getSimulatedData <- function(rmvGenesFromSbgrps=FALSE,
                             withCovars=FALSE,
                             verbose=0){
  if(verbose > 0)
    message("simulate data with R ...")

  params <- list()
  params$seed <- 1859
  set.seed(params$seed)
  
  nb.inds <- 500
  inds <- data.frame(id=paste0("ind", 1:nb.inds),
                     name=paste0("individual ", 1:nb.inds),
                     stringsAsFactors=FALSE)
  if(withCovars)
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
  pves.g <- runif(n=nb.pairs, min=0.4, max=0.8) # explained by genotype
  hets <- runif(n=nb.pairs, min=0, max=0.4)# heterogeneity
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
  
  for(i in 1:nb.pairs){
    for(s in 1:nb.subgroups){
      phenos[[s]][truth$gene[i],] <- truth[i,paste0("mu",s)] +
        truth[i,paste0("b",s)] * truth[i,paste0("sigma",s)] *
          geno.counts[truth[i,"snp"],] +
            rnorm(n=nb.inds, mean=0, sd=truth[i,paste0("sigma",s)])
    }
  }
  
  ## rmv gene9 in s3
  if(rmvGenesFromSbgrps){
    phenos[["s3"]] <- phenos[["s3"]][-9,]
    tmp <- truth
    tmp[tmp$gene == "gene9", "config"] <-
      paste(do.call(c, strsplit(tmp[tmp$gene == "gene9", "config"], ""))[-3],
            collapse="")
    tmp[tmp$gene == "gene9", c("mu3", "b3", "sigma3")] <- NA
    truth <- tmp
  }
  
  return(list(inds=inds,
              gene.coords=gene.coords,
              snp.coords=snp.coords,
              geno.counts=geno.counts,
              phenos=phenos,
              truth=truth,
              params=params,
              configs=configs))
}

writeSimulatedData <- function(data=NULL, geno.format="custom", verbose=0){
  stopifnot(! is.null(data))
  
  if(verbose > 0)
    message("write simulated data ...")
  
  write.table(x=data$gene.coords, file=gzfile("gene_coords.bed.gz"),
              quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
  if(geno.format == "custom"){
    write.table(x=data$snp.coords, file=gzfile("snp_coords.bed.gz"),
                quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
    write.table(x=data$geno.counts, file=gzfile("genotypes.txt.gz"), quote=FALSE)
    tmp <- data.frame(subgroup=names(data$phenos), file=rep("genotypes.txt.gz",
                                                     length(data$phenos)))
    write.table(x=tmp, file="list_genotypes.txt", quote=FALSE,
                row.names=FALSE, col.names=FALSE)
  } else if(geno.format == "vcf"){
    message("ERROR: saving genotypes in VCF not yet implemented")
  } else if(geno.format == "impute"){
    message("ERROR: saving genotypes in IMPUTE not yet implemented")
  }
  
  for(s in 1:length(data$phenos))
    write.table(x=data$phenos[[s]],
                file=gzfile(paste0("phenotypes_",names(data$phenos)[s],".txt.gz")),
                quote=FALSE)
  tmp <- data.frame(subgroup=names(data$phenos),
                    file=paste0("phenotypes_",names(data$phenos),".txt.gz"))
  write.table(x=tmp, file="list_phenotypes.txt", quote=FALSE,
              row.names=FALSE, col.names=FALSE)
  
  tmp <- rbind(c("id", data$inds$id),
               c("sex", data$inds$sex))
  write.table(x=tmp, file="covariates.txt", quote=FALSE,
              row.names=FALSE, col.names=FALSE)
  tmp <- data.frame(subgroup=names(data$phenos),
                    file=rep("covariates.txt", length(data$phenos)))
  write.table(x=tmp, file="list_covariates.txt", quote=FALSE,
              row.names=FALSE, col.names=FALSE)
  
  write.table(x=data$truth, file=gzfile("truth.txt.gz"),
              quote=FALSE, row.names=FALSE)
}

## Calculate the summary statistics in each tissue
calcSstatsOnSimulatedData <- function(data=NULL){
  stopifnot(! is.null(data))
  sstats <- lapply(data$phenos, function(x){
    tmp <- data.frame(ftr=rep(NA, data$params$nb.pairs), snp=NA, maf=NA, n=NA,
                      pve=NA, sigmahat=NA, betahat.geno=NA, sebetahat.geno=NA,
                      betapval.geno=NA, stringsAsFactors=FALSE)
    if("sex" %in% names(data$inds)){
      tmp$betahat.sex <- NA
      tmp$sebetahat.sex <- NA
      tmp$betapval.sex <- NA
    }
    tmp
  })
  i <- 0
  for(g in 1:nrow(data$gene.coords)){
    for(p in 1:nrow(data$snp.coords)){
      if(data$snp.coords$chr[p] != data$gene.coords$chr[g] ||
         data$snp.coords$start[p] < data$gene.coords$start[g] - data$params$len.cis ||
         data$snp.coords$start[p] > data$gene.coords$start[g] + data$params$len.cis)
        next # SNP not in cis
      i <- i + 1
      gene <- data$gene.coords$id[g]
      for(s in 1:length(data$phenos)){
        sstats[[s]]$ftr[i] <- gene
        sstats[[s]]$snp[i] <- data$snp.coords$id[p]
        sstats[[s]]$maf[i] <- sum(as.numeric(data$geno.counts[p,])) /
          (2 * length(data$phenos[[s]][g,]))
        if(gene %in% rownames(data$phenos[[s]])){
          sstats[[s]]$n[i] <- length(data$phenos[[s]][g,])
          if(! "sex" %in% names(data$inds)){
            tmp <- summary(lm(as.numeric(data$phenos[[s]][g,]) ~
                              as.numeric(data$geno.counts[p,])))
          } else
            tmp <- summary(lm(as.numeric(data$phenos[[s]][g,]) ~
                              as.numeric(data$geno.counts[p,]) +
                              as.numeric(data$inds$sex)))
          sstats[[s]]$pve[i] <- tmp$r.squared
          sstats[[s]]$sigmahat[i] <- tmp$sigma
          sstats[[s]]$betahat.geno[i] <- tmp$coefficients[2,"Estimate"]
          sstats[[s]]$sebetahat.geno[i] <- tmp$coefficients[2,"Std. Error"]
          sstats[[s]]$betapval.geno[i] <- tmp$coefficients[2,"Pr(>|t|)"]
          if("sex" %in% names(data$inds)){
            sstats[[s]]$betahat.sex[i] <- tmp$coefficients[3,"Estimate"]
            sstats[[s]]$sebetahat.sex[i] <- tmp$coefficients[3,"Std. Error"]
            sstats[[s]]$betapval.sex[i] <- tmp$coefficients[3,"Pr(>|t|)"]
          }
        } else # case where gene is absent in subgroup
          sstats[[s]]$n[i] <- 0
      }
    }
  }
  return(sstats)
}

getStdSstatsAndCorrSmallSampleSize <- function(data, sstats, g, p,
                                               correct=TRUE){
  std.sstats.corr <- do.call(rbind, lapply(sstats, function(x){
    N <- x[x$ftr == data$gene.coords$id[g] &
           x$snp == data$snp.coords$id[p],
           "n"]
    if(N == 0){
      c(NA, NA, NA)
    } else{
      sigmahat <- x[x$ftr == data$gene.coords$id[g] &
                    x$snp == data$snp.coords$id[p],
                    6]
      betahat <- x[x$ftr == data$gene.coords$id[g] &
                   x$snp == data$snp.coords$id[p],
                   7]
      sebetahat <- x[x$ftr == data$gene.coords$id[g] &
                     x$snp == data$snp.coords$id[p],
                     8]
      bhat <- betahat / sigmahat
      sebhat <- sebetahat / sigmahat
      t <- qnorm(pt(-abs(bhat/sebhat), N-2, log=TRUE), log=TRUE)
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

## Calculate the log10(ABF) of Wen & Stephens
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
        l10abfs[config,] <- rep(1, ncol(gridS))
    }
  }
  
  return(l10abfs)
}

calcRawAbfsOnSimulatedData <- function(data=NULL, sstats=NULL, grids=NULL){
  stopifnot(! is.null(data), ! is.null(sstats), ! is.null(grids))
  l10abfs.raw <- data.frame(ftr=rep(NA, 10*data$params$nb.pairs),
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
      std.sstats.corr <- getStdSstatsAndCorrSmallSampleSize(data, sstats, g, p)
      
      l10abfs.raw.const.gridL <- calcL10AbfsRawConstGridL(std.sstats.corr,
                                                          grids$gridL)
      for(m in 1:nrow(l10abfs.raw.const.gridL)){
        l10abfs.raw$ftr[i] <- data$gene.coords$id[g]
        l10abfs.raw$snp[i] <- data$snp.coords$id[p]
        l10abfs.raw$config[i] <- rownames(l10abfs.raw.const.gridL)[m]
        l10abfs.raw[i,-c(1:3)] <- l10abfs.raw.const.gridL[m,]
        i <- i + 1
      }
      
      l10abfs.raw.all.gridS <- calcL10AbfsRawAllGridS(std.sstats.corr,
                                                      grids$gridS,
                                                      data$configs)
      for(m in 1:nrow(l10abfs.raw.all.gridS)){
        l10abfs.raw$ftr[i] <- data$gene.coords$id[g]
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
  l10abfs.avg$l10abf.sin[i] <- log10WeightedSum(tmp)
  
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

calcAvgAbfsOnSimulatedData <- function(data=NULL, l10abfs.raw=NULL){
  stopifnot(! is.null(data), ! is.null(l10abfs.raw))
  l10abfs.avg <- data.frame(ftr=rep(NA, data$params$nb.pairs), snp=NA,
                            nb.subgroups=NA, nb.samples=NA)
  for(i in c("l10abf.gen", "l10abf.gen.fix", "l10abf.gen.maxh",
             "l10abf.sin", "l10abf.gen.sin", "l10abf.all")){
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
      l10abfs.avg$ftr[i] <- data$gene.coords$id[g]
      l10abfs.avg$snp[i] <- data$snp.coords$id[p]
      l10abfs.avg$nb.subgroups[i] <- 0
      l10abfs.avg$nb.samples[i] <- 0
      
      tmp <- l10abfs.raw[l10abfs.raw$ftr == l10abfs.avg$ftr[i] &
                         l10abfs.raw$snp == l10abfs.avg$snp[i] &
                         l10abfs.raw$config == "gen",
                         c(4:ncol(l10abfs.raw))]
      l10abfs.avg$l10abf.gen[i] <- log10WeightedSum(tmp)
      tmp <- l10abfs.raw[l10abfs.raw$ftr == l10abfs.avg$ftr[i] &
                         l10abfs.raw$snp == l10abfs.avg$snp[i] &
                         l10abfs.raw$config == "gen-fix",
                         c(4:ncol(l10abfs.raw))]
      l10abfs.avg$l10abf.gen.fix[i] <- log10WeightedSum(tmp)
      tmp <- l10abfs.raw[l10abfs.raw$ftr == l10abfs.avg$ftr[i] &
                         l10abfs.raw$snp == l10abfs.avg$snp[i] &
                         l10abfs.raw$config == "gen-maxh",
                         c(4:ncol(l10abfs.raw))]
      l10abfs.avg$l10abf.gen.maxh[i] <- log10WeightedSum(tmp)
      
      for(config in configs){
        tmp <- l10abfs.raw[l10abfs.raw$ftr == l10abfs.avg$ftr[i] &
                           l10abfs.raw$snp == l10abfs.avg$snp[i] &
                           l10abfs.raw$config == config,
                           c(4:13)]
        if(length(grep("-", config)) == 0 && sum(tmp == 1.0) != length(tmp)){
          l10abfs.avg$nb.subgroups[i] <- l10abfs.avg$nb.subgroups[i] + 1
          l10abfs.avg$nb.samples[i] <- l10abfs.avg$nb.samples[i] +
            ncol(data$phenos[[as.numeric(config)]])
        }
        l10abfs.avg[i, paste0("l10abf.",config)] <- log10WeightedSum(tmp)
      }
      
      l10abfs.avg <- calcAvgAbfsSinAndGenSin(l10abfs.avg, i)
      
      l10abfs.avg <- calcAvgAbfAll(l10abfs.avg, i, configs)
      
      i <- i + 1
    }
  }
  return(l10abfs.avg)
}

getResultsOnSimulatedData <- function(data=NULL, grids=NULL, verbose=0){
  stopifnot(! is.null(data), ! is.null(grids))
  if(verbose > 0)
    message("get results on simulated data with R ...")
  
  sstats <- calcSstatsOnSimulatedData(data)
  
  l10abfs.raw <- calcRawAbfsOnSimulatedData(data, sstats, grids)
  
  l10abfs.avg <- calcAvgAbfsOnSimulatedData(data, l10abfs.raw)
  
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
    for(j in 5:ncol(res$sstats[[s]]))
      tmp <- cbind(tmp, sprintf(fmt="%.06e", res$sstats[[s]][row.ids,j]))
    colnames(tmp) <- colnames(res$sstats[[s]])
    write.table(x=tmp, file=gzfile(paste0("exp_eqtlbma_sumstats_",s,".txt.gz")),
                quote=FALSE, row.names=FALSE, sep=" ", na="nan")
  }
  
  tmp <- res$l10abfs.raw[,1:3]
  for(j in 4:ncol(res$l10abfs.raw))
    tmp <- cbind(tmp, gsub("NA", "nan", sprintf(fmt="%.06e", res$l10abfs.raw[,j])))
  colnames(tmp) <- colnames(res$l10abfs.raw)
  write.table(x=tmp, file=gzfile("exp_eqtlbma_l10abfs_raw.txt.gz"),
              quote=FALSE, row.names=FALSE, sep=" ", na="nan")
  
  tmp <- res$l10abfs.avg[,1:4]
  for(j in 5:ncol(res$l10abfs.avg))
    tmp <- cbind(tmp, gsub("NA", "nan", sprintf(fmt="%.06e", res$l10abfs.avg[,j])))
  colnames(tmp) <- colnames(res$l10abfs.avg)
  write.table(x=tmp, file=gzfile("exp_eqtlbma_l10abfs_avg-grids.txt.gz"),
              quote=FALSE, row.names=FALSE, sep=" ", na="nan")
}

run <- function(){
  params <- parseArgs()
  if(params$verbose > 0)
    message(paste0("START functional_tests.R (", date(), ")"))
  setwd(params$dir.name)
  data <- getSimulatedData(params$rmvGenesFromSbgrps,
                           params$withCovars,
                           params$verbose)
  writeSimulatedData(data, geno.format="custom", verbose=params$verbose)
  grids <- getGrids()
  writeGrids(grids, params$verbose)
  res <- getResultsOnSimulatedData(data, grids, params$verbose)
  writeResultsOnSimulatedData(res, params$verbose)
  if(params$verbose > 0)
    message(paste0("END functional_tests.R (", date(), ")"))
}

##-----------------------------------------------------------------------------
## Run the program

run()
