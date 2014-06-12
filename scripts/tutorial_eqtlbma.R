#!/usr/bin/env Rscript

## `tutorial_eqtlbma.R' simulates data for the tutorial of the eQtlBma package
## Copyright (C) 2013-2014 Timothée Flutre
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

rm(list=ls())
prog.name <- "tutorial_eqtlbma.R"
prog.version <- "1.0"

R.v.maj <- as.numeric(R.version$major)
R.v.min.1 <- as.numeric(strsplit(R.version$minor, "\\.")[[1]][1])
if(R.v.maj < 3 || (R.v.maj == 3 && R.v.min.1 < 0))
    stop("require R >= 3.0 (for built-in mclapply)", call.=FALSE)

suppressPackageStartupMessages(require("parallel")) # for built-in mclapply
suppressPackageStartupMessages(require("GenomicRanges")) # from Bioconductor
suppressPackageStartupMessages(require("MASS")) # for mvrnorm

##' Display the help on stdout
##'
##' The format complies with help2man (http://www.gnu.org/s/help2man)
##' @title Help
help <- function(){
    txt <- paste0("`", prog.name, "' simulate data for the tutorial of the eQtlBma package.\n")
    txt <- paste0(txt, "\n")
    txt <- paste0(txt, "Usage: ", prog.name, " [OPTIONS] ...\n")
    txt <- paste0(txt, "\n")
    txt <- paste0(txt, "Options:\n")
    txt <- paste0(txt, " -h, --help\tdisplay the help and exit\n")
    txt <- paste0(txt, " -V, --version\toutput version information and exit\n")
    txt <- paste0(txt, " -v, --verbose\tverbosity level (0/default=1/2)\n")
    txt <- paste0(txt, "     --pkg\tabsolute path to the package\n")
    txt <- paste0(txt, "     --nsbgrps\tnb of tissues (default=3)\n")
    txt <- paste0(txt, "     --ninds\tnb of individuals (default=200)\n")
    txt <- paste0(txt, "\t\tthey are diploid (allele dose in {0,1,2})\n")
    txt <- paste0(txt, "\t\tsame individuals in all subgroups\n")
    txt <- paste0(txt, "     --ngenes\tnb of genes (default=1000)\n")
    txt <- paste0(txt, "     --agl\taverage gene length (default=10000)\n")
    txt <- paste0(txt, "     --ail\taverage intergenic length (default=50000)\n")
    txt <- paste0(txt, "     --nsnps\tnb of SNPs (default=10000)\n")
    txt <- paste0(txt, "\t\tuniformly distributed along the genome\n")
    txt <- paste0(txt, "     --anchor\tanchor for cis region (default=TSS/TSS+TES)\n")
    txt <- paste0(txt, "     --cr5\tradius of cis region in 5' (default=1000)\n")
    txt <- paste0(txt, "     --cr3\tradius of cis region in 3' (default=1000)\n")
    txt <- paste0(txt, "     --fsg\tfixed nb of cis SNPs per gene (or use --asg)\n")
    txt <- paste0(txt, "     --asg\taverage nb of cis SNPs per gene (default=50)\n")
    txt <- paste0(txt, "     --maf\tminor allele frequency (default=0.3)\n")
    txt <- paste0(txt, "     --rare\tproportion of SNPs with rare alleles (with MAF=0.02, default=0.1)\n")
    ## txt <- paste0(txt, "     --related\tsimulate genotypes for related individuals\n")
    ## txt <- paste0(txt, "\t\tusing the \"popgen\" package from J. Marchini\n")
    ## txt <- paste0(txt, "\t\tsome SNPs can be fixed\n")
    txt <- paste0(txt, "     --pi0\tprior proba for a gene to have no eQTL in any subgroup (default=0.3)\n")
    txt <- paste0(txt, "     --seed\tseed for the RNG (default=1859)\n")
    txt <- paste0(txt, "     --dir\tdirectory in which files are written (current by default)\n")
    txt <- paste0(txt, "     --ncores\tnb of cores to run in parallel (default=1)\n")
    txt <- paste0(txt, "\n")
    txt <- paste0(txt, "Examples:\n")
    txt <- paste0(txt, " Rscript ./", prog.name, " --pkg ~/src/eqtlbma\n")
    message(txt)
}

##' Display version and license information on stdout
##'
##' To comply with help2man (http://www.gnu.org/s/help2man)
##' @title Version
version <- function(){
    txt <- paste0(prog.name, " ", prog.version, "\n")
    txt <- paste0(txt, "\n")
    txt <- paste0(txt, "Copyright (C) 2013-2014 Timothée Flutre.\n")
    txt <- paste0(txt, "License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>\n")
    txt <- paste0(txt, "This is free software; see the source for copying conditions. There is NO\n")
    txt <- paste0(txt, "warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n")
    txt <- paste0(txt, "\n")
    txt <- paste0(txt, "Written by Timothée Flutre.\n")
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
        else if(args[i] == "-v" || args[i] == "--verbose"){
            params$verbose <- as.numeric(args[i+1])
            i <- i + 1
        }
        else if(args[i] == "--pkg"){
            params$pkg <- args[i+1]
            i <- i + 1
        }
        else if(args[i] == "--nsbgrps"){
            params$nb.subgroups <- as.numeric(args[i+1])
            i <- i + 1
        }
        else if(args[i] == "--ninds"){
            params$nb.inds <- as.numeric(args[i+1])
            i <- i + 1
        }
        else if(args[i] == "--ngenes"){
            params$nb.genes <- as.numeric(args[i+1])
            i <- i + 1
        }
        else if(args[i] == "--agl"){
            params$avg.gene.length <- as.numeric(args[i+1])
            i <- i + 1
        }
        else if(args[i] == "--ail"){
            params$avg.intergenic.length <- as.numeric(args[i+1])
            i <- i + 1
        }
        else if(args[i] == "--anchor"){
            params$anchor <- args[i+1]
            i <- i + 1
        }
        else if(args[i] == "--cr5"){
            params$cis.radius.5p <- as.numeric(args[i+1])
            i <- i + 1
        }
        else if(args[i] == "--cr3"){
            params$cis.radius.3p <- as.numeric(args[i+1])
            i <- i + 1
        }
        else if(args[i] == "--fsg"){
            params$fixed.nb.snps.per.gene <- as.numeric(args[i+1])
            i <- i + 1
        }
        else if(args[i] == "--asg"){
            params$avg.nb.snps.per.gene <- as.numeric(args[i+1])
            i <- i + 1
        }
        else if(args[i] == "--maf"){
            params$maf <- as.numeric(args[i+1])
            i <- i + 1
        }
        else if(args[i] == "--rare"){
            params$prop.rare <- as.numeric(args[i+1])
            i <- i + 1
        }
        else if(args[i] == "--related"){
            params$related <- TRUE
            i <- i + 1
        }
        else if(args[i] == "--pi0"){
            params$pi0 <- as.numeric(args[i+1])
            i <- i + 1
        }
        else if(args[i] == "--seed"){
            params$seed <- as.numeric(args[i+1])
            i <- i + 1
        }
        else if(args[i] == "--dir"){
            params$dir <- args[i+1]
            i <- i + 1
        }
        else if(args[i] == "--ncores"){
            params$nb.cores <- as.numeric(args[i+1])
            i <- i + 1
        }
        else{
            write(paste0(prog.name, ": invalid option -- ", args[i], "\n"), stderr())
            help()
            quit("no", status=1)
        }
    }
    
    return(params)
}

##' Check the values of the command-line parameters
##'
##' @param params list of parameters
checkParams <- function(params){
    if(is.null(params$pkg)){
        write(paste0(prog.name, ": option --pkg is required\n"), stderr())
        help()
        quit("no", status=1)
    }
    if(! file.exists(params$pkg)){
        write(paste0(prog.name, ": path ", params$pkg, " doesn't exist\n"),
              stderr())
        quit("no", status=1)
    }
    source(paste0(params$pkg, "/scripts/utils_eqtlbma.R"))
    if(params$prop.rare < 0 || params$prop.rare > 1){
        write(paste0(prog.name, ": --rare ", params$prop.rare,
                     " should be between 0 and 1\n"),
              stderr())
        quit("no", status=1)
    }
    
    if(params$related){
        require(popgen)
        library(MASS)
    }
    
    return(params)
}

##' Simulate individuals
##'
##' @title Individuals
##' @param nb.inds number of individuals
##' @param verbose verbosity level (0/default=1/2)
##' @return Data.frame with individuals in rows and name/sex in columns
simulIndividuals <- function(nb.inds, verbose=1){
    if(verbose > 0)
        message(paste0("simulate ", nb.inds, " individuals ..."))
    
    inds <- data.frame(name=paste0("ind", 1:nb.inds),
                       sex=sample(c(0, 1), nb.inds, replace=TRUE),
                       stringsAsFactors=FALSE)
    
    return(inds)
}

##' Simulate subgroups
##'
##' @title Subgroups
##' @param nb.subgroups integer
##' @param inds data.frame with individuals in rows and name/sex in columns
##' @param verbose verbosity level (0/default=1/2)
##' @return List with one data.frame per subgroup
simulSubgroups <- function(nb.subgroups, inds, verbose=1){
    if(verbose > 0)
        message(paste0("simulate ", nb.subgroups, " subgroups ..."))
    
    subgroups <- lapply(1:nb.subgroups, function(s){
        data.frame(ind=inds$name)
    })
    names(subgroups) <- paste0("tissue", 1:nb.subgroups)
    
    return(subgroups)
}

##' Simulate gene coordinates
##'
##' In the BED format, all on one chromosome
##' @title Coordinates
##' @param nb.genes number of genes
##' @param avg.gene.length average gene length
##' @param avg.intergenic.length average length of intergenic regions
##' @param verbose verbosity level (0/default=1/2)
##' @return Data.frame with gene in rows
simulGeneCoordinates <- function(nb.genes, avg.gene.length,
                                 avg.intergenic.length, verbose=1){
    if(verbose > 0)
        message(paste0("simulate coordinates for ", nb.genes, " genes ..."))
    
    ## explore negative binomial distribution for various parameters
    ## http://www.math.uah.edu/stat/bernoulli/NegativeBinomial.html
    debug <- FALSE
    if(debug){
        E.X <- function(n, p){(n * (1 - p)) / p}
        V.X <- function(n, p){(n * (1 - p)) / p^2}
        S.X <- function(n, p){(2 - p) / sqrt(n * (1 - p))}
        n <- 10^(1); p <- 10^(-3)
        E.X(n, p); V.X(n, p); S.X(n, p)
        tmp <- rnbinom(n=1000, size=n, prob=p) 
        summary(tmp); IQR(tmp)
        hist(x=tmp, main=paste0("n=", format(n,digits=2,scientific=TRUE), " p=", format(p,digits=2,scientific=TRUE), " min=", format(min(tmp),digits=2,scientific=TRUE), " max=", format(max(tmp),digits=2,scientific=TRUE)), breaks=100, xlab="gene lengths (bp)")
        abline(v=E.X(n,p), lwd=2)
        abline(v=E.X(n,p) - 2 * sqrt(V.X(n,p)), lwd=2, lty=2)
        abline(v=E.X(n,p) + 2 * sqrt(V.X(n,p)), lwd=2, lty=2)
    }
    
    ## X ~ NB(n,p) is "gene length"
    E.x <- avg.gene.length
    S.x <- 1.0
    p <- ((S.x^2 * E.x + 4) - sqrt((S.x^2 * E.x + 4)^2 - 16)) / 2
    n <- (E.x * p) / (1 - p)
    gene.lengths <- rnbinom(n=nb.genes, size=n, prob=p)
    if(verbose > 0)
        message(paste0("gene lengths: mean=", format(mean(gene.lengths), scientific=TRUE, digits=2),
                       " sd=", format(sd(gene.lengths), scientific=TRUE, digits=2),
                       " min=", format(min(gene.lengths), scientific=TRUE, digits=2),
                       " med=", format(median(gene.lengths), scientific=TRUE, digits=2),
                       " max=", format(max(gene.lengths), scientific=TRUE, digits=2)))
    
    ## X ~ NB(n,p) is "intergenic length"
    E.x <- avg.intergenic.length
    S.x <- 1.2
    p <- ((S.x^2 * E.x + 4) - sqrt((S.x^2 * E.x + 4)^2 - 16)) / 2
    n <- (E.x * p) / (1 - p)
    intergenic.lengths <- rnbinom(n=nb.genes, size=n, prob=p)
    if(verbose > 0)
        message(paste0("intergenic lengths: mean=", format(mean(intergenic.lengths), scientific=TRUE, digits=2),
                       " sd=", format(sd(intergenic.lengths), scientific=TRUE, digits=2),
                       " min=", format(min(intergenic.lengths), scientific=TRUE, digits=2),
                       " med=", format(median(intergenic.lengths), scientific=TRUE, digits=2),
                       " max=", format(max(intergenic.lengths), scientific=TRUE, digits=2)))
    
    interval.lengths <- c(rbind(intergenic.lengths, gene.lengths))
    
    gene.coords.bed <-
        data.frame(chr=rep("chr1", nb.genes),
                   start=cumsum(interval.lengths)[seq(from=1,
                       to=length(interval.lengths)-1, by=2)],
                   end=cumsum(interval.lengths)[seq(from=2,
                       to=length(interval.lengths), by=2)],
                   name=paste0("gene", 1:nb.genes),
                   score=rep(1000, nb.genes),
                   strand=rep("+", nb.genes),
                   ## strand=sample(x=c("+", "-"), size=nb.genes, replace=TRUE),
                   stringsAsFactors=FALSE)
    ## hist(x=gene.coords.bed$end - gene.coords.bed$start,
    ##      main="Gene lengths (in bp)", xlab="")
    
    return(gene.coords.bed)
}

##' Simulate SNP coordinates
##'
##' In the BED format, all on one chromosome
##' @title Coordinates
##' @param gene.coords.bed data.frame with genes in rows
##' @param anchor gene boundary(ies) for the cis region (TSS or TSS+TES)
##' @param cis.radius.5p length of the 5' half of the cis region (in bp)
##' @param cis.radius.3p length of the 3' half of the cis region (in bp)
##' @param fixed.nb.snps.per.gene scalar
##' @param avg.nb.snps.per.gene scalar, used for negative binomial
##' @return Data.frame with SNP in rows
simulSnpCoordinates <- function(gene.coords.bed, anchor, cis.radius.5p,
                                cis.radius.3p, fixed.nb.snps.per.gene=NULL,
                                avg.nb.snps.per.gene=NULL, verbose=1){
    stopifnot(anchor %in% c("TSS", "TSS+TES"),
              xor(is.null(fixed.nb.snps.per.gene),
                  is.null(avg.nb.snps.per.gene)))
    
    if(verbose > 0)
        message("simulate SNP coordinates ...")
    
    nb.genes <- nrow(gene.coords.bed)
    
    if(! is.null(fixed.nb.snps.per.gene)){
        nb.cis.snps.per.gene <- rep(x=fixed.nb.snps.per.gene, times=nb.genes)
    } else{ # draw nb of cis SNPs per gene from a negative binomial
        E.x <- avg.nb.snps.per.gene
        S.x <- 1 # skewness
        p <- ((S.x^2 * E.x + 4) - sqrt((S.x^2 * E.x + 4)^2 - 16)) / 2
        n <- (E.x * p) / (1 - p)
        nb.cis.snps.per.gene <- rnbinom(n=nb.genes, size=n, prob=p)
        if(verbose > 0)
            message(paste0("total nb of SNPs (in cis of at least one gene): ",
                           sum(nb.cis.snps.per.gene)))
        ## summary(nb.cis.snps.per.gene)
        ## sum(nb.cis.snps.per.gene == 0)
    }
    
    snp.loci <- do.call(c, lapply(1:nb.genes, function(g){
        if(gene.coords.bed$strand[g] == "+"){
            coord.5p <- max(1, gene.coords.bed$start[g] + 1 - cis.radius.5p)
            if(anchor == "TSS"){
                coord.3p <- gene.coords.bed$start[g] + 1 + cis.radius.3p
            } else # TSS+TES
                coord.3p <- gene.coords.bed$end[g] + cis.radius.3p
            sample(x=seq(from=coord.5p, to=coord.3p, by=1),
                   size=nb.cis.snps.per.gene[g])
        } else{ # "-" strand
            coord.5p <- gene.coords.bed$end[g] + cis.radius.5p
            if(anchor == "TSS"){
                coord.3p <- gene.coords.bed$end[g] - cis.radius.3p
            } else # TSS+TES
                coord.3p <- max(1, gene.coords.bed$start[g] + 1 - cis.radius.3p)
            sample(x=seq(from=coord.3p, to=coord.5p, by=1),
                   size=nb.cis.snps.per.gene[g])
        }
    }))
    
    snp.coords.bed <- data.frame(chr="chr1",
                                 start=snp.loci - 1,
                                 end=snp.loci,
                                 name=paste0("snp", 1:length(snp.loci)),
                                 stringsAsFactors=FALSE)
    
    return(snp.coords.bed)
}

##' Simulate genotypes
##'
##' @title Genotypes
##' @param snp.coords.bed data.frame with SNPs in rows
##' @param inds data.frame with individuals in rows
##' @param maf minor allele frequency (except for rare alleles)
##' @param prop.rare proportion of SNPs with rare alleles (MAF=0.02)
##' @param related boolean
##' @param verbose verbosity level (0/default=1/2/3)
##' @return Matrix of genotypes with SNPs in rows and individuals in columns
simulGenotypes <- function(snp.coords.bed, inds, maf, prop.rare, related,
                           verbose=1){
    if(verbose > 0)
        message("simulate genotypes (maf=", maf, " rare=", prop.rare, ") ...")
    
    nb.snps <- nrow(snp.coords.bed)
    nb.inds <- nrow(inds)
    
    if(! related){
        freq.hwe <- function(maf){ # Hardy-Weinberg equilibrium
            f2 <- maf^2 # proba of having 2 copies of minor allele
            f1 <- 2 * maf * (1-maf) # one copy
            f0 <- 1 - f2 - f1 # zero copy
            return(c(f0, f1, f2))
        }
        nb.snps.rare <- floor(nb.snps * prop.rare)
        nb.snps.common <- nb.snps - nb.snps.rare
        genos.dose <-
            matrix(data=replicate(nb.snps.common,
                       sample(x=0:2, size=nb.inds, replace=TRUE,
                              prob=freq.hwe(maf))),
                   nrow=nb.snps.common, ncol=nb.inds, byrow=TRUE,
                   dimnames=list(snp=snp.coords.bed$name[1:nb.snps.common],
                       ind=inds$name))
        if(nb.snps.rare > 0)
            genos.dose <-
                rbind(genos.dose,
                      matrix(data=replicate(nb.snps.rare,
                                 sample(x=0:2, size=nb.inds,
                                        replace=TRUE,
                                        prob=freq.hwe(0.02))),
                             nrow=nb.snps.rare, ncol=nb.inds,
                             byrow=TRUE,
                             dimnames=list(snp=snp.coords.bed$name[(nb.snps.common+1):nb.snps],
                                 ind=inds$name)))
    } else{
        nb.pops <- 3
        nb.inds.per.pop <- floor(nb.inds / nb.pops)
        nb.snps.per.ind <- floor(nb.snps / nb.inds)
        
        genos <- simMD(N=nb.inds.per.pop, P=nb.pops, L=nb.snps.per.ind,
                       c.vec1=runif(n=nb.pops, min=0.1, max=0.3), c.vec2=1, ac=2, beta=1)
        dimnames(genos) <- list(ind=paste0("ind", 1:dim(genos)[1]),
                                chr=c("chr.pat", "chr.mat"), # paternal and maternal
                                snp=paste0("snp", 1:dim(genos)[3]))
        genos.dose <- apply(X=genos, MARGIN=1, FUN=function(x){ 4 - colSums(x) })
        maf <- function(x){tmp <- sum(x)/(2*length(x)); ifelse(tmp <= 0.5, tmp, 1-tmp)}
        genos.dose <- genos.dose[apply(genos.dose, 1, maf) > 0.05,] # discard rare variants
        dim(genos.dose)
        N <- ncol(genos.dose) # final nb of individuals
        P <- nrow(genos.dose) # final nb of SNPs per individual
        
        ## estimate K as genomic relationship matrix (Habier et al 2007)
        ## http://cran.r-project.org/web/packages/synbreed/vignettes/IntroSyn.pdf
        ## W <- t(genos.dose)
        ## P <- apply(genos.dose, 1, function(x){rep(2*maf(x), N)})
        ## U <- ((W - P) %*% t(W - P)) / (2 * sum(maf(genos.dose[1,]) * (1 - maf(genos.dose[1,]))))
        ## K <- U
        
        ## using a pedigree
        ## library(synbreed)
        ## genos <- simul.pedigree(gener=4, ids=c(3,5,8,8))
        ## library(kinship)
        ## K <- kinship(id=, father.id=, mother.id=)
    }
    
    return(genos.dose)
}

##' Simulate gene expression levels
##'
##' Based on Bayesian multivariate linear regressions
##' @title Expression levels
##' @param subgroups list with one data.frame per subgroup
##' @param inds data.frame with individuals in rows
##' @param genos.dose matrix of genotypes with SNPs in rows and individuals in columns
##' @param gn2sn list of genes, each with a corresponding vector of SNPs
##' @param related boolean
##' @param pi0 prior proba for a gene to have no eQTL in any subgroup
##' @param nb.cores nb of cores for parallel execution (via mclapply)
##' @param verbose verbosity level (0/default=1/2)
##' @return List with expression levels and truth
simulGeneExpLevels <- function(subgroups, inds, genos.dose, gn2sn, related, pi0,
                               nb.cores=1, verbose=1){
    if(verbose > 0)
        message("simulate gene expression levels ...")
    
    nb.subgroups <- length(subgroups)
    nb.inds <- nrow(inds)
    nb.genes <- length(gn2sn)
    
    ## intercepts and other covariates (e.g. sex)
    ## use inverse-Gamma prior on variance
    Sigma.beta.c <- diag(x=1 / rgamma(nb.subgroups, 2, 1),
                         nrow=nb.subgroups, ncol=nb.subgroups)
    
    ## genotype effect sizes (a single SNP)
    ## use "meta-analysis" prior (prior grid of variances)
    oma2.plus.phi2 <- 0.8^2 # prior variance of subgroup SNP effect
    oma2.over.oma2.plus.phi2 <- 3/4 # homogeneity
    oma2 <- oma2.plus.phi2 * oma2.over.oma2.plus.phi2 # prior var of avg SNP effect
    phi2 <- oma2.plus.phi2 * (1 - oma2.over.oma2.plus.phi2) # prior var of subgroup SNP effect, given avg effect
    Sigma.beta.g <- matrix(rep(oma2, nb.subgroups^2), nb.subgroups, nb.subgroups) +
        diag(rep(phi2, nb.subgroups), nb.subgroups, nb.subgroups)
    
    ## errors: covariance between subgroups
    ## use inverse-Wishart prior
    r <- ceiling(0.1 * nb.inds) # "r small relative to sample size"
    q.i <- 2 # intercept + other covariates
    m.i <- ceiling(0.3 * nb.inds) # maybe there is a better choice?
    nu.i <- m.i - q.i - r - 1
    stopifnot(nu.i > 0)
    H.i <- diag(nb.subgroups) # maybe there is a better choice?
    cov.err.S <- solve(rWishart(n=1, df=m.i,
                                Sigma=(1/nu.i)*solve(H.i))[,,1])
    if(verbose > 0){
        message("covariance matrix of the errors (same for all genes):")
        print(cov.err.S)
    }
    
    ## errors: covariance between individuals
    if(related){ # "random" effect, but not tested!!
        K <- estimKinshipBaldingNichols(genos.dose)
        sigma.u <- 1 # prior std dev of "random" genetic effect (genetic relatedness)
        u <- do.call(cbind, lapply(1:nb.genes, function(g){
            mvrnorm(n=1, mu=rep(0,N), Sigma=sigma.u^2 * K)
        }))
        dimnames(u) <- dimnames(e)
    } else
        cov.err.I <- diag(nb.inds)
    
    truth <- data.frame(gene=do.call(c, lapply(names(gn2sn), function(n){rep(n, length(gn2sn[[n]]))})),
                        snp=do.call(c, lapply(names(gn2sn), function(n){gn2sn[[n]]})),
                        config=rep(0, sum(sapply(gn2sn, length))),
                        stringsAsFactors=FALSE)
    truth <- cbind(truth, matrix(0.0, ncol=nb.subgroups))
    colnames(truth)[seq(4,4+nb.subgroups-1)] <- paste0("beta.g.",
                                                       1:nb.subgroups)
    
    configs <- getBinaryConfigs(nb.subgroups) # first row is null config
    prior.configs <- c(0.0, # conditional on being an eQTL in at least one
                       rep((1-0.3)/(nrow(configs)-2), nrow(configs)-2),
                       0.3) # fully active config assumed most frequent
    stopifnot(all.equal(sum(prior.configs[2:length(prior.configs)]), 1.0))
    
    if(nb.cores == 1)
        pb <- txtProgressBar(min=0, max=nb.genes, style=3)
    explevels <- mclapply(1:nb.genes, function(g){
        if(nb.cores == 1)
            setTxtProgressBar(pb, g)
        gene <- names(gn2sn)[g]
        X.c <- cbind(rep(1, nb.inds),
                     inds$sex)
        B.c <- matrix(mvrnorm(n=2, mu=rep(0, nb.subgroups),
                              Sigma=Sigma.beta.c),
                      nrow=2, ncol=nb.subgroups)
        ## having a cov.err.S per gene sometimes creates E with NAs (because of singular V %x% U)
        ## cov.err.S <- solve(rWishart(n=1, df=m.i,
        ##                             Sigma=(1/nu.i)*solve(H.i))[,,1])
        E <- matvrnorm(n=1, M=matrix(0, nrow=nb.inds, ncol=nb.subgroups),
                       U=cov.err.I, V=cov.err.S)[,,1]
        if(sum(is.na(E)) > 0){
            if(nb.cores == 1) print("\n")
            print(cov.err.S)
            stop(paste0("NAs in error matrix for ", gene))
        }
        if(runif(1) < pi0){
            X <- X.c
            B <- B.c
        } else{ # at most one eQTN per gene
            eqtn <- sample(x=gn2sn[[g]], size=1) # TODO: allow prior to depend on annotations
            X.g <- matrix(genos.dose[eqtn,], ncol=1)
            B.g <- matrix(mvrnorm(n=1, mu=rep(0, nb.subgroups),
                                  Sigma=Sigma.beta.g),
                          nrow=1, ncol=nb.subgroups)
            config <- configs[sample(x=nrow(configs), size=1,
                                     prob=prior.configs),]
            B.g[1, which(config == 0)] <- 0
            truth[which(truth$gene == gene & truth$snp == eqtn),
                  "config"] <<- paste0(which(config != 0), collapse="-")
            truth[which(truth$gene == gene & truth$snp == eqtn),
                  seq(4,4+nb.subgroups-1)] <<- B.g
            X <- cbind(X.c, X.g)
            B <- rbind(B.c, B.g)
        }
        Y <- X %*% B + E
        rownames(Y) <- inds$name
        colnames(Y) <- names(subgroups)
        Y
    }, mc.cores=nb.cores)
    names(explevels) <- names(gn2sn)
    if(nb.cores == 1)
        close(pb)
    
    if(verbose > 0){
        message("true configuration proportions:")
        print(table(truth$config[truth$config != "0"]) /
              sum(truth$config != "0"))
    }
    
    return(list(explevels=explevels, truth=truth))
}

##' Get the large and small grids
##'
##' @title Grids
##' @return List
getGrids <- function(){
  grids <- list()
  grids$gridL <- makeGrid(grid.type="general")
  grids$gridS <- makeGrid(grid.type="small")
  return(grids)
}

##' Write simulated data to files
##'
##' @title Write data
##' @param dir 
##' @param inds 
##' @param subgroups 
##' @param gene.coords.bed 
##' @param snp.coords.bed 
##' @param genos.dose 
##' @param explevels.genes 
##' @param grids 
##' @param gn2sn 
##' @param verbose 
writeData <- function(dir, inds, subgroups, gene.coords.bed, snp.coords.bed,
                      genos.dose, explevels.genes, grids, gn2sn, verbose=1){
    if(verbose > 0)
        message("write data ...")
    
    nb.subgroups <- length(subgroups)
    
    p2f <- paste0(dir, "/gene_coords.bed.gz")
    write.table(x=gene.coords.bed, file=gzfile(p2f), quote=F, sep="\t",
                row.names=F, col.names=F)
    
    p2f <- paste0(dir, "/snp_coords.bed.gz")
    write.table(x=snp.coords.bed, file=gzfile(p2f), quote=F, sep="\t",
                row.names=F, col.names=F)
    p2f <- paste0(dir, "/genotypes.txt.gz")
    write.table(x=genos.dose, file=gzfile(p2f), quote=F, sep="\t", row.names=T,
                col.names=T)
    tmp <- data.frame(subgroup=names(subgroups),
                      file=rep(p2f, length(subgroups)))
    p2f <- paste0(dir, "/list_genotypes.txt")
    write.table(x=tmp, file=p2f, quote=F, sep="\t", row.names=F, col.names=F)
    
    tmp <- lapply(1:nb.subgroups, function(s){
        explevels.s <-
            do.call(rbind, lapply(explevels.genes$explevels, function(Y.g){
                Y.g[,s]
            }))
        p2f <- paste0(dir, "/explevels_", names(subgroups)[s],
                      ".txt.gz")
        write.table(x=explevels.s, file=gzfile(p2f), quote=F, sep="\t",
                    row.names=T, col.names=T)
    })
    tmp <- data.frame(subgroup=names(subgroups),
                      file=paste0(dir, "/explevels_", names(subgroups),
                          ".txt.gz"))
    p2f <- paste0(dir, "/list_explevels.txt")
    write.table(x=tmp, file=p2f, quote=F, sep="\t", row.names=F, col.names=F)
    
    p2f <- paste0(dir, "/covariates.txt.gz")
    tmp <- matrix(inds$sex, nrow=1, ncol=nrow(inds),
                  dimnames=list(covariates=c("sex"), inds=inds$name))
    write.table(x=tmp, file=gzfile(p2f), quote=F, sep="\t", row.names=T,
                col.names=T)
    tmp <- data.frame(subgroup=names(subgroups),
                      file=paste0(dir, "/covariates.txt.gz"))
    p2f <- paste0(dir, "/list_covariates.txt")
    write.table(x=tmp, file=p2f, quote=F, sep="\t", row.names=F, col.names=F)
    
    p2f <- paste0(dir, "/grid_phi2_oma2_general.txt")
    write.table(x=grids$gridL, file=p2f, quote=F, sep="\t", row.names=F,
                col.names=F)
    p2f <- paste0(dir, "/grid_phi2_oma2_with-configs.txt")
    write.table(x=grids$gridS, file=p2f, quote=F, sep="\t", row.names=F,
                col.names=F)
    
    p2f <- paste0(dir, "/truth.txt.gz")
    write.table(x=explevels.genes$truth, file=gzfile(p2f), quote=F, sep="\t",
                row.names=F, col.names=T)
    
    pdf("hist_maf.pdf")
    plotHistMinAllelFreq(genos.dose, main="Tutorial of eQtlBma",
                         breaks=50, col="grey", border="white")
    dev.off()
    
    pdf("hist_cis-snps-per-gene.pdf")
    plotHistCisSnpsPerGene(gn2sn)
    dev.off()
}

run <- function(params){
    set.seed(params$seed)
    
    inds <- simulIndividuals(params$nb.inds)
    
    subgroups <- simulSubgroups(params$nb.subgroups, inds)
    
    gene.coords.bed <- simulGeneCoordinates(params$nb.genes,
                                            params$avg.gene.length,
                                            params$avg.intergenic.length,
                                            params$verbose)
    
    snp.coords.bed <- simulSnpCoordinates(gene.coords.bed,
                                          params$anchor,
                                          params$cis.radius.5p,
                                          params$cis.radius.3p,
                                          params$fixed.nb.snps.per.gene,
                                          params$avg.nb.snps.per.gene,
                                          params$verbose)
    
    gn2sn <- findCisSnpsPerGene(gene.coords.bed,
                                snp.coords.bed,
                                params$anchor,
                                params$cis.radius.5p,
                                params$cis.radius.3p,
                                params$verbose)
    
    genos.dose <- simulGenotypes(snp.coords.bed,
                                 inds,
                                 params$maf,
                                 params$prop.rare,
                                 params$related,
                                 params$verbose)
    
    explevels.genes <- simulGeneExpLevels(subgroups,
                                          inds,
                                          genos.dose,
                                          gn2sn,
                                          params$related,
                                          params$pi0,
                                          params$nb.cores,
                                          params$verbose)
    
    grids <- getGrids()
    
    writeData(params$dir, inds, subgroups, gene.coords.bed, snp.coords.bed,
              genos.dose, explevels.genes, grids, gn2sn)
}

main <- function(){
    params <- list(verbose=1,
                   pkg=NULL,
                   nb.subgroups=3,
                   nb.inds=200,
                   nb.genes=10^3,
                   avg.gene.length=10^4,
                   avg.intergenic.length=5*10^4,
                   anchor="TSS",
                   cis.radius.5p=10^3,
                   cis.radius.3p=10^3,
                   fixed.nb.snps.per.gene=NULL,
                   avg.nb.snps.per.gene=50,
                   maf=0.3,
                   prop.rare=0.1,
                   related=FALSE,
                   pi0=0.3,
                   seed=1859,
                   dir=getwd(),
                   nb.cores=1)
    
    params <- parseCmdLine(params)
    
    params <- checkParams(params)
    
    if(params$verbose > 0){
        message(paste0("START ", prog.name, " ",
                       format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
        message(paste0("cwd: ", getwd()))
    }
    
    system.time(run(params))
    
    if(params$verbose > 0){
        message(paste0("END ", prog.name, " ",
                       format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
        ## print(object.size(x=lapply(ls(), get)), units="Kb") # return an error I don't understand
    }
}

if(! interactive())
    main()
