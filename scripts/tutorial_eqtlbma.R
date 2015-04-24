#!/usr/bin/env Rscript

## `tutorial_eqtlbma.R' simulates data for the tutorial of the eQtlBma package
## Copyright (C) 2013-2015 Timothée Flutre
## License: GPL-3+
## Persons: Timothée Flutre [cre,aut]

rm(list=ls())
prog.name <- "tutorial_eqtlbma.R"
prog.version <- "1.3.1a" # same as in configure.ac

R.v.maj <- as.numeric(R.version$major)
R.v.min.1 <- as.numeric(strsplit(R.version$minor, "\\.")[[1]][1])
if(R.v.maj < 3 || (R.v.maj == 3 && R.v.min.1 < 0))
    stop("require R >= 3.0 (for built-in mclapply)", call.=FALSE)

suppressPackageStartupMessages(library("parallel")) # for built-in mclapply
suppressPackageStartupMessages(library("GenomicRanges")) # from Bioconductor
suppressPackageStartupMessages(library("MASS")) # for mvrnorm

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
    txt <- paste0(txt, "     --ninds\tnb of individuals per subgroup (default=200)\n")
    txt <- paste0(txt, "\t\tindividuals are diploid\n")
    txt <- paste0(txt, "\t\ta single number, say 200, means same individuals in all subgroups\n")
    txt <- paste0(txt, "\t\tseveral numbers, say 200-150-200, means different individuals between subgroups (must agree with --nsbgrps)\n")
    txt <- paste0(txt, "     --ngenes\tnb of genes (default=1000)\n")
    txt <- paste0(txt, "     --nchrs\tnb of chromosome(s) (default=1)\n")
    txt <- paste0(txt, "     --agl\taverage gene length (default=10000)\n")
    txt <- paste0(txt, "     --ail\taverage intergenic length (default=50000)\n")
    txt <- paste0(txt, "     --anchor\tanchor for cis region (default=TSS/TSS+TES)\n")
    txt <- paste0(txt, "     --cr5\tradius of cis region in 5' (default=1000)\n")
    txt <- paste0(txt, "     --cr3\tradius of cis region in 3' (default=1000)\n")
    txt <- paste0(txt, "     --fsg\tfixed nb of cis SNPs per gene (or use --asg)\n")
    txt <- paste0(txt, "     --asg\taverage nb of cis SNPs per gene (default=50)\n")
    txt <- paste0(txt, "     --maf\tminor allele frequency (default=0.3)\n")
    txt <- paste0(txt, "     --rare\tproportion of SNPs with rare alleles (with MAF=0.02, default=0.0)\n")
    txt <- paste0(txt, "     --pi0\tprior proba for a gene to have no eQTL in any subgroup (default=0.3)\n")
    txt <- paste0(txt, "     --coverr\terror covariance between subgroups (default=1)\n")
    txt <- paste0(txt, "\t\t0: the SxS covariance matrix is diagonal (usually the case if different individuals between subgroups), same for all genes\n")
    txt <- paste0(txt, "\t\t1: the SxS covariance matrix is unconstrained (usually the case if same individuals in all subgroups), same for all genes\n")
    ## txt <- paste0(txt, "\t\t2: the SxS covariance matrix is unconstrained, different for each gene\n"
    txt <- paste0(txt, "     --seed\tseed for the RNG (default=1859)\n")
    txt <- paste0(txt, "     --dir\tdirectory in which files are written (current by default)\n")
    txt <- paste0(txt, "     --ncores\tnb of cores to run in parallel (default=1)\n")
    txt <- paste0(txt, "\n")
    txt <- paste0(txt, "Examples:\n")
    txt <- paste0(txt, " Rscript ./", prog.name, " --pkg ~/src/eqtlbma\n")
    txt <- paste0(txt, "\n")
    txt <- paste0(txt, "Report bugs to <eqtlbma-users@googlegroups.com>.")
    write(txt, stdout())
}

##' Display version and license information on stdout
##'
##' To comply with help2man (http://www.gnu.org/s/help2man)
##' @title Version
version <- function(){
    txt <- paste0(prog.name, " ", prog.version, "\n")
    txt <- paste0(txt, "\n")
    txt <- paste0(txt, "Copyright (C) 2013-2015 Timothée Flutre.\n")
    txt <- paste0(txt, "License GPL-3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>\n")
    txt <- paste0(txt, "This is free software; see the source for copying conditions. There is NO\n")
    txt <- paste0(txt, "warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n")
    txt <- paste0(txt, "\n")
    txt <- paste0(txt, "Written by Timothée Flutre [cre,aut].")
    write(txt, stdout())
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
            params$nb.inds <- as.numeric(strsplit(args[i+1], "-")[[1]])
            i <- i + 1
        }
        else if(args[i] == "--ngenes"){
            params$nb.genes <- as.numeric(args[i+1])
            i <- i + 1
        }
        else if(args[i] == "--nchrs"){
            params$nb.chrs <- as.numeric(args[i+1])
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
        else if(args[i] == "--pi0"){
            params$pi0 <- as.numeric(args[i+1])
            i <- i + 1
        }
        else if(args[i] == "--coverr"){
            params$coverr <- as.numeric(args[i+1])
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
        help()
        quit("no", status=1)
    }
    source(paste0(params$pkg, "/scripts/utils_eqtlbma.R"))

    if(length(params$nb.inds) > 1 &
       length(params$nb.inds) != params$nb.subgroups){
        write(paste0(prog.name, ": --ninds should correspond to --nsbgrps\n"),
              stderr())
        help()
        quit("no", status=1)
    }

    if(! is.null(params$fixed.nb.snps.per.gene))
        params$avg.nb.snps.per.gene <- NULL

    if(! params$coverr %in% c(0,1,2)){
        write(paste0(prog.name, ": --coverr should be 0, 1 or 2\n"),
              stderr())
        help()
        quit("no", status=1)
    }
    if(params$coverr != 0 & length(params$nb.inds) > 1){
        write(paste0(prog.name, ": different individuals between",
                     " subgroups should go with --coverr 0\n"),
              stdout())
        help()
        quit("no", status=1)
    }

    if(params$prop.rare < 0 || params$prop.rare > 1){
        write(paste0(prog.name, ": --rare ", params$prop.rare,
                     " should be between 0 and 1\n"),
              stderr())
        help()
        quit("no", status=1)
    }

    return(params)
}

##' Simulate individuals
##'
##' @title Individuals
##' @param nb.inds number of individuals
##' @param nb.subgroups integer
##' @param verbose verbosity level (0/default=1/2)
##' @return List with one data.frame per subgroup with individuals in rows and name/sex in columns
simulIndividuals <- function(nb.inds, nb.subgroups, verbose=1){
    tot.nb.inds <- ifelse(length(nb.inds) == 1, nb.inds, sum(nb.inds))
    if(verbose > 0){
        msg <- paste0("simulate ", tot.nb.inds, " individuals and ",
                      nb.subgroups, " subgroup",
                      ifelse(nb.subgroups > 1, "s", ""))
        if(length(nb.inds) == 1){
            msg <- paste0(msg, " (same individuals in all subgroups) ...")
        } else
            msg <- paste0(msg, " (different individuals between subgroups) ...")
        write(msg, stdout())
    }

    inds <- list("tissue1"=data.frame(name=paste0("ind", 1:nb.inds[1]),
                     sex=sample(c(0, 1), nb.inds[1], replace=TRUE),
                     stringsAsFactors=FALSE))
    if(nb.subgroups > 1){
      if(length(nb.inds) == 1){
        for(s in 2:nb.subgroups)
          inds[[as.character(s)]] <- inds[[1]]
        names(inds) <- paste0("tissue", 1:nb.subgroups)
      } else{
        for(s in 2:nb.subgroups)
          inds[[as.character(s)]] <-
            data.frame(name=paste0("ind", (sum(nb.inds[1:(s-1)])+1):
                                     (sum(nb.inds[1:(s-1)])+nb.inds[s])),
                       sex=sample(c(0, 1), nb.inds[s], replace=TRUE),
                       stringsAsFactors=FALSE)
        names(inds) <- paste0("pop", 1:nb.subgroups)
      }
    }

    return(inds)
}

##' Simulate gene coordinates
##'
##' In the BED format, all on one chromosome
##' @title Coordinates
##' @param nb.genes number of genes
##' @param nb.chrs number of chromosomes
##' @param avg.gene.length average gene length
##' @param avg.intergenic.length average length of intergenic regions
##' @param verbose verbosity level (0/default=1/2)
##' @return Data.frame with gene in rows
simulGeneCoordinates <- function(nb.genes, nb.chrs, avg.gene.length,
                                 avg.intergenic.length, verbose=1){
    if(verbose > 0)
        write(paste0("simulate coordinates for ", nb.genes, " gene",
                     ifelse(nb.genes > 1, "s", ""), " (",
                     nb.chrs, " chr", ifelse(nb.chrs > 1, "s", ""), ") ..."),
              stdout())

    chr.names <- sprintf(paste0("chr%0", floor(log10(nb.chrs))+1, "i"),
                         1:nb.chrs)
    tmp <- round(seq(0, nb.genes, length.out=nb.chrs + 1))
    nb.genes.chrs <- rev(rev(tmp)[-length(tmp)] - rev(tmp)[-1])
    names(nb.genes.chrs) <- chr.names

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
    gene.length.chrs <- lapply(nb.genes.chrs, function(ngc){
        rnbinom(n=ngc, size=n, prob=p)
    })
    names(gene.length.chrs) <- chr.names
    gene.lengths <- do.call(c, gene.length.chrs)
    if(verbose > 0)
        write(paste0("gene lengths: mean=", format(mean(gene.lengths), scientific=TRUE, digits=2),
                     " sd=", format(sd(gene.lengths), scientific=TRUE, digits=2),
                     " min=", format(min(gene.lengths), scientific=TRUE, digits=2),
                     " med=", format(median(gene.lengths), scientific=TRUE, digits=2),
                     " max=", format(max(gene.lengths), scientific=TRUE, digits=2)),
              stdout())

    ## X ~ NB(n,p) is "intergenic length"
    E.x <- avg.intergenic.length
    S.x <- 1.2
    p <- ((S.x^2 * E.x + 4) - sqrt((S.x^2 * E.x + 4)^2 - 16)) / 2
    n <- (E.x * p) / (1 - p)
    intergenic.length.chrs <- lapply(nb.genes.chrs, function(ngc){
        rnbinom(n=ngc, size=n, prob=p)
    })
    names(intergenic.length.chrs) <- chr.names
    intergenic.lengths <- do.call(c, intergenic.length.chrs)
    if(verbose > 0)
        write(paste0("intergenic lengths: mean=", format(mean(intergenic.lengths), scientific=TRUE, digits=2),
                     " sd=", format(sd(intergenic.lengths), scientific=TRUE, digits=2),
                     " min=", format(min(intergenic.lengths), scientific=TRUE, digits=2),
                     " med=", format(median(intergenic.lengths), scientific=TRUE, digits=2),
                     " max=", format(max(intergenic.lengths), scientific=TRUE, digits=2)),
              stdout())

    interval.length.chrs <- lapply(chr.names, function(chr){
        c(rbind(intergenic.length.chrs[[chr]], gene.length.chrs[[chr]]))
    })
    names(interval.length.chrs) <- chr.names

    gene.coords.bed <- do.call(rbind, lapply(chr.names, function(chr){
        idx <- which(chr.names == chr)
        data.frame(chr=rep(chr, nb.genes.chrs[chr]),
                   start=cumsum(interval.length.chrs[[chr]])[seq(from=1,
                       to=length(interval.length.chrs[[chr]])-1, by=2)],
                   end=cumsum(interval.length.chrs[[chr]])[seq(from=2,
                       to=length(interval.length.chrs[[chr]]), by=2)],
                   name=sprintf(paste0("gene%0",
                       floor(log10(nb.genes.chrs[chr]))+1, "i"),
                       (1+c(0, cumsum(nb.genes.chrs))[idx]):
                           (c(0, cumsum(nb.genes.chrs))[idx+1])),
                   score=rep(1000, nb.genes.chrs[chr]),
                   strand=rep("+", nb.genes.chrs[chr]),
                   ## strand=sample(x=c("+", "-"), size=nb.genes.chrs[c],
                   ##     replace=TRUE),
                   stringsAsFactors=FALSE)
    }))
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
        write("simulate SNP coordinates ...", stdout())

    nb.genes <- nrow(gene.coords.bed)
    nb.genes.chrs <- tapply(gene.coords.bed$name,
                            list(as.factor(gene.coords.bed$chr)),
                            length)
    nb.chrs <- length(nb.genes.chrs)
    chr.names <- names(nb.genes.chrs)

    if(! is.null(fixed.nb.snps.per.gene)){
        nb.cis.snps.per.gene <- lapply(names(nb.genes.chrs), function(chr){
            rep(x=fixed.nb.snps.per.gene, times=nb.genes.chrs[chr])
        })
    } else{ # draw nb of cis SNPs per gene from a negative binomial
        E.x <- avg.nb.snps.per.gene
        S.x <- 1 # skewness
        p <- ((S.x^2 * E.x + 4) - sqrt((S.x^2 * E.x + 4)^2 - 16)) / 2
        n <- (E.x * p) / (1 - p)
        nb.cis.snps.per.gene <- lapply(names(nb.genes.chrs), function(chr){
            rnbinom(n=nb.genes.chrs[chr], size=n, prob=p)
        })
    }
    names(nb.cis.snps.per.gene) <- chr.names
    nb.cis.snps.genes <- do.call(c, nb.cis.snps.per.gene)
    nb.cis.snps.per.chr <- sapply(nb.cis.snps.per.gene, sum)
    if(verbose > 0){
        write(paste0("total nb of SNPs (in cis of at least one gene): ",
                     sum(nb.cis.snps.genes)), stdout())
        ## summary(nb.cis.snps.per.gene)
        write(paste0("nb of gene(s) with no cis SNPs: ",
                     sum(nb.cis.snps.genes == 0)), stdout())
    }

    snp.loci <- lapply(names(nb.genes.chrs), function(chr){
        gene.coords.bed.chr <- gene.coords.bed[gene.coords.bed$chr == chr,]
        do.call(c, lapply(1:nb.genes.chrs[chr], function(g){
            if(gene.coords.bed.chr$strand[g] == "+"){
                coord.5p <-
                    max(1, gene.coords.bed.chr$start[g] + 1 - cis.radius.5p)
                if(anchor == "TSS"){
                    coord.3p <-
                        gene.coords.bed.chr$start[g] + 1 + cis.radius.3p
                } else # TSS+TES
                    coord.3p <- gene.coords.bed.chr$end[g] + cis.radius.3p
                sample(x=seq(from=coord.5p, to=coord.3p, by=1),
                       size=nb.cis.snps.per.gene[[chr]][g])
            } else{ # "-" strand
                coord.5p <- gene.coords.bed.chr$end[g] + cis.radius.5p
                if(anchor == "TSS"){
                    coord.3p <- gene.coords.bed.chr$end[g] - cis.radius.3p
                } else # TSS+TES
                    coord.3p <-
                        max(1, gene.coords.bed.chr$start[g] + 1 - cis.radius.3p)
                sample(x=seq(from=coord.3p, to=coord.5p, by=1),
                       size=nb.cis.snps.per.gene[[chr]][g])
            }
        }))
    })
    names(snp.loci) <- chr.names

    snp.coords.bed <- do.call(rbind, lapply(chr.names, function(chr){
        idx <- which(chr.names == chr)
        data.frame(chr=rep(chr, nb.cis.snps.per.chr[chr]),
                   start=snp.loci[[chr]] - 1,
                   end=snp.loci[[chr]],
                   name=sprintf(paste0("snp%0",
                       floor(log10(length(snp.loci[[chr]])))+1, "i"),
                       (1+c(0, cumsum(nb.cis.snps.per.chr))[idx]):
                           (c(0, cumsum(nb.cis.snps.per.chr))[idx+1])),
                   stringsAsFactors=FALSE)
    }))

    return(snp.coords.bed)
}

##' Simulate genotypes
##'
##' @title Genotypes
##' @param snp.coords.bed data.frame with SNPs in rows
##' @param nb.inds vector of integer(s) (1 element if same individuals)
##' @param inds list with, per subgroup, a data.frame with individuals in rows
##' @param maf minor allele frequency (except for rare alleles)
##' @param prop.rare proportion of SNPs with rare alleles (MAF=0.02)
##' @param verbose verbosity level (0/default=1/2/3)
##' @return Matrix of genotypes with SNPs in rows and individuals in columns
simulGenotypes <- function(snp.coords.bed, nb.inds, inds, maf, prop.rare,
                           verbose=1){
    if(verbose > 0)
        write(paste0("simulate genotypes (maf=", maf, " rare=", prop.rare,
                     ") ..."), stdout())

    nb.snps <- nrow(snp.coords.bed)
    tot.nb.inds <- ifelse(length(nb.inds) == 1, nb.inds, sum(nb.inds))
    tot.ind.names <- c()
    if(length(nb.inds) == 1){
        tot.ind.names <- inds[[1]]$name
    } else
        tot.ind.names <- unlist(lapply(inds, function(x){x$name}))

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
                   sample(x=0:2, size=tot.nb.inds, replace=TRUE,
                          prob=freq.hwe(maf))),
               nrow=nb.snps.common, ncol=tot.nb.inds, byrow=TRUE,
               dimnames=list(snp=snp.coords.bed$name[1:nb.snps.common],
                   ind=tot.ind.names))
    if(nb.snps.rare > 0)
        genos.dose <-
            rbind(genos.dose,
                  matrix(data=replicate(nb.snps.rare,
                             sample(x=0:2, size=tot.nb.inds,
                                    replace=TRUE,
                                    prob=freq.hwe(0.02))),
                         nrow=nb.snps.rare, ncol=tot.nb.inds,
                         byrow=TRUE,
                         dimnames=list(snp=snp.coords.bed$name[(nb.snps.common+1):nb.snps],
                             ind=tot.ind.names)))

    return(genos.dose)
}

##' Simulate gene expression levels
##'
##' Based on Bayesian multivariate linear regressions
##' @title Expression levels
##' @param nb.inds vector of integer(s) (1 element if same individuals)
##' @param inds data.frame with individuals in rows
##' @param genos.dose matrix of genotypes with SNPs in rows and individuals in columns
##' @param gn2sn list of genes, each with a corresponding vector of SNPs
##' @param pi0 prior proba for a gene to have no eQTL in any subgroup
##' @param nb.cores nb of cores for parallel execution (via mclapply)
##' @param verbose verbosity level (0/default=1/2)
##' @return List with expression levels and truth
simulGeneExpLevels <- function(nb.inds, inds, genos.dose, gn2sn,
                               pi0, coverr, nb.cores=1, verbose=1){
    if(verbose > 0)
        write("simulate gene expression levels ...", stdout())

    nb.subgroups <- length(inds)
    tot.nb.inds <- ifelse(length(nb.inds) == 1, nb.inds, sum(nb.inds))
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

    ## errors: variance-covariance within-between subgroups
    if(coverr == 0 | coverr == 1){
        r <- ceiling(0.1 * tot.nb.inds) # "r small relative to sample size"
        q.i <- 2 # intercept + other covariates
        m.i <- ceiling(0.3 * tot.nb.inds) # maybe there is a better choice?
        nu.i <- m.i - q.i - r - 1
        stopifnot(nu.i > 0)
        H.i <- diag(nb.subgroups) # maybe there is a better choice?
        cov.err.S <- solve(rWishart(n=1, df=m.i,
                                    Sigma=(1/nu.i)*solve(H.i))[,,1])
        if(coverr == 0) # keps only the values on the diagonal
            cov.err.S <- cov.err.S * (matrix(0,nb.subgroups,nb.subgroups) +
                                      diag(nb.subgroups))
    }
    if(verbose > 0){
        write("covariance matrix of the errors (same for all genes):",
              stdout())
        print(cov.err.S)
    }

    ## errors: covariance between individuals (assumed unrelated)
    cov.err.I <- diag(tot.nb.inds)

    truth <- data.frame(gene=do.call(c, lapply(names(gn2sn), function(n){rep(n, length(gn2sn[[n]]))})),
                        snp=do.call(c, lapply(names(gn2sn), function(n){gn2sn[[n]]})),
                        config=rep(0, sum(sapply(gn2sn, length))),
                        stringsAsFactors=FALSE)
    truth <- cbind(truth, matrix(NA, ncol=nb.subgroups))
    colnames(truth)[seq(4, 4+nb.subgroups-1)] <- paste0("sigma.",
                                                        1:nb.subgroups)
    col.idx.sigmas <- grep("sigma", colnames(truth))
    truth <- cbind(truth, matrix(0.0, ncol=nb.subgroups))
    colnames(truth)[seq(4+nb.subgroups, 4+2*nb.subgroups-1)] <-
        paste0("beta.g.", 1:nb.subgroups)
    col.idx.betas <- grep("beta", colnames(truth))

    configs <- getBinaryConfigs(nb.subgroups) # first row is null config
    prior.configs <- c()
    if(nb.subgroups == 1){
      prior.configs <- c(0.0, 1.0)
    } else
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

        if(length(nb.inds) == 1){ # same individuals in all subgroups
            X.c <- cbind(rep(1, nb.inds[1]),
                         inds[[1]]$sex)
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
            truth[which(truth$gene == gene), col.idx.sigmas] <<-
                rep(sqrt(diag(cov.err.S)), each=sum(truth$gene == gene))

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
                      col.idx.betas] <<- B.g
                X <- cbind(X.c, X.g)
                B <- rbind(B.c, B.g)
            }

            Y <- X %*% B + E
            rownames(Y) <- inds[[1]]$name
            colnames(Y) <- names(inds)
        } else{ # different individuals between subgroups
            config <- rep(0, nb.subgroups)
            if(runif(1) >= pi0){
                eqtn <- sample(x=gn2sn[[g]], size=1)
                B.g <- matrix(mvrnorm(n=1, mu=rep(0, nb.subgroups),
                                      Sigma=Sigma.beta.g),
                              nrow=1, ncol=nb.subgroups)
                config <- configs[sample(x=nrow(configs), size=1,
                                         prob=prior.configs),]
                B.g[1, which(config == 0)] <- 0
                truth[which(truth$gene == gene & truth$snp == eqtn),
                      "config"] <<- paste0(which(config != 0), collapse="-")
                truth[which(truth$gene == gene & truth$snp == eqtn),
                      col.idx.betas] <<- B.g
            }

            Y <- lapply(1:nb.subgroups, function(s){
                X.c.s <- cbind(rep(1, nb.inds[s]),
                               inds[[s]]$sex)
                B.c.s <- matrix(rnorm(n=2, mean=0,
                                      sd=sqrt(diag(Sigma.beta.c)[s])))

                E.s <- matrix(rnorm(n=nb.inds[s], mean=0,
                                    sd=sqrt(diag(cov.err.S)[s])))
                truth[which(truth$gene == gene), paste0("sigma.",s)] <<-
                    sqrt(diag(cov.err.S)[s])

                if(config[s] == 0){
                    X.s <- X.c.s
                    B.s <- B.c.s
                } else{ # eQTN active in this subgroup
                    X.g.s <- matrix(genos.dose[eqtn,inds[[s]]$name], ncol=1)
                    B.g.s <- B.g[1,s]
                    X.s <- cbind(X.c.s, X.g.s)
                    B.s <- rbind(B.c.s, B.g.s)
                }

                Y.s <- X.s %*% B.s + E.s
                rownames(Y.s) <- inds[[s]]$name
                colnames(Y.s) <- names(inds)[s]
                Y.s
            })
        }

        Y
    }, mc.cores=nb.cores)
    names(explevels) <- names(gn2sn)
    if(nb.cores == 1)
        close(pb)

    if(verbose > 0){
        null.genes <- tapply(truth$config, list(truth$gene), function(configs){
            all(configs == "0")
        })
        write(paste0("true pi0: ", sum(null.genes) / length(null.genes)),
              stdout())
        write("true configuration proportions:", stdout())
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
##' @param nb.inds
##' @param inds
##' @param gene.coords.bed
##' @param snp.coords.bed
##' @param genos.dose
##' @param explevels.genes
##' @param grids
##' @param gn2sn
##' @param verbose
writeData <- function(dir, nb.inds, inds, gene.coords.bed, snp.coords.bed,
                      genos.dose, explevels.genes, grids, gn2sn, verbose=1){
    if(verbose > 0)
        write("write data ...", stdout())

    nb.subgroups <- length(inds)

    p2f <- paste0(dir, "/gene_coords.bed.gz")
    write.table(x=gene.coords.bed, file=gzfile(p2f), quote=F, sep="\t",
                row.names=F, col.names=F)

    p2f <- paste0(dir, "/snp_coords.bed.gz")
    write.table(x=snp.coords.bed, file=gzfile(p2f), quote=F, sep="\t",
                row.names=F, col.names=F)

    if(length(nb.inds) == 1){
        p2f <- paste0(dir, "/genotypes.txt.gz")
        write.table(x=genos.dose, file=gzfile(p2f), quote=F, sep="\t",
                    row.names=T, col.names=T)
        tmp <- data.frame(subgroup=names(inds),
                          file=rep(p2f, length(inds)))
    } else{ # different individuals between subgroups
        for(s in 1:nb.subgroups){
            p2f <- paste0(dir, "/genotypes_", names(inds)[s], ".txt.gz")
            write.table(x=genos.dose[,inds[[s]]$name], file=gzfile(p2f),
                        quote=F, sep="\t", row.names=T, col.names=T)
        }
        tmp <- data.frame(subgroup=names(inds),
                          file=paste0(dir, "/genotypes_", names(inds),
                              ".txt.gz"))
    }
    p2f <- paste0(dir, "/list_genotypes.txt")
    write.table(x=tmp, file=p2f, quote=F, sep="\t", row.names=F, col.names=F)

    tmp <- lapply(1:nb.subgroups, function(s){
        explevels.s <-
            do.call(rbind, lapply(explevels.genes$explevels, function(Y.g){
                if(length(nb.inds) == 1){
                    Y.g[,s]
                } else # different individuals between subgroups
                    as.vector(Y.g[[s]])
            }))
        colnames(explevels.s) <- inds[[s]]$name
        p2f <- paste0(dir, "/explevels_", names(inds)[s], ".txt.gz")
        write.table(x=explevels.s, file=gzfile(p2f), quote=F, sep="\t",
                    row.names=T, col.names=T)
    })
    tmp <- data.frame(subgroup=names(inds),
                      file=paste0(dir, "/explevels_", names(inds), ".txt.gz"))
    p2f <- paste0(dir, "/list_explevels.txt")
    write.table(x=tmp, file=p2f, quote=F, sep="\t", row.names=F, col.names=F)

    if(length(nb.inds) == 1){
        p2f <- paste0(dir, "/covariates.txt.gz")
        tmp <- matrix(inds[[1]]$sex, nrow=1, ncol=nrow(inds[[1]]),
                      dimnames=list(covariates=c("sex"), inds=inds[[1]]$name))
        write.table(x=tmp, file=gzfile(p2f), quote=F, sep="\t", row.names=T,
                    col.names=T)
        tmp <- data.frame(subgroup=names(inds),
                          file=rep(p2f, length(inds)))
    } else{ # different individuals between subgroups
        for(s in 1:nb.subgroups){
            p2f <- paste0(dir, "/covariates_", names(inds)[s], ".txt.gz")
            tmp <- matrix(inds[[s]]$sex, nrow=1, ncol=nrow(inds[[s]]),
                          dimnames=list(covariates=c("sex"),
                              inds=inds[[s]]$name))
            write.table(x=tmp, file=gzfile(p2f), quote=F, sep="\t", row.names=T,
                        col.names=T)
        }
        tmp <- data.frame(subgroup=names(inds),
                          file=paste0(dir, "/covariates_", names(inds),
                              ".txt.gz"))
    }
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

    if(length(nb.inds) == 1){
        if(nrow(genos.dose) > 1){
            pdf("hist_maf.pdf")
            plotHistMinAllelFreq(genos.dose, main="Tutorial of eQtlBma",
                                 breaks=50, col="grey", border="white")
            dev.off()
        }
    } else{ # different individuals between subgroups
        for(s in 1:nb.subgroups){
            if(nrow(genos.dose[,inds[[s]]$name]) > 1){
                pdf(paste0("hist_maf_", names(inds)[s], ".pdf"))
                plotHistMinAllelFreq(genos.dose[,inds[[s]]$name],
                                     main=paste0("Tutorial of eQtlBma (",
                                         names(inds)[s], ")"),
                                     breaks=50, col="grey", border="white")
                dev.off()
            }
        }
    }

    if(length(gn2sn) > 1){
        pdf("hist_cis-snps-per-gene.pdf")
        plotHistCisSnpsPerGene(gn2sn)
        dev.off()
    }
}

run <- function(params){
    set.seed(params$seed)

    inds <- simulIndividuals(params$nb.inds, params$nb.subgroups)

    gene.coords.bed <- simulGeneCoordinates(params$nb.genes,
                                            params$nb.chrs,
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
                                 params$nb.inds,
                                 inds,
                                 params$maf,
                                 params$prop.rare,
                                 params$verbose)

    explevels.genes <- simulGeneExpLevels(params$nb.inds,
                                          inds,
                                          genos.dose,
                                          gn2sn,
                                          params$pi0,
                                          params$coverr,
                                          params$nb.cores,
                                          params$verbose)

    grids <- getGrids()

    writeData(params$dir, params$nb.inds, inds, gene.coords.bed,
              snp.coords.bed, genos.dose, explevels.genes, grids, gn2sn)
}

main <- function(){
    params <- list(verbose=1,
                   pkg=NULL,
                   nb.subgroups=3,
                   nb.inds=200,
                   nb.genes=10^3,
                   nb.chrs=1,
                   avg.gene.length=10^4,
                   avg.intergenic.length=5*10^4,
                   anchor="TSS",
                   cis.radius.5p=10^3,
                   cis.radius.3p=10^3,
                   fixed.nb.snps.per.gene=NULL,
                   avg.nb.snps.per.gene=50,
                   maf=0.3,
                   prop.rare=0.0,
                   pi0=0.3,
                   coverr=1,
                   seed=1859,
                   dir=getwd(),
                   nb.cores=1)

    params <- parseCmdLine(params)

    params <- checkParams(params)

    if(params$verbose > 0){
        start.time <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
        write(paste0("START ", prog.name, " ", prog.version, " ", start.time),
              stdout())
        args <- commandArgs(trailingOnly=TRUE)
        write(paste("cmd-line:", prog.name, paste(args, collapse=" ")),
              stdout())
        write(paste0("cwd: ", getwd()), stdout())
    }

    system.time(run(params))

    if(params$verbose > 0){
        end.time <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
        difft <- as.numeric(
            difftime(as.POSIXct(end.time, format="%Y-%m-%d %H:%M:%S"),
                     as.POSIXct(start.time, format="%Y-%m-%d %H:%M:%S"),
                     units="days"))
        difft.d <- floor(difft)
        difft.h <- floor(((difft - difft.d) * 24) %% 24)
        difft.m <- floor(((difft - difft.d - difft.h/24) * 24*60) %% (24 * 60))
        difft.s <- floor(((difft - difft.d - difft.h/24 - difft.m/(24*60)) *
                          24*60*60) %% (24 * 60 * 60))
        run.length <- sprintf("%02i:%02i:%02i", difft.h, difft.m, difft.s)
        write(paste0("END ", prog.name, " ", prog.version, " ", end.time,
                     " (", run.length, ")"),
              stdout())
        ## print(object.size(x=lapply(ls(), get)), units="Kb") # return an error I don't understand
    }
}

if(! interactive())
    main()
