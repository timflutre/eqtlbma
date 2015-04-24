## `utils_eqtlbma.R' contains utility functions for the eQtlBma package
## Copyright (C) 2013-2015 Timothée Flutre
## License: GPL-3+
## Persons: Timothée Flutre [cre,aut]

utils_eqtlbma.version <- "1.3.1a" # same as in configure.ac

##' Transform gene expression levels into a N(0,1) via quantile normalization
##'
##'
##' Performed gene by gene, across samples
##' @title Quantile normalization
##' @param mat matrix with genes in rows and samples in columns
##' @param break.ties.rand break ties randomly (default=TRUE)
##' @param seed to make it reproducible
##' @return Matrix with transformed genes in rows and samples in columns
transformGeneExpInStdNormal <- function(mat, break.ties.rand=TRUE,
                                        seed=1859){
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

##' Remove confounders from a matrix of gene expression levels
##'
##' Can be PCs, PEER factors, etc;  linear regression per gene
##' @title Remove "expression" confounders
##' @param X matrix with samples in rows and genes in columns,
##'  will be centered and scaled before confounders are removed
##' @param confounders matrix with samples in rows and confounders in columns
##' @return Matrix of residuals
removeConfoundersFromGeneExp <- function(X, confounders){
    stopifnot(is.matrix(X),
              is.matrix(confounders),
              nrow(X) == nrow(confounders))
    if(nrow(X) > ncol(X))
        warning("input matrix doesn't seem to have samples in rows and genes in columns")

    res <- lm.fit(x=confounders, y=scale(X, center=TRUE, scale=TRUE))

    return(t(res$residuals))
}

##' Make the grid used to compute Bayes Factors
##'
##' @title Grid
##' @param grid.type "general" indicates the meta-analysis grid (large),
##' otherwise the configuration grid (small) is returned
##' @param no.het if TRUE, the grid is built without heterogeneity
##' @return Matrix with the grid
makeGrid <- function(grid.type="general", no.het=FALSE){
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

##' Plot the histogram of the minor allele frequency per SNP
##'
##' Missing values (encoded as NA) are discarded.
##' @title Histogram of minor allele frequencies
##' @param genos.dose matrix or data.frame of allele dose ([0;1]) with
##' SNPs in rows and samples in columns
##' @param maf vector of minor allele frequencies
##' @param main string for the main title (default="")
##' @param xlim default=c(0,0.5)
##' @return vector of minor allele frequencies (to avoid recomputing it)
plotHistMinAllelFreq <- function(genos.dose=NULL, maf=NULL, main="",
                                 xlim=c(0,0.5), ...){
    stopifnot(! is.null(genos.dose) || ! is.null(maf))

    if(! is.null(genos.dose) & is.null(maf)){
        if(nrow(genos.dose) < ncol(genos.dose))
            warning("input matrix doesn't seem to have SNPs in rows and samples in columns")
        genos.dose <- as.matrix(genos.dose)
        maf <- apply(genos.dose, 1, function(x){
            x <- x[complete.cases(x)] # discard missing values (encoded as NA)
            tmp <- sum(x) / (2 * length(x))
            ifelse(tmp <= 0.5, tmp, 1 - tmp)
        })
    }

    tmp <- hist(x=maf, xlab="Minor allele frequency", ylab="Number of SNPs",
                main=main, xlim=xlim, ...)

    invisible(return(maf))
}

##' Identify SNPs in cis of each gene
##'
##' Use findOverlaps() from GenomicRanges, TSS +- cis.radius
##' @title SNPs in cis of genes
##' @param gene.coords.bed data.frame in BED format
##' @param snp.coords.bed data.frame in BED format
##' @param anchor gene boundary(ies) for the cis region (TSS or TSS+TES)
##' @param cis.radius.5p length of the 5' half of the cis region (in bp)
##' @param cis.radius.3p length of the 3' half of the cis region (in bp)
##' @param verbose verbosity level (0/default=1/2)
##' @return List with gene names as components, each having a corresponding
##' vector of SNP names (reported genes have at least one cis SNP)
findCisSnpsPerGene <- function(gene.coords.bed, snp.coords.bed, anchor,
                               cis.radius.5p, cis.radius.3p, verbose=1){
    library("GenomicRanges") # from Bioconductor
    stopifnot("chr" %in% colnames(gene.coords.bed),
              "start" %in% colnames(gene.coords.bed),
              "end" %in% colnames(gene.coords.bed),
              "strand" %in% colnames(gene.coords.bed),
              "chr" %in% colnames(snp.coords.bed),
              "start" %in% colnames(snp.coords.bed),
              "end" %in% colnames(snp.coords.bed),
              anchor %in% c("TSS", "TSS+TES"),
              is.numeric(cis.radius.5p),
              is.numeric(cis.radius.3p))

    ## convert SNP coordinates from data.frame to GRanges structure
    snp.coords.gr <-
        GRanges(seqnames=Rle(snp.coords.bed$chr),
                ranges=IRanges(start=as.numeric(snp.coords.bed$start) + 1,
                    end=as.numeric(snp.coords.bed$end)))
    names(snp.coords.gr) <- snp.coords.bed$name

    ## idem for genes, but now defined as cis regions
    cis.regions.gr <-
        GRanges(seqnames=Rle(gene.coords.bed$chr),
                ranges=IRanges(start=as.numeric(gene.coords.bed$start) + 1,
                    end=as.numeric(gene.coords.bed$end)))
    names(cis.regions.gr) <- gene.coords.bed$name
    strand(cis.regions.gr) <- gene.coords.bed$strand

    idx.strand.p <- which(strand(cis.regions.gr) == "+")
    idx.strand.m <- which(strand(cis.regions.gr) == "-")
    if(anchor == "TSS"){
        end(cis.regions.gr[idx.strand.p]) <- start(cis.regions.gr[idx.strand.p])
        start(cis.regions.gr[idx.strand.m]) <- end(cis.regions.gr[idx.strand.m])
    }

    start(cis.regions.gr[idx.strand.p]) <-
        ifelse(start(cis.regions.gr[idx.strand.p]) - cis.radius.5p > 0,
               start(cis.regions.gr[idx.strand.p]) - cis.radius.5p, 1)
    end(cis.regions.gr[idx.strand.p]) <- end(cis.regions.gr[idx.strand.p]) +
        cis.radius.3p

    end(cis.regions.gr[idx.strand.m]) <- end(cis.regions.gr[idx.strand.m]) + cis.radius.5p
    start(cis.regions.gr[idx.strand.m]) <-
        ifelse(start(cis.regions.gr[idx.strand.m]) - cis.radius.3p > 0,
               start(cis.regions.gr[idx.strand.m]) - cis.radius.3p, 1)

    idx.gs.pairs <- as.matrix(findOverlaps(cis.regions.gr, snp.coords.gr))

    gn2sn <- sapply(split(idx.gs.pairs[,2],
                          names(cis.regions.gr)[idx.gs.pairs[,1]]),
                    function(i) names(snp.coords.gr)[i])

    if(verbose > 0){
        message("nb of gene(s) with at least one cis SNP: ", length(gn2sn))
        message("nb of gene-SNP pairs: ", sum(sapply(gn2sn, length)))
    }

    return(gn2sn)
}

##' Plot the histogram of the number of cis SNPs per gene
##'
##' @title Histogram of cis SNPs
##' @param gn2sn list with one component per gene corresponding to a vector
##' of SNP names in cis of that gene; see findCisSnpsPerGene()
plotHistCisSnpsPerGene <- function(gn2sn){
    stopifnot(is.list(gn2sn))

    nb.cis.snps.per.gene <- sapply(gn2sn, length)
    tmp <- hist(x=nb.cis.snps.per.gene,
                main=paste0(length(gn2sn), " genes and ",
                    length(unique(do.call(c, gn2sn))), " SNPs"),
                xlab="number of cis SNPs per gene",
                col="grey", border="white")

    x <- tmp$breaks[ceiling((3/4)*length(tmp$breaks))]
    y <- max(tmp$counts) / 2
    text(x=x, y=1.2*y, labels=paste0("minimum: ", min(nb.cis.snps.per.gene), " SNPs"))
    text(x=x, y=y, labels=paste0("median: ", median(nb.cis.snps.per.gene), " SNPs"))
    text(x=x, y=0.8*y, labels=paste0("maximum: ", max(nb.cis.snps.per.gene), " SNPs"))
}

##' Get all binary configurations
##'
##' Done recursively
##' @title Configurations
##' @param nb.subgroups
##' @param verbose
##' @return Matrix with configurations in rows and subgroups in columns
getBinaryConfigs <- function(nb.subgroups=1, verbose=0){
    stopifnot(is.numeric(nb.subgroups), length(nb.subgroups) == 1)
    if(nb.subgroups >= 20)
        warning("nb of subgroups may be too high, better to keep it below 20!",
                call.=FALSE, immediate.=TRUE)
    if(verbose > 0)
        cat("list all configurations for the effect of a SNP in each subgroup\n")
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

##' Simulate from a matrix-variate Normal distribution
##'
##' @param n the number of samples required (default=1)
##' @param M a matrix giving the means of the variables
##' @param U a positive-definite symmetric matrix specifying the among-row
##' covariance matrix of the variables
##' @param V a positive-definite symmetric matrix specifying the among-column
##' covariance matrix of the variables
##' @return Array with samples in the third dimension
matvrnorm <- function(n=1, M, U, V){
    library(MASS)
    stopifnot(nrow(M) == nrow(U),
              ncol(M) == nrow(V),
              nrow(U) == ncol(U),
              nrow(V) == ncol(V))

    tmp <- lapply(1:n, function(i){
        matrix(data=mvrnorm(n=1, mu=c(M), Sigma=V %x% U),
               nrow=nrow(M), ncol=ncol(M))
    })

    return(array(data=do.call(c, tmp), dim=c(nrow(M), ncol(M), n)))
}

##' Calculate the probability to be active in s subgroups
##' by summing configuration probabilities
##'
##' @title Activity probabilities
##' @param configs data.frame with column "id" and column "proba"
##' ("id" should be "1" or "1-2" or "1-3-4", etc)
##' @param plot.it plot the probas (default=FALSE)
##' @param main main title for the plot
##' @return List with probas per subgroup
calcActivityProbasPerSubgroup <- function(configs, plot.it=FALSE,
                                          main=NULL){
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

##' Calculate the pairwise eQTL sharing by marginalizing other subgroups
##'
##' @title Pairwise eQTL sharing
##' @param configs data.frame with configuration probabilities
##' @param reformat if TRUE, configs should have a column "id" ("1", "1-3-4", etc)
##' and column "proba"; else, configs should have one column
##' per subgroup (with 0 or 1 per row) and a column  "proba"
##' @param subgroups vector of names (otherwise "subgroup.s" will be used)
##' @return Matrix with element i,j being
##' Pr(eQTL in subgroup j | eQTL in subgroup i)
calcMarginalPairwiseEqtlSharing <- function(configs, reformat=TRUE,
                                            subgroups=NULL){
    stopifnot(is.data.frame(configs),
              is.logical(reformat))

    S <- log2(nrow(configs) + 1)
    if(is.null(subgroups))
        subgroups <- paste0("subgroup.", 1:S)

    ## reformat the data
    if(reformat){
        stopifnot("id" %in% colnames(configs),
                  "proba" %in% colnames(configs))
        tmp <- list()
        for(s in subgroups)
            tmp[[s]] <- rep(0, nrow(configs))
        tmp[["proba"]] <- configs$proba
        tmp <- do.call(cbind, tmp)
        for(i in 1:nrow(configs))
            tmp[i, as.numeric(strsplit(configs$id[i], "-")[[1]])] <- 1
    } else
        tmp <- configs

    ## compute the marginals
    pi1.marginal <- matrix(nrow=S, ncol=S,
                           dimnames=list(subgroups, subgroups))
    diag(pi1.marginal) <- 1
    for(j in 1:S){ # for each column
        for(i in 1:S){ # for each row
            num <- which(tmp[,i] == 1 & tmp[,j] == 1)
            denom <- which(tmp[,i] == 1)
            pi1.marginal[i,j] <- sum(tmp[num,"proba"]) / sum(tmp[denom,"proba"])
        }
    }

    return(pi1.marginal)
}

##' After the EM has been fitted on each pair of subgroups, this function
##' can be used to reformat the results into a pairwise "eQTL sharing" matrix
##'
##' @param dat data.frame with 4 columns subgroup1|subgroup2<sep>proba<sep>subgroup2|subgroup1<sep>proba
##' and S(S-1)/2 rows where S is the nb of subgroups
##' @return Matrix so that mat[i,j] corresponds to Pr(active in subgroup j | active in subgroup i)
reformatPairwiseEqtlSharing <- function(dat){
    stopifnot(ncol(dat) == 4)

    S <- (1 + sqrt(1 + 4 * 2 * nrow(dat))) / 2
    message(paste0("nb of subgroups: ", S))

    subgroups <- sort(unique(do.call(c, lapply(dat[,1], function(x){
        strsplit(x, "\\|")[[1]]
    }))))
    message(paste(subgroups, collapse="-"))

    mat <- matrix(nrow=S, ncol=S,
                  dimnames=list(subgroups, subgroups))

    diag(mat) <- 1
    for(k in 1:nrow(dat)){
        j <- strsplit(dat[k,1], "\\|")[[1]][1]
        i <- strsplit(dat[k,1], "\\|")[[1]][2]
        mat[i,j] <- as.numeric(dat[k,2])
        mat[j,i] <- as.numeric(dat[k,4])
    }

    return(mat)
}

##' Estimate pi0 (proba for a null hypothesis to be true)
##' via the EBF procedure
##'
##' Wen, arXiv:1311.3981
##' @title EBF procedure
##' @param log10.bfs vector containing the log10(BF) of each test
##' @param verbose verbosity level (0/default=1)
##' @return Numeric
estimatePi0WithEbf <- function(log10.bfs, verbose=1){
    stopifnot(is.numeric(log10.bfs), is.vector(log10.bfs))
    if(verbose > 0)
        message(paste0("nb of tests: ", length(log10.bfs)))

    tmp <- log10.bfs[order(log10.bfs)] # sort in increasing order

    d0 <- which(cumsum(10^tmp) / seq_along(10^tmp) >= 1)[1]
    if(verbose > 0)
        message(paste0("cutoff at the ", d0, "-th BF"))

    pi0.ebf <- d0 / length(tmp)
    if(verbose > 0)
        message(paste0("estimate pi0-hat = ",
                       format(x=pi0.ebf, scientific=TRUE, digits=6)))

    return(pi0.ebf)
}

##' Estimate pi0 (proba for a null hypothesis to be true)
## via the QBF procedure
##'
##' Wen, arXiv:1311.3981
##' @title QBF procedure
##' @param log10.bfs matrix with tests in rows and two columns, the true log10(BF)
##' and the gamma-quantile log10(BF) under the null
##' @param gamma level of the quantile (e.g. 0.5 for the median)
##' @param verbose verbosity level (0/default=1)
##' @return Numeric
estimatePi0WithQbf <- function(log10.bfs, gamma, verbose=1){
    stopifnot(is.numeric(log10.bfs), is.matrix(log10.bfs),
              ncol(log10.bfs) == 2)
    if(verbose > 0)
        message(paste0("nb of tests: ", nrow(log10.bfs)))

    pi0.qbf <- sum(log10.bfs[,1] <= log10.bfs[,2]) / (nrow(log10.bfs) * gamma)
    if(verbose > 0)
        message(paste0("estimate pi0-hat = ",
                       format(x=pi0.qbf, scientific=TRUE, digits=6)))

    return(pi0.qbf)
}

##' Call significant tests by controlling the Bayesian FDR
##'
##' procedure from Newton et al (Biostatistics, 2004) also described
##' in Muller et al (JASA, 2006)
##' @title Bayesian FDR control
##' @param log10.bfs vector containing the log10(BF) of each test
##' @param pi0 estimate of the proba for a null hypothesis to be true
##' @param fdr.level threshold below which a null is rejected (default=0.05)
##' @param verbose verbosity level (0/default=1)
##' @return Vector of booleans, TRUE if null is rejected (thus called significant)
controlBayesFdr <- function(log10.bfs, pi0, fdr.level=0.05, verbose=1){
    stopifnot(is.numeric(log10.bfs), is.vector(log10.bfs))
    if(verbose > 0)
        message(paste0(length(log10.bfs), " tests and pi0-hat = ",
                       format(x=pi0, scientific=TRUE, digits=6)))

    ## compute the posterior probability of each test being null
    post.null <- pi0 / (pi0 + (1-pi0) * 10^log10.bfs)

    ## find the cutoff for which their cumulative mean is >= fdr.level
    idx <- order(post.null)
    post.null <- post.null[idx]
    for(L in 1:length(post.null))
        if(mean(post.null[1:L]) >= fdr.level)
            break
    log10.bf.L <- log10.bfs[idx[L]]
    if(verbose > 0)
        message(paste0(L, " significant tests, at cutoff log10(BF)=", log10.bf.L))

    significant <- log10.bfs >= log10.bf.L

    return(significant)
}

##' Compute log_{10}(\sum_i w_i 10^x_i) stably
##'
##' @param x vector of log10(values)
##' @param weights optional vector of weights (equal weights if not specified)
##' @return Numeric
log10WeightedSum <- function(x, weights=NULL){
    stopifnot(is.numeric(x), is.vector(x))

    if(! is.null(weights)){
        stopifnot(is.vector(weights))
    } else
        weights <- rep(1/length(x), length(x))

    max <- max(x)
    max + log10(sum(weights * 10^(x - max)))
}

##' Calculate the asymptotic Bayes Factor proposed by Wakefield
##' in Genetic Epidemiology 33:79-86 (2009)
##'
##' http://dx.doi.org/10.1002/gepi.20359
##' @title Asymptotic Bayes factor from Wakefield
##' @param theta.hat MLE of the additive genetic effect
##' @param V variance of theta.hat
##' @param W variance of the prior on theta
##' @param log10 to return the log10 of the ABF (default=TRUE)
##' @return Numeric
calcAsymptoticBayesFactorWakefield <- function(theta.hat, V, W, log10=TRUE){
    z2 <- theta.hat^2 / V # Wald statistic

    log10.ABF <- 0.5 * log10(V) - 0.5 * log10(V + W) +
        (0.5 * z2 * W / (V + W)) / log(10)

    if(log10)
        return(log10.ABF)
    else
        return(10^log10.ABF)
}

##' Calculate the exact Bayes Factor proposed by Servin and Stephens
##' in PLoS Genetics 3,7 (2007)
##'
##' http://dx.doi.org/10.1371/journal.pgen.0030114
##' @title Bayes factor from Servin and Stephens
##' @param G vector of genotypes
##' @param Y vector of phenotypes
##' @param sigma.a variance of the prior on the additive genetic effect
##' @param sigma.d variance of the prior on the dominance genetic effect
##' @param log10 to return the log10 of the ABF (default=TRUE)
##' @return Numeric
calcExactBayesFactorServinStephens <- function(G, Y, sigma.a, sigma.d, log10=TRUE){
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
