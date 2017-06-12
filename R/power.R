##' @importFrom stats qt
##' @importFrom stats quantile
##' @importFrom parallel mclapply
##' @importFrom matrixStats rowVars
##' @import ggplot2
NULL

##' Estimate power for microarray studies by resampling existing dataset
##'
##' Uses Tibshirani method (2006).
##'
##' Note that this function uses mclapply to parallelise sampling.  Setting `options(mc.cores=...)` can help speed things up if you have multiple cores available.  See ?mclapply for details.
##' @title tibs.power
##' @param norm matrix of normalized gene expression values from a reference expt. columns=samples, rows=genes.
##' @param ngenes.diff number of differentially expressed genes
##' @param fold.diff fold difference at differentially expressed genes
##' @param nsim number of simulations
##' @param sample.sizes number of samples in simulated experiment, can be a vector
##' @param v variance of expression for each gene, supply this if you don't have norm
##' @param p.thr p value threshold for calling significance, can be a vector to specify multiple values.  You should use only one of p.thr and n.thr.
##' @param n.thr threshold for calling the top n genes significant, can be a vector to specify multiple values.  You should use only one of p.thr and n.thr.
##' @param sample.sizes.stage2 number of pairs of samples in an optional stage 2 of the simulated experiment, can be a vector to specify multiple values.
##' @param p.thr.stage2 p value threshold for calling significance in an optional stage 2 of the simulated experiment, can NOT be a vector
##' @param p.stage2.method by default, the p.thr.stage2 will be adjusted using Bonferroni conditional on the number of genes rejected at stage one.  Set p.stage2.method="fixed" to use this p.thr.stage2 unadjusted.
##' @export
##' @return list summarising the median, 5th and 95th centiles of the false discovery rate (FDR) and false negative rate (FNR)
##' @references Tibshirani, R. (2006). A simple method for assessing sample sizes in microarray experiments. BMC Bioinformatics 7, 106.
##' @author Chris Wallace
##' @examples
##' ## Suppose we have a reference expression dataset
##' norm <- matrix(rnorm(2000),200,10)
##' ## We want to examine power for the following scenario:
##' tibs.power(norm,
##' ngenes.diff=10, # 10 genes are truly differentially expressed
##' fold.diff=2, # fold difference is 2 at those differentially expressed genes
##' nsim=10, # average 10 simulated experiments (in reality, should use many more)
##' sample.sizes=c(4,8), # first stage will use either 4 or 8 pairs
##' n.thr=50, # take the top 50 genes through to stage 2
##' sample.sizes.stage2=c(20,5), # stage 2 will use either 20 or 5 pairs
##' p.thr.stage2=0.05, # alpha at stage 2 is 0.05 ...
##' p.stage2.method="bonf") # ... and is corrected for the number of genes tested (50) by Bonferroni
#'tibs.power(norm, ngenes.diff = 50, fold.diff=2, nsim=10, sample.sizes=8, p.thr=0.05/c(200,500,1000,2000))
tibs.power <- function(norm=NULL,ngenes.diff,fold.diff,nsim,
                       sample.sizes,
                       v=NULL,
                       sample.ratio=1,
                       p.thr=numeric(0),
                       n.thr=numeric(0),
                       sample.sizes.stage2=numeric(0),
                       sample.ratio.stage2=1,
                       p.thr.stage2=numeric(0),
                       p.stage2.method=c("bonf","fixed"))  {
    p.stage2.method <- match.arg(p.stage2.method)

    if(!xor(length(p.thr),length(n.thr)))
        stop("please give either a threshold for p values, or a threshold for the top n signals")
    if(length(sample.sizes.stage2) & !length(p.thr.stage2))
        stop("please give a threshold for 2nd stage p values")
    ## if(length(sample.sizes.stage2)>1 || length(p.thr.stage2)>1)
    ##     stop("only single values for sample.sizes.stage2 and p.thr.stage2 can be given at this stage")    
    
    ## set up
    if(!is.null(norm))
        v <- apply(norm,1,var)
    ngenes <- length(v)
    vars <- list(ngenes.diff=ngenes.diff,fold.diff=fold.diff,ngroup1=sample.sizes,sample.ratio=sample.ratio)   
    if(length(p.thr)) {
        vars <- c(vars,list(alpha=p.thr))
    } else {
        vars <- c(vars,list(n=n.thr))
    }
    twostage <- FALSE
    splitvars <- c("ngenes.diff","fold.diff","ngroup1","sample.ratio")
    if(length(sample.sizes.stage2)) {
        twostage <- TRUE
        vars <- c(vars,list(npairs.2=sample.sizes.stage2))
        splitvars <- c(splitvars,"npairs.2")
    }
    ret <- do.call("expand.grid",vars)
    ret$fdr.95 <- ret$fdr.5 <- ret$fdr.med <- ret$fnr.95 <- ret$fnr.5 <- ret$fnr.med <- NA
    ret$ngroup2 <- round(ret$ngroup1 * ret$sample.ratio)
    retlist <- split(ret, ret[,splitvars])
    message("running ",nsim," simulations for each of ",length(retlist)," scenarios.")
    
    pb <- txtProgressBar(min = 0, max = length(retlist), style = 3)
                                        #done <- mclapply(seq_along(retlist), function(i) {
    i=0
    retlist <- mclapply(retlist, function(ret) {
        i=i +1
        setTxtProgressBar(pb, i); 
        ## ret <- retlist[[i]]
        ## constant within ret
        ng1 <- ret$ngroup1[1]
        ng2 <- ret$ngroup2[1]
        ngenes.diff <- ret$ngenes[1]
        fold.diff <- ret$fold.diff[1]

        fdrfnr <- vector("list",nsim)

        if(twostage)
            nstar2 <- ret$npairs.2[1]

        
        for(B in 1:nsim) {
            genes.diff <- sample(1:ngenes,ngenes.diff)

            ## stage one
            dstar <- simindep(v,ng1,ng2,genes.diff,fold.diff)
            rejected <- if(length(p.thr)) {
                            rejected.p(dstar,p.thr=ret$alpha,df=ng1+ng2-2)
                        } else {
                            rejected.n(dstar,ret$n)
                        }

            ## optional stage 2
            if(twostage) {
                ret$alpha.2 <- switch(p.stage2.method,
                                      "fixed"=p.thr.stage2,
                                      "bonf"=sapply(rejected, function(r) p.thr.stage2/sum(r)))
                d2 <- simindep(v,nstar2,genes.diff,fold.diff)
                dstar2 <- lapply(rejected, function(r) ifelse(r,d2,NA)) # only measured if rejected at stage 1
                rejected <- mapply(rejected.p, dstar2, ret$alpha.2, nstar2)
            }
            fdrfnr[[B]] <- ff(rejected, genes.diff)
        }

        fdrfnr <- do.call("rbind",fdrfnr)
        np <- ncol(fdrfnr)/2
        q <- apply(fdrfnr,2,quantile,c(0.05,0.5,0.95),na.rm=TRUE)
        ret$fdr.5 <- q[1,1:np]
        ret$fdr.med <- q[2,1:np]
        ret$fdr.95 <- q[3,1:np]
        ret$fnr.5 <- q[1,-c(1:np)]
        ret$fnr.med <- q[2,-c(1:np)]
        ret$fnr.95 <- q[3,-c(1:np)]
        return(ret)
        ## retlist[[i]] <- ret
    })
    close(pb)
    do.call("rbind",retlist)
}

simpaired <- function(norm,nstar,genes.diff,fold.diff) {
    vec <- 1:ncol(norm)
    pairs <- t(sapply(as.list(1:nstar),function(x) sample(vec,2)))
    sim.data <- norm[,pairs[,1]]-norm[,pairs[,2]]
    d.sd <- sqrt(rowVars(sim.data) / (nstar-1) )
    d <- rowMeans(sim.data)/d.sd
    d[genes.diff] <- d[genes.diff] + log2(fold.diff) / (d.sd[genes.diff])
    return(d)
}

simindep <- function(v,nstar1,nstar2=nstar1,genes.diff,fold.diff) {
    d.sd <- sqrt(v) * sqrt(1/nstar1 + 1/nstar2)
    d <- rt(length(v),df=nstar1+nstar2-2)
    d[genes.diff] <- d[genes.diff] + log2(fold.diff) / d.sd[genes.diff] 
    return(d)
}

rejected.p <- function(dstar,p.thr,df) {
    ## c.thr <- qt(p.thr/2,df,lower.tail=FALSE) 
    rejected <- lapply(p.thr, function(p) {
        thr <- qt(p/2,df,lower.tail=FALSE) 
        ifelse(is.na(dstar), FALSE, abs(dstar)>thr)
    })
    return(rejected)
}
rejected.n <- function(dstar,n.thr) {
    r <- rank(abs(dstar))
    mxr <- max(r)
    rejected <- lapply(n.thr, function(n) {
        r > mxr - n
    })
    return(rejected)
}
ff <- function(rejected,genes.diff) {
    nr <- unlist(lapply(rejected, sum, na.rm=TRUE))
    V <- unlist(lapply(rejected, function(r) sum(r[-genes.diff],na.rm=TRUE)))
    T <- unlist(lapply(rejected, function(r) sum(!r[genes.diff],na.rm=TRUE)))
    fdr <- V/nr
    fnr <- T/length(genes.diff)
    if(any(nr==0)) {
        wh <- which(nr==0)
        fdr[wh] <- 0
        fnr[wh] <- 1
    }
    c(fdr,fnr)
}

##' Plot FDR or FNR results of power.tibs()
##'
##' @title plot_fnr and plot_fdr
##' @param x object returned by power.tibs()
##' @export
##' @return plot of fnr or fdr against sample size (number of pairs) is drawn on current graphics device
##' @author Chris Wallace
plot_fdr <- function(x) {
     plotter(x,what="fdr")
}

##' @export
##' @rdname plot_fdr
plot_fnr <- function(x) {
    plotter(x,what="fnr")
}

plotter <- function(x,what=c("fdr","fnr")) {
        if("alpha" %in% colnames(x)) {
            x$lt <- as.factor(x$alpha)
            lt.nm <- "alpha"
        } else {
            x$lt <- as.factor(x$n)
            lt.nm <- "n genes"
        }

        what=match.arg(what)
        if(what=="fdr") {
            p <- ggplot(x, aes_(x=~npairs,col=~lt,y=~fdr.med,ymin=~fdr.5,ymax=~fdr.95)) + ggtitle("False Discovery Rate") + ylab("FDR")
        } else {
            p <- ggplot(x, aes_(x=~npairs,col=~lt,y=~fnr.med,ymin=~fnr.5,ymax=~fnr.95)) + ggtitle("False Negative Rate") + ylab("FNR")
        }
 
        p + geom_pointrange() + geom_path(aes_(lty=~lt)) + xlab("Sample Size") + facet_wrap(ngenes.diff ~ fold.diff, labeller=label_both) + scale_colour_discrete(lt.nm) + scale_linetype_discrete(lt.nm)
}
    
    
##' Estimate power for microarray studies by resampling existing dataset
##'
##' Uses Tibshirani method (2006).
##'
##' Note that this function uses mclapply to parallelise sampling.  Setting `options(mc.cores=...)` can help speed things up if you have multiple cores available.  See ?mclapply for details.
##' @title tibs.power
##' @param norm matrix of normalized gene expression values from a reference expt. columns=samples, rows=genes.
##' @param ngenes.diff number of differentially expressed genes
##' @param fold.diff fold difference at differentially expressed genes
##' @param nsim number of simulations
##' @param sample.sizes number of pairs of samples in simulated experiment, can be a vector
##' @param p.thr p value threshold for calling significance, can be a vector to specify multiple values.  You should use only one of p.thr and n.thr.
##' @param n.thr threshold for calling the top n genes significant, can be a vector to specify multiple values.  You should use only one of p.thr and n.thr.
##' @param sample.sizes.stage2 number of pairs of samples in an optional stage 2 of the simulated experiment, can be a vector to specify multiple values.
##' @param p.thr.stage2 p value threshold for calling significance in an optional stage 2 of the simulated experiment, can NOT be a vector
##' @param p.stage2.method by default, the p.thr.stage2 will be adjusted using Bonferroni conditional on the number of genes rejected at stage one.  Set p.stage2.method="fixed" to use this p.thr.stage2 unadjusted.
##' @export
##' @return list summarising the median, 5th and 95th centiles of the false discovery rate (FDR) and false negative rate (FNR)
##' @references Tibshirani, R. (2006). A simple method for assessing sample sizes in microarray experiments. BMC Bioinformatics 7, 106.
##' @author Chris Wallace
##' @examples
##' ## Suppose we have a reference expression dataset
##' norm <- matrix(rnorm(2000),200,10)
##' ## We want to examine power for the following scenario:
##' tibs.power(norm,
##' ngenes.diff=10, # 10 genes are truly differentially expressed
##' fold.diff=2, # fold difference is 2 at those differentially expressed genes
##' nsim=10, # average 10 simulated experiments (in reality, should use many more)
##' sample.sizes=c(4,8), # first stage will use either 4 or 8 pairs
##' n.thr=50, # take the top 50 genes through to stage 2
##' sample.sizes.stage2=c(20,5), # stage 2 will use either 20 or 5 pairs
##' p.thr.stage2=0.05, # alpha at stage 2 is 0.05 ...
##' p.stage2.method="bonf") # ... and is corrected for the number of genes tested (50) by Bonferroni

tibs.paired <- function(norm,ngenes.diff,fold.diff,nsim,
                        sample.sizes,
                        sample.ratio=1,
                       p.thr=numeric(0),
                       n.thr=numeric(0),
                       sample.sizes.stage2=numeric(0),
                        sample.ratio2=1,
                       p.thr.stage2=numeric(0),
                       p.stage2.method=c("bonf","fixed"))  {
    p.stage2.method <- match.arg(p.stage2.method)

    if(!xor(length(p.thr),length(n.thr)))
        stop("please give either a threshold for p values, or a threshold for the top n signals")
    if(length(sample.sizes.stage2) & !length(p.thr.stage2))
        stop("please give a threshold for 2nd stage p values")
    ## if(length(sample.sizes.stage2)>1 || length(p.thr.stage2)>1)
    ##     stop("only single values for sample.sizes.stage2 and p.thr.stage2 can be given at this stage")    
    
    ## set up
    ngenes <- nrow(norm)
    vars <- list(ngenes.diff=ngenes.diff,fold.diff=fold.diff,npairs=sample.sizes)    
    if(length(p.thr)) {
        vars <- c(vars,list(alpha=p.thr))
    } else {
        vars <- c(vars,list(n=n.thr))
    }
    twostage <- FALSE
    splitvars <- c("ngenes.diff","fold.diff","npairs")
    if(length(sample.sizes.stage2)) {
        twostage <- TRUE
        vars <- c(vars,list(npairs.2=sample.sizes.stage2))
        splitvars <- c(splitvars,"npairs.2")
    }
    ret <- do.call("expand.grid",vars)
    ret$fdr.95 <- ret$fdr.5 <- ret$fdr.med <- ret$fnr.95 <- ret$fnr.5 <- ret$fnr.med <- NA
    retlist <- split(ret, ret[,splitvars])
    message("running ",nsim," simulations for each of ",length(retlist)," scenarios.")
    
    retlist <- mclapply(retlist, function(ret) {

        ## constant within ret
        nstar <- ret$npairs[1]
        ngenes.diff <- ret$ngenes[1]
        fold.diff <- ret$fold.diff[1]

        fdrfnr <- vector("list",nsim)

        if(twostage)
            nstar2 <- ret$npairs.2[1]
        
        for(B in 1:nsim) {
            genes.diff <- sample(1:ngenes,ngenes.diff)

            ## stage one
            ## system.time(dstar <- simd(norm,nstar,genes.diff,fold.diff))
            ## system.time(dstar <- simd2(norm,nstar,genes.diff,fold.diff))
            ## system.time(dstar <- simd3(norm,nstar,genes.diff,fold.diff))
            dstar <- simd(norm,nstar,genes.diff,fold.diff)
            rejected <- if(length(p.thr)) {
                            rejected.p(dstar,ret$alpha,nstar)
                        } else {
                            rejected.n(dstar,ret$n)
                        }

            ## optional stage 2
            if(twostage) {
                ret$alpha.2 <- switch(p.stage2.method,
                                      "fixed"=p.thr.stage2,
                                      "bonf"=sapply(rejected, function(r) p.thr.stage2/sum(r)))
                d2 <- simd(norm,nstar2,nstar2*sample.ratio2,genes.diff,fold.diff)
                dstar2 <- lapply(rejected, function(r) ifelse(r,d2,NA)) # only measured if rejected at stage 1
                rejected <- mapply(rejected.p, dstar2, ret$alpha.2, nstar2)
            }
            fdrfnr[[B]] <- ff(rejected, genes.diff)
        }
        fdrfnr <- do.call("rbind",fdrfnr)
        np <- ncol(fdrfnr)/2
        q <- apply(fdrfnr,2,quantile,c(0.05,0.5,0.95),na.rm=TRUE)
        ret$fdr.5 <- q[1,1:np]
        ret$fdr.med <- q[2,1:np]
        ret$fdr.95 <- q[3,1:np]
        ret$fnr.5 <- q[1,-c(1:np)]
        ret$fnr.med <- q[2,-c(1:np)]
        ret$fnr.95 <- q[3,-c(1:np)]
        return(ret)
    })
    do.call("rbind",retlist)
}


