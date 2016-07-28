
##' Estimate power for microarray studies by resampling existing dataset
##'
##' Uses Tibshirani method (2006)
##' @title tibs.power
##' @param norm matrix of normalized gene expression values from a reference expt
##' @param ngenes.diff number of differentially expressed genes
##' @param fold.diff fold difference at differentially expressed genes
##' @param nsim number of simulations
##' @param sample.sizes number of pairs of samples in simulated experiment, can be a vector
##' @param p.thr p value threshold for calling significance, can be a vector
##' @export
##' @return list summarising the median, 5th and 95th centiles of the false discovery rate (FDR) and false negative rate (FNR)
##' @references Tibshirani, R. (2006). A simple method for assessing sample sizes in microarray experiments. BMC Bioinformatics 7, 106.
##' @author Chris Wallace
tibs.power <- function(norm,ngenes.diff,fold.diff,nsim,
                       sample.sizes,
                       p.thr) { #c(seq(4,10),seq(12,20,2),seq(24,40,4)) # n pairs
                                        #p.thr <- c(1e-4,1e-5,1e-6,1e-8)
    fdr.5 <- fdr.95 <- fdr.med <- fnr.5 <- fnr.95 <- fnr.med <- matrix(NA,length(sample.sizes),length(p.thr))
    i <- j <- 0
    ngenes <- nrow(norm)
    for(nstar in sample.sizes) {
        
        cat("\nsample size",nstar,"\n")
        ## pick nstar "pairs" from norm
        pairs <- t(sapply(as.list(1:nstar),function(x) sample(1:12,2)))
        sim.data <- norm[,pairs[,1]]-norm[,pairs[,2]]
        d.sd <- apply(sim.data,1,sd)
        d <- apply(sim.data,1,mean)/d.sd
        
        i <- i +1
        fdr <- fnr <- matrix(NA,nsim,length(p.thr)) # i=sim, j=p
        for(B in 1:nsim) {
            genes.diff <- sample(1:ngenes,ngenes.diff)
            dstar <- d
            dstar[genes.diff] <- d[genes.diff] + log2(fold.diff) / (d.sd[genes.diff]/sqrt(nstar))
            
            j <- 0
            for(p in p.thr) {
                
                cat("\tp=",p)
                j <- j+1
                c.thr <- qt(p/2,nstar-1,lower.tail=FALSE)    
                
                rejected <- abs(dstar)>c.thr
                nr <- sum(rejected)
                if(nr==0) {
                    fdr[B,j] <- NA
                    fnr[B,j] <- NA
                    next
                }
                
                V <- sum(rejected[-genes.diff])
                T <- sum(!rejected[genes.diff])
                fdr[B,j] <- V/nr
                fnr[B,j] <- T/ngenes.diff
            }
        }
        
        q <- apply(fdr,2,quantile,c(0.05,0.5,0.95),na.rm=TRUE)
        fdr.5[i,] <- q[1,]
        fdr.med[i,] <- q[2,]
        fdr.95[i,] <- q[3,]
        q <- apply(fnr,2,quantile,c(0.05,0.5,0.95),na.rm=TRUE)
        fnr.5[i,] <- q[1,]
        fnr.med[i,] <- q[2,]
        fnr.95[i,] <- q[3,]
    }
    return(list(fdr.5=fdr.5,fdr.med=fdr.med,fdr.95=fdr.95,fnr.5=fnr.5,fnr.med=fnr.med,fnr.95=fnr.95))
}
    
##' Plot FDR or FNR results of power.tibs()
##'
##' @title tibs.plot
##' @param x vector giving position of points on x axis
##' @param y matrix, each column to be plotted against x
##' @param low optional matrix of lower interval about y
##' @param high optional matrix of upper interval about y
##' @param ... other params passed to plot()
##' @export
##' @return plot is drawn on current graphics device
##' @author Chris Wallace
    tibs.plot <- function(x,y,low=NULL,high=NULL,...) { # x a vector, y is a matrix, each column to be plotted against x
        R <- range(y,na.rm=TRUE)
        if(!is.null(low))
            R <- range(c(low,R),na.rm=TRUE)
        if(!is.null(high))
            R <- range(c(high,R),na.rm=TRUE)
        
        plot(NA,xlim=range(x),ylim=R,...)
        if(!is.null(low))
            for(j in 1:ncol(y))
                lines(x,low[,j],col=j,lty=3)
        if(!is.null(high))
            for(j in 1:ncol(y))
                lines(x,high[,j],col=j,lty=3)
        for(j in 1:ncol(y))
            lines(x,y[,j],col=j,pch=j,type="b")
        
    }
    
    
