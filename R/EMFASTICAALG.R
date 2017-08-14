###########################
## EMFASTICAALG.R
##
## Written by Xiaotian Zhu
###########################

##
## accessory functions
##

save_rng <- function(savefile="tempfile1"){# function for saving radom seed info
  oldseed <- oldRNGkind <- NULL
  rm(oldseed)
  rm(oldRNGkind)
  if (exists(".Random.seed"))  {
    oldseed <- get(".Random.seed", .GlobalEnv)
  } else stop("don't know how to save before set.seed() or r*** call")
  oldRNGkind <- RNGkind()
  save("oldseed","oldRNGkind",file=savefile)
  invisible(savefile)
}

restore_rng <- function(savefile) {# function for loading radom seed info
  oldseed <- oldRNGkind <- NULL
  rm(oldseed)
  rm(oldRNGkind)
  load(savefile)
  do.call("RNGkind",as.list(oldRNGkind))  
  assign(".Random.seed", oldseed, .GlobalEnv)
}



#' ESTIMATEDMEMBER
#' 
#' A function calculates estimated class membership from an EMFASTICAALG
#' object.
#' 
#' 
#' @param rstICAMIX An EMFASTICAALG object.
#' @return %% ~Describe the value returned %% If it is a LIST, use
#' \item{estimatedmember }{A factor with levels representing the estimated
#' classes.} %% \item{comp2 }{Description of 'comp2'} %% ...
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' ## An Example that runs the NSMM-ICA algorithm on Cohen's tone data
#' data(tonedata, package="mixtools")
#' 
#' b <- EMFASTICAALG(tonedata, 2, h=0, tol=1e-8)
#' estimatedgroup <- ESTIMATEDMEMBER(b) # estimated species info
#' 
#' 
#' @export ESTIMATEDMEMBER
ESTIMATEDMEMBER <-
  function(rstICAMIX){# function calculates estimated class membership from an EMFASTICAALG object
    postmembership <- rstICAMIX$MembershipProbs == apply(rstICAMIX$MembershipProbs, 1, max) # estimated membership matrix
    ncluster <- ncol(rstICAMIX$MembershipProbs)
    estimatedmember <- factor(rep((1:ncluster),ncol(rstICAMIX$InputData))[as.vector(t(postmembership))]) # estimated membership info
    estimatedmember
  }



#' CLASSDIFFRATE
#' 
#' A function calculates classification difference rate between two factors. It
#' is used in interpreting info stored in EMFASTICA object.
#' 
#' 
#' @param factor1 First factor.
#' @param factor2 Second factor of the same length as the First factor.
#' @return %% ~Describe the value returned %% If it is a LIST, use \item{answer
#' }{The percentage of instances when factor1[i] is not equal factor2[i].} %%
#' \item{comp2 }{Description of 'comp2'} %% ...
#' @keywords ~kwd1 ~kwd2
#' @examples
#' ## An example evaluates the classification difference rate
#' ## between two classification results in the form of factors
#' fac1<-factor(c(1,4,2,3,1,1,3,3,1,2,2,1))
#' fac2<-factor(c(3,1,2,2,1,2,4,3,2,3,1,1))
#' CLASSDIFFRATE(fac1, fac2)
#' 
#' @export CLASSDIFFRATE
CLASSDIFFRATE <-
  function(factor1, factor2){# function calculates classification difference rate between two factors
    length(which(as.numeric(levels(factor1))[factor1] != as.numeric(levels(factor2))[factor2]))/length(factor1)
  }



#' ATRANSDENSITY
#' 
#' A function evaluates density of linearly transformed random vector on a
#' given grid. It is used in processing EMFASTICAALG object to obtain density
#' estimation of the mixture components.
#' 
#' 
#' @param grid A matrix whose columns store the grid points.
#' @param A Matrix for the linear transformation.
#' @param f Density function before the linear transformation.
#' @return %% ~Describe the value returned %% If it is a LIST, use
#' \item{answer}{Matrix of the same size as grid, with each element being the
#' evaluated linear transformed density at the corresponding grid point.} %%
#' \item{comp2 }{Description of 'comp2'} %% ...
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' ## An example that evaluates the 2-D uniform distribution on a linear transformation of [1,3]x[1,3]
#' ## f1ind is the density of the uniform distribution on [1,3]^r
#' f1ind <- function(grid){# mixture component 1 original signal density function
#'   n <- ncol(grid)
#'   r <- nrow(grid)
#'   answer <- rep(1,n)
#'   for(i in 1:n){
#'     for(j in 1:r){
#'       answer[i] <- answer[i] * (grid[j,i] >= 1 & grid[j,i] <= 3) / 2
#'     }
#'   }
#'   answer
#' }
#' 
#' A <- matrix(c(6, 9, -12, 15), 2, 2, byrow = FALSE)
#' 
#' gridpoints <- t(as.matrix(expand.grid(seq(-32,12,2),seq(18,80,2))))
#' 
#' f1trans <- ATRANSDENSITY(gridpoints, A, f1ind)
#' 
#' plot(t(gridpoints),col=(f1trans>0))
#' 
#' @export ATRANSDENSITY
ATRANSDENSITY <- function(grid, A, f){
  # function evaluates density of linearly transformed random vector
  # on a given grid of which the columns store the grid points
  n <- ncol(grid)
  grid_ <- solve(A) %*% grid
  answer <- f(grid_)
  answer <- answer / abs(det(A))
  answer
}


##
## main function
##
## an R wrapper for carrying out NSMM-ICA
## on nonparametric multivariate ICA mixture data
## 
##


#' EMFASTICAALG
#' 
#' An R wrapper for carrying out NSMM-ICA on nonparametric multivariate ICA
#' mixture data.
#' 
#' 
#' @param DataMatrix A matrix of which the rows are data entries. Its dimension
#' is \code{n} by \code{r}.
#' @param numCluster Predetermined number of mixing components \code{m}.
#' @param h Bandwidth. If \code{h} is set equal zero (default), iterative
#' bandwidth selection will be used.
#' @param maxiter Maximum number of iterations. Default is \code{300}.
#' @param icaiter Maximum number of ICA iterations in each step. Default is
#' \code{150}.
#' @param tol Threshold that defines convergence (of the outer loop). Default
#' is \code{1}\code{e}-\code{6}.
#' @param verb \code{TRUE} (default) or \code{FALSE} indicating whether to
#' print out info at each iteration.
#' @param combine \code{TRUE} (default) or \code{FALSE} indicating whether to
#' implement the ICA step.
#' @param seednum Seed number (default is \code{82196}) used in kmeans before
#' 1st iteration.
#' 
#' @return %% ~Describe the value returned %% If it is a LIST, use 
#' The returned value is an \code{EMFASTICAALG} object which consists of a list of items:
#' \item{$InputData }{A matrix of which the columns are data entries. Its
#' dimension is \code{r} by \code{n}.} \item{$Lambdas }{A matrix where rows
#' store estimated mixing weights from each iteration.} \item{$WMtrs }{List of
#' \code{r} by \code{r} unmixing matrices for each of the m clusters.}
#' \item{$WUnmixZ }{List of unmixing matrices for whitened data for each of the
#' m clusters.} \item{$OriginalSignals }{List of Recovered ICA components for
#' each of the \code{m} clusters.} \item{$ProductDensity }{\code{m} by \code{n}
#' matrix where each row stores the estimated density value of the observed
#' data points for each of the m clusters.} \item{$MembershipProbs }{\code{n}
#' by \code{m} matrix where each row stores the component membership
#' probabilities of the corresponding data point.} \item{$ObjValue }{Vector
#' holding values of data loglikelihood.} \item{$ICABandWidth }{Matrix holding
#' choices of bandwidth for original signals.} \item{$call }{The function call
#' that results in the returned object.} \item{$time }{Computing time elapsed
#' in second.} %% \item{comp2 }{Description of 'comp2'} %% ...
#' 
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' ## An Example that runs the NSMM-ICA algorithm on Cohen's tone data
#' data(tonedata, package="mixtools")
#' b <- EMFASTICAALG(tonedata, 2, h=0, tol=1e-8)  
#' 
#' @export EMFASTICAALG
EMFASTICAALG <-
  function (DataMatrix, numCluster, h = 0, 
            maxiter = 300, icaiter = 150, tol = 1e-6, 
            verb = TRUE, combine = TRUE, seednum = 82196, ...){
    RNGstore <- save_rng() # get current random seed info
    
    t0 <- proc.time() # get current time
    
    n <- nrow(DataMatrix) # number of observation
    r <- ncol(DataMatrix) # data dimension
    
    #if(h == 0){
    #  h = 0.9 * (n*r)^{-0.2} * min( sd(DataMatrix), IQR(DataMatrix)/1.34 )
    #}
    
    MemberProbs <- matrix(0, n, numCluster)
    set.seed(seednum)
    preCluster <- kmeans(DataMatrix, numCluster)
    for(i in 1:numCluster)
    {
      MemberProbs[preCluster$cluster==i, i] <- 1
    }
    res <- .Call('icamix_EMInterwovenFastICA', PACKAGE = 'icamix', t(DataMatrix), MemberProbs, h, maxiter+1, icaiter, tol, verb, combine)
    for (i in 1:numCluster){
      res$WMtrs[[i]] <- matrix(res$WMtrs[[i]],r)
      res$WUnmixZ[[i]] <- matrix(res$WUnmixZ[[i]],r)
      res$OriginalSignals[[i]] <- matrix(res$OriginalSignals[[i]],r)
    }
    #for (i in 1:r){
    #  tempresult$Densities[[i]] <- matrix(tempresult$Densities[[i]], numCluster)
    #}
    t1 <- proc.time()
    computingtime <- (t1 - t0)[3]
    cat("\n computing time:", computingtime, " s\n")
    restore_rng(RNGstore)
    file.remove(RNGstore)
    res$call <- match.call()
    res$time <- computingtime
    class(res) <- "EMFASTICAALG"
    res
  }


##
## methods functions
##


#' print.EMFASTICAALG
#' 
#' \code{\link[base]{print}} method for class \code{EMFASTICAALG}.
#' 
#' 
#' @param x An \code{EMFASTICA} object.
#' 
#' @return Returned (invisibly) is the full value of \code{x} itself. %%
#' ~Describe the value returned %% If it is a LIST, use %% \item{comp1
#' }{Description of 'comp1'} %% \item{comp2 }{Description of 'comp2'} %% ...
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' ## An Example that runs the NSMM-ICA algorithm on Cohen's tone data
#' data(tonedata, package="mixtools")
#' 
#' b <- EMFASTICAALG(tonedata, 2, h=0, tol=1e-8)
#' print(b)
#' 
#' @export print.EMFASTICAALG
print.EMFASTICAALG <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\nNumber of Observations: ")
  cat(ncol(x$InputData))
  cat("\nAttributes per Observation: ")
  cat(nrow(x$InputData))
  cat("\nNumber of Mixture Components: ")
  cat(ncol(x$Lambdas))
  cat("\nNumber of Iterations: ")
  cat(nrow(x$Lambdas)-1)
  cat("\nComputing Time Elapsed (sec): ")
  cat(as.numeric(x$time))
  cat("\nData Likelihood from Last Interation:\n")
  cat(x$ObjValue[nrow(x$ObjValue),1])
  cat("\nEstimated Lambdas:\n")
  cat(x$Lambdas[nrow(x$Lambdas),])
  invisible(x)
}



#' summary.EMFASTICAALG
#' 
#' \code{\link[base]{summary}} method for class \code{EMFASTICAALG}.
#' 
#' 
#' @param object An EMFASTICAALG object.
#' 
#' @return The returned value is a "summary.EMFASTICAALG" object which consists
#' a list: %% ~Describe the value returned %% If it is a LIST, use
#' \item{$inputData }{A matrix of which the columns are data entries. Its
#' dimension is \code{r} by \code{n}.} \item{$originSig }{List of Recovered ICA
#' components for each of the \code{m} clusters.} \item{$call }{The function
#' call which results in the corresponding EMFASTICAALG object.} \item{$time
#' }{Computing time elaped in second.} \item{$numIter }{Total number of
#' iterations.} \item{$lastObj }{Objective function value from the last
#' iteration.} \item{$compMeans }{Means of each mixture component.}
#' \item{$compVars }{Covariances of each mixture component.} \item{$numObs
#' }{Total number of observations.} \item{$numAtr }{Dimension of data points.}
#' \item{$numCls }{Number of mixture components.} \item{$estWts }{Estimated
#' mixing weights.} \item{$estCls }{A factor whose levels represent estimated
#' class membership.} \item{$dataLklhd }{A vector storing data loglikelihood
#' from each iteration.} %% ...
#' @examples
#' 
#' ## An Example that runs the NSMM-ICA algorithm on Cohen's tone data
#' data(tonedata, package="mixtools")
#' 
#' b <- EMFASTICAALG(tonedata, 2, h=0, tol=1e-8)
#' summary(b)
#' 
#' @export summary.EMFASTICAALG
summary.EMFASTICAALG <- function(object, ...)
{
  # extract simple info from EMFASTICAALG result
  n <- ncol(object$InputData)
  r <- nrow(object$InputData)
  m <- ncol(object$Lambdas)
  t <- nrow(object$Lambdas)-1
  obj <- object$ObjValue[nrow(object$ObjValue),1]
  lambdas <- object$Lambdas[nrow(object$Lambdas),]
  estimatedclass <- ESTIMATEDMEMBER(object) # estimated class info
  datalikelihood <- object$ObjValue[-1,]
  
  # mixture component means and variances
  compmeans <- list()
  compvars <- list()
  for(j in 1:m){
    compmeans[[j]] <- object$InputData %*% object$MembershipProbs[,j] / sum(object$MembershipProbs[,j])
    centeredData <- object$InputData - compmeans[[j]] %*% rep(1,n)
    compvars[[j]] <- centeredData %*% diag(object$MembershipProbs[,j]) %*% t(centeredData) / ( (sum(object$MembershipProbs[,j]))^2 - sum((object$MembershipProbs[,j])^2) ) * sum(object$MembershipProbs[,j])
  }
  
  # compile the summary list and return
  res <- list(inputData=object$InputData, originSig=object$OriginalSignals, call=object$call, time=as.numeric(object$time), numIter=t, lastObj=obj, compMeans=compmeans, compVars=compvars, numObs=n, numAtr=r, numCls=m, estWts=lambdas, estCls=estimatedclass, dataLklhd=datalikelihood)
  class(res) <- "summary.EMFASTICAALG"
  res
}



#' print.summary.EMFASTICAALG
#' 
#' \code{\link[base]{print}} method for class \code{summary.EMFASTICAALG}.
#' 
#' 
#' @param x An \code{summary.EMFASTICA} object.
#' 
#' @return Returned (invisibly) is the full value of \code{x} itself. %%
#' ~Describe the value returned %% If it is a LIST, use %% \item{comp1
#' }{Description of 'comp1'} %% \item{comp2 }{Description of 'comp2'} %% ...
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' ## An Example that runs the NSMM-ICA algorithm on Cohen's tone data
#' data(tonedata, package="mixtools")
#' 
#' b <- EMFASTICAALG(tonedata, 2, h=0, tol=1e-8)
#' print(summary(b))
#' 
#' @export print.summary.EMFASTICAALG
print.summary.EMFASTICAALG <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\n")
  cat(x$numObs)
  cat(" Observations, ")
  cat(x$numAtr)
  cat(" Attributes, and ")
  cat(x$numCls)
  cat(" Classes.\n")
  cat("Computation Lasts for ")
  cat(x$time)
  cat(" Seconds and ")
  cat(x$numIter)
  cat(" Iterations.")
  cat("\nData Likelihood from Last Interation:\n")
  cat(x$lastObj)
  cat("\nEstimated Lambdas:\n")
  cat(x$estWts)
  cat("\nHard Classification: ")
  print(table(x$estCls))
  cat("\nMean and Covariance of Mixture Component:\n")
  for(j in 1:(x$numCls)){
    cat(paste("Component #",j,"\n Mean:\n"))
    print(x$compMeans[[j]])
    cat(paste(" Covariance:\n"))
    print(x$compVars[[j]])
  }
  invisible(x)
}



#' plot.summary.EMFASTICAALG
#' 
#' \code{plot} method for class \code{summary.EMFASTICAALG}.
#' 
#' 
#' @param x An \code{summary.EMFASTICAALG} object.
#' @param vec1 An integer vector of length two specifying the coordinates with
#' respect to which the data is scatter plotted.
#' @param vec2 An integer vector of length two specifying the coordinates with
#' respect to which the original signal for each mixture component is scatter
#' plotted.
#' 
#' @return Returned (invisibly) is the full value of \code{x} itself. %%
#' ~Describe the value returned %% If it is a LIST, use %% \item{comp1
#' }{Description of 'comp1'} %% \item{comp2 }{Description of 'comp2'} %% ...
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' ## An Example that runs the NSMM-ICA algorithm on Cohen's tone data
#' data(tonedata, package="mixtools")
#' 
#' b <- EMFASTICAALG(tonedata, 2, h=0, tol=1e-8)
#' plot(summary(b))
#' 
#' 
#' @export plot.summary.EMFASTICAALG
plot.summary.EMFASTICAALG <- function(x, vec1=c(1:2), vec2=c(1:2), ...){
  # plot data likelihood vs interation
  if(length(x$dataLklhd) > 1){
    plot.ts(x$dataLklhd,xlab="iteration",ylab="data loglikelihood",main="Objective Function Values")
  }
  
  # plot data with color according to NSMM-ICA classification
  plot((t(x$inputData))[,vec1], asp = 1, col=as.numeric(levels(x$estCls))[x$estCls]+1, pch = as.numeric(levels(x$estCls))[x$estCls], main = "NSMM-ICA Classification",cex.lab=1.1)
  
  # plot recovered independent signals
  par(mfrow = c(1,x$numCls))
  for(clsind in 1:(x$numCls)){
    keyind <- rank(levels(x$estCls))[clsind]
    plot(t(x$originSig[[keyind]])[x$estCls==clsind,vec2], asp = 1, col = clsind + 1, pch = clsind, xlab = 'X', ylab = 'Y', main = paste("Comp #", clsind," IndSig"))
  }
  par(mfrow=c(1,1)) # reset par
  invisible(x)
}



#' plot.EMFASTICAALG
#' 
#' \code{plot} method for class \code{EMFASTICAALG}.
#' 
#' 
#' @param x An \code{EMFASTICAALG} object.
#' @param vec1 An integer vector of length two specifying the coordinates with
#' respect to which the data is scatter plotted.
#' @param vec2 An integer vector of length two specifying the coordinates with
#' respect to which the original signal for each mixture component is scatter
#' plotted.
#' 
#' @return Returned (invisibly) is the full value of \code{x} itself. %%
#' ~Describe the value returned %% If it is a LIST, use %% \item{comp1
#' }{Description of 'comp1'} %% \item{comp2 }{Description of 'comp2'} %% ...
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' ## An Example that runs the NSMM-ICA algorithm on Cohen's tone data
#' data(tonedata, package="mixtools")
#' 
#' b <- EMFASTICAALG(tonedata, 2, h=0, tol=1e-8)
#' plot(b)
#' 
#' @export plot.EMFASTICAALG
plot.EMFASTICAALG <- function(x, vec1=c(1:2), vec2=c(1:2), ...){
  plot(summary(x), vec1, vec2)
  invisible(x)
}









