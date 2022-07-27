## multiscale

#' Detect and estimate multiple change points in the sparse regression model using multiple bandwidths
#'
#' @param X design matrix
#' @param y response vector
#' @param Gset integer vector of bandwidths; default smallest is \code{round(30 + ncol(X)/100)}
#' @param lambda  regularisation parameter; either a numeric, or one of \code{"min","1se"} (see \link[glmnet]{cv.glmnet})
#' @param family response type, one of \code{"gaussian","binomial","poisson"}
#' @param threshold numeric test rejection threshold; see  \link[moseg]{moseg} for default choice
#' @param grid.resolution controls number of subsamples to take
#' @param nu numeric localisation tuning parameter
#' @param do.plot Boolean - return plots
#' @param do.scale Boolean - scale \code{X,y}
#' @param ncores number of parallel cores
#' @param ... optional arguments to \code{glmnet}
#'
#' @return List containing
#' \itemize{
#'   \item{\code{cps}}{ integer vector of estimated change points}
#'   \item{\code{plot}}{ multiscale plot}
#'   \item{\code{moseg.G}}{ list of `moseg` objects corresponding to `Gset` in ascending order}
#' }
#' @export
#'
#' @seealso  \link[moseg]{moseg.ms.cv}
#' @examples
#' eqX <- eqdata[,-c(1,9)]
#' eq_mosum <- moseg.ms(as.matrix(eqX), eqdata[,9], c(60,90,120), ncores = 2)
moseg.ms <- function(X, y, Gset, lambda = c("min","1se"), family = c("gaussian","binomial","poisson"),
                               threshold = NULL, grid.resolution = NULL, nu = 0.5, do.plot = TRUE, do.scale = TRUE,
                               ncores = NULL, ...){
  n <- dim(X)[1]
  p <- dim(X)[2]
  X <- as.matrix(X)
  if(do.scale) {
    X <- scale(X)
    y <- scale(y)
  }
  if(is.null(Gset)) {
    G <- getG(p, n)
    Gset <- c(G, 4/3*G, 5/3*G)
  }
  ## validate inputs
  if ( !( all(Gset > 0) )) stop("All entries of Gset must be positive")
  if ( !( all(Gset <= n/2) )) stop("All entries of Gset must be at most n/2")
  if(is.character(lambda))  lambda <-  match.arg(lambda, c("min","1se"))
  family <-  match.arg(family, c("gaussian","binomial","poisson"))
  if ( !(family %in% c("gaussian","binomial","poisson"))) stop("family must be \"gaussian\", \"binomial\", or \"poisson\" ")
  # if ( !(threshold.constant >= 0)) stop("threshold.constant must be at least 0")
  if ( !(nu > 0)) stop("nu must be positive")
  # if ( !(threshold.log.constant >= 0)) stop("threshold.log.constant must be at least 0")
  Gset <- sort(Gset) #into ascending order
  Glen <- length(Gset)
  anchors <-  c()
  moseg.G <- as.list(1:Glen)
  Reject <- 0
  for (ii in 1:Glen){
    moseg.G[[ii]] <-  moseg(X, y, G= Gset[ii], lambda =  lambda, family = family, threshold =  threshold,
                                grid.resolution = grid.resolution, nu = nu,
                                do.refinement = F, do.plot =F, ncores = ncores, ...)
    if(length(moseg.G[[ii]]$cps)>0) Reject <- 1
  }
  refined_cps <- NULL
  if(Reject){
    interval <- matrix(0, n, n) #matrix(0, Glen, n) #detection set
    # anchor sets
    anchors <- moseg.G[[1]]$cps
    if(!is.null(anchors)){
      for (ii in 1:length(anchors)) interval[ii,(anchors[ii]-Gset[1]+1):(anchors[ii]+Gset[1])] <- 1
      Gout <- rep(Gset[1], length(anchors))
    }
    if(Glen > 1){
      for (ii in 2:Glen){
        K <- moseg.G[[ii]]$cps
        if(length(K)>0){
          if(is.null(anchors)){
            anchors <- K

            for (jj in 1:length(K)) {
              interval_jj <-(K[jj]- floor(nu* Gset[ii])+1):(K[jj]+floor(nu* Gset[ii]))
              interval[jj,interval_jj] <- 1
              Gout <- rep(Gset[ii], length(anchors))
            }
          } else{
            for (cp in K) {
              cp_window_min <- max(cp-floor(nu* Gset[ii])+1 , 1)
              cp_window_max <- min((cp+floor(nu* Gset[ii])), n)
              cp_window <- cp_window_min:cp_window_max
              if(sum(interval[,cp_window]) == 0) {
                anchors <- append(anchors, cp)
                Gout <- append(Gout, Gset[ii])
                interval[length(anchors),cp_window] <- 1
              }
            }
          }
        }
      }
    }
    q <- length(anchors)
    rank.anchors <- rank(anchors)
    anchors <- anchors[rank.anchors] #sort anchors
    Gout_min <- Gout_max <- Gout[rank.anchors] #sort Gout by anchors
    cps <- anchors
    # cluster sets
    if(Glen > 1){
      for (ii in 2:Glen){
        K <- moseg.G[[ii]]$cps
        if(is.null(cps)) cps <- K
        else for (cp in K) {
          cond1 <- FALSE
          cp_window <- max(cp-Gset[ii]+1,1):min(cp+Gset[ii],n)
          rs <- rowSums(interval[,cp_window] ) > 0
          if(sum(rs)==1) { #intersection nonempty for exactly one anchor
            window_len <- (Gset[ii]+floor(Gset[ii]/2))
            cluster_window <- (cp - window_len +1):(cp+window_len)
            cluster_window <- cluster_window[cluster_window>0 & cluster_window <= n ]
            cond1 <- sum(interval[-which(rs),cluster_window]) == 0 #intersection empty for all other anchors
            if(cond1) Gout_max[rs] <- Gset[ii]
          }
        }
      }
    }
    if(do.plot)  par(mfrow = c(max(q%/%6,1) , max(q%%6,1) ))
    refined_cps <- anchors
    Gstar <- floor(Gout_min*3/4 + Gout_max/4)
    if(is.null(lambda) | lambda %in% c("1se","min") ) lambda <- moseg.G[[1]]$lambda
    for (k in 1:q) {

      G_window <- #(anchors[k] - Gstar[k]):(anchors[k]-1)
        (anchors[k]- floor(Gout_min[k]/2) - Gstar[k]):(anchors[k]- floor(Gout_min[k]/2))
      if(any(G_window<1)) G_window <- 1: Gstar[k]
      lmod <- glmnet(X[G_window,], y[G_window], family = family, ...)

      G_window <-  #(anchors[k]):(anchors[k]+Gstar[k]-1)
        (anchors[k]+ floor(Gout_min[k]/2)):(anchors[k]+ floor(Gout_min[k]/2)+Gstar[k]-1)
      if(any(G_window>n)) G_window <- (n-Gstar[k]+1):n
      rmod <- glmnet(X[G_window,], y[G_window], family = family, ...)

      rf <- refinement(X, y, anchors[k], Gstar[k], lambda, lmod, rmod, family=family)
      refined_cps[k] <- rf$cp

      if(do.plot){
        plot.ts(rf$objective, ylab = "Q", xaxt = "n")
        axis(1, at=1:(2*Gstar[k]), labels= (anchors[k]-Gstar[k]+1):(anchors[k]+Gstar[k]))
        abline(v = which.min(rf$objective), col = "purple")
      }

    }
    if(do.plot) par(mfrow = c(1,1))
  }
  if(do.plot){
    par(mfrow = c(Glen,1))
    for (ii in 1:Glen) {
      plot.ts(moseg.G[[ii]]$mosum, ylab="Detector") # plot series
      abline(h = moseg.G[[ii]]$threshold, col = "blue") #add threshold
      if(Reject) {
        abline(v = moseg.G[[ii]]$cps, col = "red")  #add estimated cps
        if(ii == Glen) abline(v = refined_cps, col = "purple")
      }
    }
    par(mfrow = c(1,1))
    pl <- recordPlot()
  } else pl <- NULL
  out <- list(anchors= anchors, refined_cps = refined_cps, plot=pl, moseg.G =moseg.G)
  attr(out, "class") <- "moseg.ms"
  return(out)
}

