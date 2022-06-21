## cv


#' Detect and estimate multiple change points in the sparse regression model,
#' selecting the number of change points using sample splitting
#'
#' @param X design matrix
#' @param y response vector
#' @param G integer bandwidth; defaults to \code{round(30 + ncol(X)/100)}
#' @param lambda vector of numeric regularisation parameters
#' @param max.cps maximum number of change points to consider
#' @param family response type, one of \code{"gaussian","binomial","poisson"}
#' @param path.length number of \code{lambda} values to consider
#' @param grid.resolution controls number of subsamples to take
#' @param nu numeric localisation tuning parameter
#' @param do.plot Boolean - return plots
#' @param do.scale Boolean - centre and scale \code{X,y}
#' @param do.refinement Boolean - perform location refinement
#' @param ncores number of parallel cores
#' @param ... optional arguments to \code{glmnet}
#' @return List containing
#' \itemize{
#'   \item{\code{mosum}}{ numeric vector of mosum detector}
#'   \item{\code{cps}}{ integer vector of estimated change points}
#'   \item{\code{refined_cps}}{ integer vector of refined change points}
#'   \item{\code{lambda}}{ selected regularisation parameter}
#'   \item{\code{threshold}}{ implied threshold}
#'   \item{\code{detectors}}{ mosum detector series for each \code{lambda} value}
#'   \item{\code{cv}}{ matrix of cross-validation errors}
#'   \item{\code{model_list}}{ list of fitted piecewise models}
#'   \item{\code{plots}}{ list of detector and refinement plots}
#'   \item{\code{family}}{ input}
#' }
#' @importFrom graphics abline legend matplot
#' @export
#'
#' @examples
#' eqX <- eqdata[,-c(1,9)]
#' eq_thr <- mosumsr.cv(as.matrix(eqX), as.matrix(eqdata[,9]), G=120, max.cps = 3, ncores = 2)
mosumsr.cv <- function(X, y, G = NULL, lambda = NULL, max.cps = NULL, family = c("gaussian","binomial","poisson"),
                       path.length = 5, grid.resolution = 1/G, nu = 0.5, do.plot = TRUE, do.scale = TRUE, do.refinement = TRUE,
                       ncores = NULL, ...){
  n <- dim(X)[1]
  p <- dim(X)[2]
  X <- as.matrix(X)
  if(do.scale) {
    X <- scale(X)
    y <- scale(y)
  }
  if(is.null(G)) G <- 30 + p/100
  family <- match.arg(family, c("gaussian","binomial","poisson"))
  if(is.null(lambda)){
    lambda.max <- max( abs(t(y - mean(y)*(1-mean(y)) ) %*% X ) )/n #/2
    lambda <- round(exp(seq(log(lambda.max), log(lambda.max * .0001), length.out = path.length)), digits = 10)
  }
  path.length <- length(lambda)
  max.cps <- min(max.cps, floor(n/G))
  out_cv <- matrix(Inf,path.length, max.cps+1)
  ranks <- refined_cps <- cps <- out_list <- list()


  ms <- get.cv.detectors(X,y,G,lambda,family, grid.resolution = grid.resolution, ncores = ncores, ...)
  for (ll in 1:path.length) {
    cps[[ll]] <- get_local_maxima(ms$mosum[,ll], 0, G, nu)
    q <- length(cps[[ll]])
    refined_cps[[ll]] <- cps[[ll]]
    if(do.refinement) for (k in 1:q) {
      L_ind <- max(1,cps[[ll]][k]-G- floor(G/2))
      G_window <- (L_ind):(L_ind+G-1)# (cps[[ll]][k]-1- floor(G/2))
      L_mod <- glmnet(X[G_window,], y[G_window], family = family, ...)
      R_ind <- min(n-G,cps[[ll]][k]+ floor(G/2))
      G_window <- (R_ind):(R_ind+G)
      R_mod <- glmnet(X[G_window,], y[G_window], family = family, ...)
      rf <- refinement(X, y, cps[[ll]][k], G, lambda[ll], L_mod,  R_mod, family=family)
      refined_cps[[ll]][k] <- rf$cp
    }
    ranks[[ll]] <- rank(-ms$mosum[cps[[ll]],ll])

    out_list[[ll]] <- fit.cv(X,y,refined_cps[[ll]],ranks[[ll]],max.cps, lambda[ll])

    out_cv[ll,1:length(out_list[[ll]]$error)] <- out_list[[ll]]$error
  }
  min.point <- which(out_cv==min(out_cv), arr.ind=TRUE)[1,]
  out_q <- (0:max.cps)[min.point[2]]
  out_lambda <- lambda[min.point[1]]
  out_cps <- cps[[min.point[1]]]
  out_cps <- out_cps[ranks[[min.point[1]]]<=out_q]
  out_refined_cps <- refined_cps[[min.point[1]]]
  out_refined_cps <- out_refined_cps[ranks[[min.point[1]]]<=out_q]

  if(out_q==0) out_refined_cps <- out_cps <- NULL

  pl <- list()
  thr <- Inf
  if(do.plot) {
    par(mfrow=c(2,1))
    matplot(0:max.cps, t(out_cv), type = 'b', col = 1+1:path.length, pch = 1+1:path.length,
            xlab = 'q', ylab = 'CV', main = 'CV for change point number estimation', log = "y")
    abline(v = out_q)
    legend('topleft', legend = lambda, col = 1+1:path.length, pch = 1+1:path.length, lty = 1, cex = .6)

    if(grid.resolution == 1/G) plot.ts(ms$mosum[,min.point[1]], ylab="Detector", xlab = "Time") # plot test statistic
    else plot(which(ms$mosum[,min.point[1]]>0), ms$mosum[which(ms$mosum>0),min.point[1]], ylab="Detector",
              xlab = "Time", ylim = c(0, max(ms$mosum[,min.point[1]])), xlim = c(1,n) )
    if(out_q>0) { # add estimated cps
      thr <- ms$mosum[out_cps[which(ranks[[min.point[1]]] == out_q)],min.point[1]]
      abline(h = thr, col = "blue") #add threshold
      abline(v = out_cps, col = "red")
      abline(v = out_refined_cps, col = "purple")
    }
    par(mfrow =c(1,1))
    pl$mosum <- recordPlot()
  }

  out <- list(cps = out_cps, refined_cps = out_refined_cps,
              threshold = thr, lambda = out_lambda,
              detectors = ms, cv= out_cv, model_list = out_list,
              plots = pl, family = family)
  attr(out, "class") <- "mosumsr"
  return(out)
}
# try_cv <- mosumsr.cv(as.matrix(eqX), as.matrix(eqdata[,9]), G=240, ncores = 4, path.length = 10, max.cps = 4)

#' @title get mosum detectors
#' @keywords internal
get.cv.detectors <- function(X, y, G, lambda, family = c("gaussian","binomial","poisson"), grid.resolution = 1/G, ncores = NULL, ...){
  n <- dim(X)[1]
  p <- dim(X)[2]
  X <- as.matrix(X)
  family <- match.arg(family, c("gaussian","binomial","poisson"))

  path.length <- length(lambda)
  ##parallel
  if(is.null(ncores)) ncores <- max(1, detectCores() -1)
  cl <- makeCluster(ncores)
  registerDoParallel(cl)

  model_list <- list() #models
  coeffs <- array(0, c(n-G+1, p+1, path.length) ) #coefficients
  mosum <- matrix(0, n, path.length) #detector

  T_Grid <- indexLE(n, G, kap=grid.resolution)
  model_list_unsorted <- foreach (t = T_Grid, .packages = "glmnet") %dopar%{
    G_window <- t:(t+G-1)
    glmnet(X[G_window,], y[G_window], family = family, ...)
  }
  for (ii in 1:length(T_Grid) ) {
    t <- T_Grid[ii]
    model_list[[t]] <- model_list_unsorted[[ii]]
  }
  for (ii in 1:length(T_Grid) ) {
    t <- T_Grid[ii]
    coeffs[t,,] <- as.vector(coef.glmnet( model_list[[t]], lambda))
  }
  for (ll in 1:length(lambda)) {
    for (t in T_Grid[T_Grid<=n-2*G+1])  {
      mosum[G+t,ll] <- sqrt(G/2) * norm(coeffs[t+G,,ll] - coeffs[t,,ll], "2")
    }
  }
  stopCluster(cl)
  return(list(coeffs=coeffs, mosum=mosum))
}
# try_cv_detectors <- get.cv.detectors(as.matrix(eqX), as.matrix(eqdata[,9]), G=120, c(.01,.1), ncores = 2)




#' @title fit models with cross-validation
#' @keywords internal
fit.cv <- function(X, y, cps, ranks, max.cps, lambda, family =  c("gaussian","binomial","poisson"), do.plot = TRUE, do.scale = TRUE, ...){
  min.len <- 3
  cps <- sort((cps)) #unique
  q <- min(length(cps),max.cps)
  cps <- cps[ranks <= q]
  ranks <- ranks[ranks <= q]
  family <- match.arg(family, c("gaussian","binomial","poisson"))
  #type <- match.arg(type, c("link","response", "coefficients", "nonzero", "class"))
  n <- nrow(X)
  p <- ncol(X)
  X <- as.matrix(X)
  if(do.scale) {
    X <- scale(X)
    y <- scale(y)
  }

  # if(is.character(lambda))  lambda <-  match.arg(lambda, c("min","1se"))
  starts <- c(0, cps); ends <- c(cps, n)
  cps2 <- floor(cps/2)
  starts2 <- c(0, cps2); ends2 <- c(cps2, floor(n/2) )

  out <- list()
  preds <- matrix(0, n, q+1)
  coeffs <- matrix(0, p+1,n)
  if(!is.null(colnames(X))){
    rownames(coeffs)[1:p+1] <- colnames(X)
    rownames(coeffs)[1] <- "intercept"
  }
  error <- rep(0, q+1)

  evens <- which(as.logical(1:n %%2) )
  odds <- which(!(1:n %% 2) )
  X.odd <- X[odds,]; y.odd <- y[odds]; y.even <- y[evens]

  ii <- q+1
  ii_cps <- cps2[ranks<ii]
  out[[ii]] <- list()
  if(min(ends2 - starts2) > min.len  ){
    for (jj in 1:length(ii_cps)) {

      out[[ii]][[jj]] <- glmnet(X.odd[(starts2[jj]+1):ends2[jj],], y.odd[(starts2[jj]+1):ends2[jj]],nfolds = 10,family=family,...)
      preds[(starts[jj]+1):ends[jj],ii] <- predict(out[[ii]][[jj]], X[(starts[jj]+1):ends[jj],], lambda )
    }
  } else warning("Segment too short. Returning Inf")
  for (ii in q:1 ) {
    ii_cps <- cps2[ranks<ii]
    out[[ii]] <- out[[ii+1]]
    preds[,ii] <- preds[,ii+1]
    if( min(ends2 - starts2) > min.len ){
      added_cp <- ranks[ii]#which.max(ranks) #which( !(cps2[ranks<ii+1] %in%  ii_cps)) #cps2[ranks==ii+1]
      out[[ii]][[added_cp]] <- out[[ii]][[added_cp+1]] <- #remove one cp
        glmnet(X.odd[(starts2[added_cp]+1):ends2[added_cp+1],], y.odd[(starts2[added_cp]+1):ends2[added_cp+1]],nfolds = 10,family=family,...)
      preds[(starts[added_cp]+1):ends[added_cp+1],ii] <- predict(out[[ii]][[added_cp]], X[(starts[added_cp]+1):ends[added_cp+1],], lambda )
    } else warning("Segment too short. Returning Inf")

  }
  for (ii in (q+1):1 ) {
    error[ii] <- sum(likelihood(y.even, preds[evens,ii], family))
  }

  out_list <- list(model = out, error = error, preds = preds, coeffs = t(coeffs))
  return(out_list)
}
# try_fit_cv <- fit.cv(as.matrix(eqX), cps = c(100,200), ranks = c(1,2), max.cps = 2, as.matrix(eqdata[,9]),lambda = .1)













#' Detect and estimate multiple change points in the sparse regression model using multiple bandwidths,
#' selecting the number of change points using sample splitting
#'
#'
#' @param X design matrix
#' @param y response vector
#' @param Gset integer vector of bandwidths; default smallest is \code{round(30 + ncol(X)/100)}
#' @param lambda  regularisation parameter; either a numeric, or one of \code{"min","1se"} (see \link[glmnet]{cv.glmnet})
#' @param family response type, one of \code{"gaussian","binomial","poisson"}
#' @param threshold numeric test rejection threshold; see reference for default choice
#' @param grid.resolution controls number of subsamples to take
#' @param nu numeric localisation tuning parameter
#' @param do.plot Boolean - return plots
#' @param do.scale Boolean - scale \code{X,y}
#' @param ncores number of parallel cores
#' @param ... optional arguments to \code{glmnet}
#'
#' @return List containing
#' \itemize{
#'   \item{\code{anchors}}{ integer vector of estimated change points}
#'   \item{\code{refined_cps}}{ integer vector of refined change points}
#'   \item{\code{mosumsr.G}}{ list of `mosumsr` objects corresponding to `Gset` in ascending order}
#' }
#' @export
#'
#' @examples
#' eqX <- eqdata[,-c(1,9)]
#' eq_mosum <- mosumsr.multiscale.cv(as.matrix(eqX), eqdata[,9], c(60,90,120), ncores = 2)
mosumsr.multiscale.cv <- function(X, y, Gset = NULL, lambda = NULL, family = c("gaussian","binomial","poisson"),
                                  threshold = NULL, grid.resolution = 1/Gset, nu = 0.5, do.plot = TRUE, do.scale = TRUE,
                                  ncores = NULL, ...){
  n <- dim(X)[1]
  p <- dim(X)[2]
  X <- as.matrix(X)
  if(do.scale) {
    X <- scale(X)
    y <- scale(y)
  }
  if(is.null(Gset)) {
    G <- 30 + p/100
    Gset <- c(G, 2*G, 3*G)
  }
  ## validate inputs
  if ( !( all(Gset > 0) )) stop("All entries of Gset must be positive")
  if ( !( all(Gset <= n/2) )) stop("All entries of Gset must be at most n/2")
  # if(is.character(lambda))  lambda <-  match.arg(lambda, c("min","1se"))
  family <-  match.arg(family, c("gaussian","binomial","poisson"))
  if ( !(family %in% c("gaussian","binomial","poisson"))) stop("family must be \"gaussian\", \"binomial\", or \"poisson\" ")
  # if ( !(threshold.constant >= 0)) stop("threshold.constant must be at least 0")
  if ( !(nu > 0)) stop("nu must be positive")
  # if ( !(threshold.log.constant >= 0)) stop("threshold.log.constant must be at least 0")
  Gset <- sort(Gset) #into ascending order
  Glen <- length(Gset)
  anchors <-  c()
  mosumsr.G <- as.list(1:Glen)
  Reject <- 0
  for (ii in 1:Glen){
    mosumsr.G[[ii]] <-  mosumsr.cv(X, y, G= Gset[ii], lambda =  lambda, family = family,
                                   grid.resolution = grid.resolution[ii], nu = nu,
                                   do.refinement = FALSE, do.plot = do.plot, ncores = ncores, ...)
    if(length(mosumsr.G[[ii]]$cps)>0) Reject <- 1
  }
  refined_cps <- NULL
  if(Reject){
    interval <- matrix(0, n, n) #matrix(0, Glen, n) #detection set
    # anchor sets
    anchors <- mosumsr.G[[1]]$cps
    if(!is.null(anchors)){
      for (ii in 1:length(anchors)) {
        L_ <- max(1, anchors[ii]-Gset[1]+1)
        R_ <- min(n, anchors[ii]+Gset[1])
        interval[ii,L_:R_] <- 1
      }
      Gout <- rep(Gset[1], length(anchors))
    }
    if(Glen > 1){
      for (ii in 2:Glen){
        K <- mosumsr.G[[ii]]$cps
        if(length(K)>0){
          if(is.null(anchors)){
            anchors <- K
            for (jj in 1:length(K)) {
              L_ <- max(1, (K[jj]- floor(nu* Gset[ii])+1))
              R_ <- min(n, (K[jj]+floor(nu* Gset[ii])))
              interval_jj <- L_:R_
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
        K <- mosumsr.G[[ii]]$cps
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
    if(is.null(lambda)) lambda <- mosumsr.G[[1]]$lambda #| lambda %in% c("1se","min")
    else lambda <- lambda[1]
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
  # if(do.plot){
  #   par(mfrow = c(2*Glen,1))
  #   for (ii in 1:Glen) {
  #     mosumsr.G[[ii]]$plots
  #     # plot.ts(mosumsr.G[[ii]]$detectors$mosum[,], ylab="Detector") # plot series
  #     # abline(h = mosumsr.G[[ii]]$threshold, col = "blue") #add threshold
  #     # if(Reject) {
  #     #   abline(v = mosumsr.G[[ii]]$cps, col = "red")  #add estimated cps
  #     #   if(ii == Glen) abline(v = refined_cps, col = "purple")
  #     # }
  #   }
  #   par(mfrow = c(1,1))
  #   pl <- recordPlot()
  # } else pl <- NULL

  out <- list(anchors= anchors, refined_cps = refined_cps,  mosumsr.G =mosumsr.G)
  attr(out, "class") <- "mosumsr.ms"
  return(out)
}
# try_ms_cv <- mosumsr.multiscale.cv(as.matrix(eqX), as.matrix(eqdata[,9]), Gset=c(120,240), ncores = 4, path.length = 10, max.cps = 4)
