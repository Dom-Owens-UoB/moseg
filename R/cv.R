## cv


#' Detect and estimate multiple change points in the sparse regression model,
#' selecting the number of change points using sample splitting
#'
#' @param X design matrix
#' @param y response vector
#' @param G integer bandwidth; default is chosen based on dimensions of \code{X}
#' @param lambda vector of numeric regularisation parameters
#' @param max.cps maximum number of change points to consider
#' @param family response type, one of \code{"gaussian","binomial","poisson"}
#' @param loss l-norm for CV loss function, one of \code{"1", "2"}
#' @param folds number of folds for CV
#' @param path.length number of \code{lambda} values to consider
#' @param grid.resolution controls number of subsamples to take
#' @param nu numeric localisation tuning parameter
#' @param do.plot Boolean - return plots
#' @param do.scale Boolean - centre and scale \code{X,y}
#' @param do.refinement Boolean - perform location refinement
#' @param ncores number of parallel cores
#' @param ... optional arguments to \code{glmnet}
#' @return \code{moseg.cv} object containing
#' \itemize{
#'   \item{\code{mosum}}{ numeric vector of mosum detector}
#'   \item{\code{cps}}{ integer vector of estimated change points}
#'   \item{\code{refined.cps}}{ integer vector of refined change points}
#'   \item{\code{G}}{ input}
#'   \item{\code{lambda}}{ selected regularisation parameter}
#'   \item{\code{threshold}}{ implied threshold}
#'   \item{\code{detectors}}{ mosum detector series for each \code{lambda} value}
#'   \item{\code{cv}}{ matrix of cross-validation errors}
#'   \item{\code{model_list}}{ list of fitted piecewise models}
#'   \item{\code{plots}}{ list of detector and refinement plots}
#'   \item{\code{family}}{ input}
#' }
#' @import glmnet
#' @import doParallel
#' @importFrom graphics abline legend matplot
#' @export
#'
#' @examples
#' eqX <- eqdata[,-c(1,9)]
#' eq_thr <- moseg.cv(as.matrix(eqX), as.matrix(eqdata[,9]), G=120, max.cps = 3, ncores = 2)
moseg.cv <- function(X, y, G = NULL, lambda = NULL, max.cps = NULL, family = c("gaussian","binomial","poisson"), loss = c("1","2"), folds = 2,
                       path.length = 5, grid.resolution = 1/G, nu = 0.5, do.plot = TRUE, do.scale = TRUE, do.refinement = TRUE,
                       ncores = NULL, ...){
  n <- dim(X)[1]
  p <- dim(X)[2]
  X <- as.matrix(X)
  if(do.scale) {
    X <- scale(X)
    y <- scale(y)
  }
  if(is.null(G)) G <- getG(p, n)
  family <- match.arg(family, c("gaussian","binomial","poisson"))
  loss <- match.arg(loss, c("1","2"))
  if(is.null(lambda)){
    lambda.max <- max( abs(t(y)  %*% X ) )/G
    lambda <- round(exp(seq(log(lambda.max), log(lambda.max * .001), length.out = path.length)), digits = 10)
  }
  path.length <- length(lambda)
  max.cps <- min(max.cps, floor(n/(nu*G)))
  out_cv <- matrix(Inf,path.length, max.cps+1)
  ranks <- refined.cps <- cps <- out_list <- list()


  ms <- get.cv.detectors(X,y,G,lambda,family, grid.resolution = grid.resolution, ncores = ncores, ...)
  for (ll in 1:path.length) {
    cps[[ll]] <- get_local_maxima(ms$mosum[,ll], 0, G, nu)
    q <- length(cps[[ll]])
    refined.cps[[ll]] <- cps[[ll]]
    limits <- c(0, cps[[ll]][-q] + floor(diff(cps[[ll]])/2), n)
    if(do.refinement) for (k in 1:q) {
      L_ind <- max(1,cps[[ll]][k]-G- floor(G/2))
      G_window <- (L_ind):(L_ind+G-1)# (cps[[ll]][k]-1- floor(G/2))
      L_mod <- glmnet(X[G_window,], y[G_window], family = family, ...)
      R_ind <- min(n-G,cps[[ll]][k]+ floor(G/2))
      G_window <- (R_ind):(R_ind+G)
      R_mod <- glmnet(X[G_window,], y[G_window], family = family, ...)
      rf <- refinement(X, y, cps[[ll]][k], G, lambda[ll], L_mod,  R_mod,
                       L_min = limits[k], U_max = limits[k+1],
                       family=family)
      refined.cps[[ll]][k] <- rf$cp
    }
    ranks[[ll]] <- rank(-ms$mosum[cps[[ll]],ll])

    out_list[[ll]] <- fit.cv(X,y,refined.cps[[ll]],ranks[[ll]],max.cps, lambda[ll], family = family, loss = loss, folds = folds)

    out_cv[ll,1:length(out_list[[ll]]$error)] <- out_list[[ll]]$error
  }
  min.point <- which(out_cv==min(out_cv), arr.ind=TRUE)[1,]
  out_q <- (0:max.cps)[min.point[2]]
  out_lambda <- lambda[min.point[1]]
  out_cps <- cps[[min.point[1]]]
  out_cps <- out_cps[ranks[[min.point[1]]]<=out_q]
  out_refined.cps <- refined.cps[[min.point[1]]]
  out_refined.cps <- out_refined.cps[ranks[[min.point[1]]]<=out_q]

  if(out_q==0) {
    out_refined.cps <- out_cps <- NULL
    thr <- Inf
  } else {
    thr <- ms$mosum[out_cps[which(ranks[[min.point[1]]] == out_q)],min.point[1]]
  }
  out <- list(cps = out_cps, refined.cps = out_refined.cps,
              G=G, threshold = thr, lambda = out_lambda,
              detectors = ms, cv = out_cv, model_list = out_list,
              family = family)
  attr(out, "pl") <- list(max.cps=max.cps, lambda=lambda, min.point=min.point,ranks=ranks)
  attr(out, "class") <- "moseg.cv"
  pl <- list()
  if(do.plot){
    par(mfrow=c(2,1))
    plot(out, type = "cv")
    plot(out, type = "mosum")
    par(mfrow =c(1,1))
    pl$mosum <- recordPlot()

  }
  out$plots <- pl
  return(out)
}


#' @title Plot the moseg detector
#' @method plot moseg.cv
#' @description Plotting method for S3 objects of class \code{moseg.cv}.
#' @param x \code{moseg.cv} object
#' @param type plot type, one of \code{"cv","mosum"}
#' @param ... unused
#'
#' @return A cv or detector plot
#' @export
plot.moseg.cv <- function(x, type = c("cv","mosum"), ...){
  type <- match.arg(type, c("cv","mosum"))
  pln <- attributes(x)$pl
  max.cps <- pln$max.cps
  lambda <- pln$lambda
  path.length <- length(lambda)
  min.point <- pln$min.point
  ranks <- pln$ranks

  if(type == "cv") {
    matplot(0:max.cps, t(x$cv), type = 'b', col = 1+1:path.length, pch = 1+1:path.length,
            xlab = 'q', ylab = 'CV', main = 'CV for change point number estimation', log = "y")
    abline(v = length(x$cps))
    options(scipen = 2)
    legend('bottomright', legend = lambda, col = 1+1:path.length, pch = 1+1:path.length, lty = 1, cex = .6, title = "lambda")
  }
  if(type == "mosum"){
    if(x$detectors$grid.resolution == 1/x$G) plot.ts(x$detectors$mosum[,min.point[1]], ylab="Detector", xlab = "Time") # plot test statistic
    else plot(which(x$detectors$mosum[,min.point[1]]>0), x$detectors$mosum[which(x$detectors$mosum>0),min.point[1]], ylab="Detector",
              xlab = "Time", ylim = c(0, max(x$detectors$mosum[,min.point[1]])), xlim = c(1,length(x$detectors$mosum[,1])))
    if(length(x$cps)>0) { # add estimated cps
      thr <- x$detector$mosum[x$cps[which(ranks[[min.point[1]]] == x$cps)],min.point[1]]
      abline(h = thr, col = "blue") #add threshold
      abline(v = x$cps, col = "red")
      abline(v = x$refined.cps, col = "purple")
    }
  }
}

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
  return(list(coeffs=coeffs, mosum=mosum, grid.resolution=grid.resolution))
}




#' @title fit models with cross-validation
#' @keywords internal
fit.cv <- function(X, y, cps, ranks, max.cps, lambda, family =  c("gaussian","binomial","poisson"), loss = c("1","2"), folds = 2,
                   do.plot = TRUE, do.scale = TRUE, ...){
  loss <- match.arg(loss, c("1","2"))
  min.len <- 3
  rankmat <- cbind(cps, ranks)
  rankmat2 <- rankmat[1,, drop = FALSE]
  if(nrow(rankmat) > 1) for (ii in 1:(nrow(rankmat)-1)) {
    if(rankmat[ii+1,1] - rankmat[ii,1] > min.len) rankmat2 <- rbind(rankmat2, rankmat[ii+1,])
  }
  cps <- rankmat2[,1] #sort((cps)) #unique
  ranks <- rank(rankmat2[,2])
  q <- min(length(cps),max.cps)
  # cps <- cps[ranks <= q]
  # cps <- sort((cps))
  # ranks <- ranks[ranks <= q]
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
  cps2 <- floor(cps/max(2,folds) )
  starts2 <- c(0, cps2); ends2 <- c(cps2, floor(n/max(2,folds)) )

  out <- list()
  preds <- matrix(0, n, q+1)
  coeffs <- matrix(0, p+1,n)
  if(!is.null(colnames(X))){
    rownames(coeffs)[1:p+1] <- colnames(X)
    rownames(coeffs)[1] <- "intercept"
  }
  error <- rep(0, q+1)




  for (fold in 1:folds) {
    test <- which(as.logical( (1:n + fold - 1) %% max(2,folds) == 0) )
    train <- which(!1:n  %in% test)#which(!(1:n %% 2) )
    X.test <- X[test,]; y.test <- y[test]; X.train <- X[train,]; y.train <- y[train]

    ii <- q+1
    ii_cps <- cps2[ranks<ii]
    out[[ii]] <- list()
    # if(fold == 1) {
    #   X.train <- X.odd; y.train <- y.odd; y.test <- y.even;  test <- evens
    # }
    # if (fold == 2){
    #   X.train <- X.even; y.train <- y.even; y.test <- y.odd;  test <- odds
    # }

    if(min(ends2 - starts2) > min.len  ){
      for (jj in 1:length(ii_cps)) {
        out[[ii]][[jj]] <- glmnet(X.train[(starts2[jj]+1):ends2[jj],], y.train[(starts2[jj]+1):ends2[jj]],nfolds = 10,family=family,...)
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
          glmnet(X.train[(starts2[added_cp]+1):ends2[added_cp+1],], y.train[(starts2[added_cp]+1):ends2[added_cp+1]],nfolds = 10,family=family,...)
        preds[(starts[added_cp]+1):ends[added_cp+1],ii] <- predict(out[[ii]][[added_cp]], X[(starts[added_cp]+1):ends[added_cp+1],], lambda )
      } else warning("Segment too short. Returning Inf")

    }
    if(loss == "1") for (ii in (q+1):1 ) {
      error[ii] <- error[ii] + sum(abs(y.test - preds[test,ii])) ##
    }

    if(loss == "2") for (ii in (q+1):1 ) {
      error[ii] <- error[ii] + sum(likelihood(y.test, preds[test,ii], family)) ##
    }
  }
  out_list <- list(error=error) #list(model = out, error = error, preds = preds, coeffs = t(coeffs))
  return(out_list)
}












#' Detect and estimate multiple change points in the sparse regression model using multiple bandwidths,
#' selecting the number of change points using sample splitting
#'
#'
#' @param X design matrix
#' @param y response vector
#' @param Gset integer vector of bandwidths; default is chosen based on dimensions of \code{X}
#' @param lambda  regularisation parameter; either a numeric, or one of \code{"min","1se"} (see \link[glmnet]{cv.glmnet})
#' @param family response type, one of \code{"gaussian","binomial","poisson"}
#' @param loss l-norm for CV loss function, one of \code{"1", "2"}
#' @param folds number of folds for CV
#' @param threshold numeric test rejection threshold; see reference for default choice
#' @param grid.resolution controls number of subsamples to take
#' @param nu numeric localisation tuning parameter
#' @param do.plot Boolean - return plots
#' @param do.scale Boolean - scale \code{X,y}
#' @param ncores number of parallel cores
#' @param ... optional arguments to \code{glmnet}
#'
#' @return \code{moseg.ms.cv} object containing
#' \itemize{
#'   \item{\code{anchors}}{ integer vector of estimated change points}
#'   \item{\code{refined.cps}}{ integer vector of refined change points}
#'   \item{\code{moseg.G}}{ list of `moseg.cv` objects corresponding to `Gset` in ascending order}
#' }
#' @import glmnet
#' @import doParallel
#' @export
#'
#' @examples
#' eqX <- eqdata[,-c(1,9)]
#' eq_mosum <- moseg.ms.cv(as.matrix(eqX), eqdata[,9], c(60,90,120), ncores = 2)
moseg.ms.cv <- function(X, y, Gset = NULL, lambda = NULL, family = c("gaussian","binomial","poisson"), loss = c("1","2"), folds = 1,
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
    G <- getG(p, n)
    Gset <- c(G, 4/3*G, 5/3*G)
  }
  ## validate inputs
  if ( !( all(Gset > 0) )) stop("All entries of Gset must be positive")
  if ( !( all(Gset <= n/2) )) stop("All entries of Gset must be at most n/2")
  # if(is.character(lambda))  lambda <-  match.arg(lambda, c("min","1se"))
  family <-  match.arg(family, c("gaussian","binomial","poisson"))
  # if ( !(threshold.constant >= 0)) stop("threshold.constant must be at least 0")
  if ( !(nu > 0)) stop("nu must be positive")
  # if ( !(threshold.log.constant >= 0)) stop("threshold.log.constant must be at least 0")
  loss <- match.arg(loss, c("1","2"))
  Gset <- sort(Gset) #into ascending order
  Glen <- length(Gset)
  anchors <-  c()
  moseg.G <- as.list(1:Glen)
  Reject <- 0
  for (ii in 1:Glen){
    moseg.G[[ii]] <-  moseg.cv(X, y, G= Gset[ii], lambda =  lambda, family = family, loss = loss, folds = folds,
                                   grid.resolution = grid.resolution[ii], nu = nu,
                                   do.refinement = FALSE, do.plot = FALSE, ncores = ncores, ...)
    if(length(moseg.G[[ii]]$cps)>0) Reject <- 1
  }
  refined.cps <- NULL
  if(Reject){
    interval <- matrix(0, n, n) #matrix(0, Glen, n) #detection set
    # anchor sets
    anchors <- moseg.G[[1]]$cps
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
        K <- moseg.G[[ii]]$cps
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
    anchors <- anchors[invPerm(rank.anchors)] #sort anchors
    Gout_min <- Gout_max <- Gout[invPerm(rank.anchors)] #sort Gout by anchors
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
    # anchors <- sort(anchors) #into ascending order
    refined.cps <- anchors
    Gstar <- floor(Gout_min*3/4 + Gout_max/4)
    if(is.null(lambda)) lambda <- moseg.G[[1]]$lambda #| lambda %in% c("1se","min")
    else lambda <- lambda[1]
    limits <- c(0, anchors[-q] + floor(diff(anchors)/2), n)
    for (k in 1:q) {

      G_window <- #(anchors[k] - Gstar[k]):(anchors[k]-1)
        (anchors[k]- floor(Gout_min[k]/2) - Gstar[k]):(anchors[k]- floor(Gout_min[k]/2))
      G_window <- G_window[G_window>0]
      #if(any(G_window<1)) G_window <- 1: Gstar[k]
      lmod <- glmnet(X[G_window,], y[G_window], family = family, ...)

      G_window <-  #(anchors[k]):(anchors[k]+Gstar[k]-1)
        (anchors[k]+ floor(Gout_min[k]/2)):(anchors[k]+ floor(Gout_min[k]/2)+Gstar[k]-1)
      # if(any(G_window>n)) G_window <- (n-Gstar[k]+1):n
      G_window <- G_window[G_window<=n]
      rmod <- glmnet(X[G_window,], y[G_window], family = family, ...)

      rf <- refinement(X, y, anchors[k], Gstar[k], lambda, lmod, rmod,
                       L_min = limits[k], U_max = limits[k+1],
                       family=family)
      refined.cps[k] <- rf$cp

      if(do.plot){
        plot.ts(rf$objective, ylab = "Q", xaxt = "n")
        axis(1, at= 1:(rf$U - rf$L + 1), labels=rf$L:rf$U ) #(anchors[k]-Gstar[k]+1):(anchors[k]+Gstar[k])
        abline(v = which.min(rf$objective), col = "purple")
      }

    }
    par(mfrow = c(1,1))
  }



  out <- list(anchors= anchors, refined.cps = refined.cps,  moseg.G =moseg.G)
  attr(out, "class") <- "moseg.ms.cv"
  if(do.plot){
    par(mfrow = c(Glen,1))
    plot(out)
    pl <- recordPlot()
  } else pl <- NULL
  out$plot <- pl
  return(out)
}

#' @title Plot the multiscale moseg.cv detector
#' @method plot moseg.ms.cv
#' @description Plotting method for S3 objects of class \code{moseg.ms.cv}.
#' @param x \code{moseg.ms} object
#' @param type plot type, one of \code{"cv","mosum"}
#' @param ... unused
#'
#' @return cv or detector plots
#' @export
plot.moseg.ms.cv <- function(x, type = c("cv","mosum"), ...){
  Glen <- length(x$moseg.G)
  par(mfrow = c(Glen,1))
  for (ii in 1:Glen) {
    plot(x$moseg.G[[ii]], type = type)
    if(ii == Glen) abline(v = x$refined.cps, col = "purple")
  }
  par(mfrow = c(1,1))
}
