## moseg

#' @title default bandwidth
#' @keywords internal
getG <- function(p, n){
  err <- .012
  coeffs <- c(-0.01538316,0.02009430,0.02669467)
  out <- exp(err/coeffs[1] - (coeffs[2]/coeffs[1]) * log(log(p))   - (coeffs[3]/coeffs[1]) * log(log(n)) )
  round(out)
}

#' @title Create grid
#' @keywords internal
indexLE <- function(n, G, kap=1){ #define grid
  a <- 1:n
  R <- max(1, floor(kap*G))
  b <- a[union(seq.int(1, n-G+1, R),n-G+1)]
  return(b)
}

#' @title Locate change points
#' @keywords internal
get_local_maxima <- function(Tn, D_n, G, nu = 1/4) { ## eta check location
  n <- length(Tn)
  cps <- c()
  nu <- min(nu, 1)
  window <- floor(nu*G)
  for(t in (G+1):(n-G+1)){
    if( Tn[t] > D_n & Tn[t] == max(Tn[(t - window):(t+window)]) ) cps <- append(cps,t)
  }
  return(cps)
}

#' Detect and estimate multiple change points in the sparse regression model
#'
#' @param X design matrix
#' @param y response vector
#' @param G integer bandwidth; default is chosen based on dimensions of \code{X}
#' @param lambda regularisation parameter; either a numeric, or one of \code{"min","1se"} (see \link[glmnet]{cv.glmnet})
#' @param threshold numeric test rejection threshold; see details for default
#' @param n.cps chosen number of change points to return; if specified, overrides \code{threshold}
#' @param grid.resolution controls number of subsamples to take
#' @param nu numeric localisation tuning parameter
#' @param family response type, one of \code{"gaussian","binomial","poisson"}
#' @param do.refinement Boolean - perform location refinement
#' @param do.plot Boolean - return plots
#' @param do.scale Boolean - centre and scale \code{X,y}
#' @param ncores number of parallel cores
#' @param ... optional arguments to \code{glmnet}
#' @return List containing
#' \itemize{
#'   \item{\code{mosum}}{ numeric vector of mosum detector}
#'   \item{\code{cps}}{ integer vector of estimated change points}
#'   \item{\code{refined.cps}}{ integer vector of refined change points}
#'   \item{\code{G}}{ input}
#'   \item{\code{lambda}}{ input}
#'   \item{\code{threshold}}{ input}
#'   \item{\code{grid.resolution}}{ input}
#'   \item{\code{family}}{ input}
#'   \item{\code{plots}}{ list of detector and refinement plots}
#'   \item{\code{model_list}}{ list of fitted \code{glmnet} models at each grid point}
#' }
#' @details The default threshold is chosen as a product of exponents of \code{lambda} and \code{G}.
#'  For a more accurate, but slightly slower, procedure, see \link[moseg]{moseg.cv}.
#' @import glmnet
#' @import doParallel
#' @export
#'
#' @examples
#' eqX <- eqdata[,-c(1,9)]
#' eq_mosum <- moseg(as.matrix(eqX), eqdata[,9], 120, grid.resolution = 1/10, ncores = 2)
moseg <- function(X, y, G = NULL, lambda = c("min","1se"), family = c("gaussian","binomial","poisson"), threshold = NULL, n.cps = NULL,
                    grid.resolution = 1/G, nu = 0.5, do.refinement = TRUE, do.plot = TRUE, do.scale = TRUE,
                    ncores = NULL, ...){
  n <- dim(X)[1]
  p <- dim(X)[2]
  X <- as.matrix(X)
  if(do.scale) {
    X <- scale(X)
    y <- scale(y)
  }
  if(is.null(G)) G <- getG(p, n)
  G <- round(G)
  ## validate inputs
  if ( !(G > 0)) stop("G must be positive")
  if ( !(G <= n/2)) stop("G must be at most n/2")
  if(is.character(lambda))  lambda <-  match.arg(lambda, c("min","1se"))
  family <-  match.arg(family, c("gaussian","binomial","poisson"))
  if ( !(family %in% c("gaussian","binomial","poisson"))) stop("family must be \"gaussian\", \"binomial\", or \"poisson\" ")
  if(!is.null(grid.resolution) ) {
    if ( !(grid.resolution <= 1 & grid.resolution>= 1/G) ) stop("grid.resolution must be between 1/G and 1")
  } else grid.resolution <- 1/G
  if ( !(nu > 0)) stop("nu must be positive")

  ##parallel
  if(is.null(ncores)) ncores <- max(1, detectCores() -1)
  cl <- makeCluster(ncores)
  registerDoParallel(cl)

  pl <- list() #plots
  model_list <- list() #models
  coeffs <- matrix(0, nrow = n-G+1, ncol = p+1) #coefficients
  mosum <- rep(0, n) #detector
  ## Threshold
  if(is.null(lambda) || is.character(lambda) ) {
    model_list[[1]] <- cv.glmnet(X[1:G,], y[1:G], family = family, ...)
    model_list[[n-G+1]] <-  cv.glmnet(X[(n-G+1):n,], y[(n-G+1):n], family = family, ... )
    if(lambda == "1se") lambda <- (model_list[[1]]$lambda.1se + model_list[[n-G+1]]$lambda.1se  )/2
    if(lambda == "min") lambda <- (model_list[[1]]$lambda.min + model_list[[n-G+1]]$lambda.min  )/2
  } else {
    model_list[[1]] <- model_list[[n-G+1]] <- 0
  }
  ## Step 1
  T_Grid <- indexLE(n, G, kap=grid.resolution)
  if(is.list(model_list[[1]])) T_Grid_trim <- T_Grid[-c(1,length(T_Grid))] else T_Grid_trim <- T_Grid
  model_list_unsorted <- foreach (t = T_Grid_trim, .packages = "glmnet") %dopar%{
    G_window <- t:(t+G-1)
    glmnet(X[G_window,], y[G_window], family = family, ...)
  }
  for (ii in 1:length(T_Grid_trim) ) {
    t <- T_Grid_trim[ii]
    model_list[[t]] <- model_list_unsorted[[ii]]
  }
  for (ii in 1:length(T_Grid) ) {
    t <- T_Grid[ii]
    coeffs[t,] <- as.vector(coef.glmnet( model_list[[t]], lambda))
  }
  for (t in T_Grid[T_Grid<=n-2*G+1])  {
    mosum[G+t] <- sqrt(G/2) * norm(coeffs[t+G, ] - coeffs[t,], "2")
  }
  ## Locate
  cps <- c()
  if(!is.null(n.cps)) threshold <- 0 else
  if(is.null(threshold)) threshold <- max(mosum[c(G+1,n-G+1)]) *prod(c(G,lambda) ^c(0.5436,0.1113 ) )  #prod(c(G,lambda) ^c(0.2382334,-0.6929399 ) ) #
  cps <- get_local_maxima(mosum, threshold, G, nu= max(nu, grid.resolution/G) )

  if( is.null(cps) ){
    q <- 0
    refined.cps <- c()
  } else {
    ## Step 2
    q <- length(cps)
    refined.cps <- cps
    limits <- c(0, cps[-q] + floor(diff(cps)/2), n)
    if(do.plot) par(mfrow = c(max(q%/%6,1) , max(q%%6,1) ))
    if(do.refinement) for (k in 1:q) {
      L_ind <- max(1,cps[k]-G- floor(G/2))
      if(!is.list(model_list[[L_ind]]) ){ # - floor(G/2)
        G_window <- (L_ind):(cps[k]-1- floor(G/2))
        model_list[[L_ind]] <- glmnet(X[G_window,], y[G_window], family = family, ...)
      }
      R_ind <- min(n-G+1,cps[k]+ floor(G/2))
      if(!is.list(model_list[[R_ind]]) ){#+ floor(G/2)
        G_window <- (R_ind):(R_ind+G)
        model_list[[R_ind]] <- glmnet(X[G_window,], y[G_window], family = family, ...)
      }

      rf <- refinement(X, y, cps[k], G, lambda, model_list[[L_ind]],  model_list[[R_ind ]],
                       L_min = limits[k], U_max = limits[k+1],
                       family=family)
      refined.cps[k] <- rf$cp

      if(do.plot){
        plot.ts(rf$objective, ylab = "Q", xaxt = "n")
        axis(1, at= 1:(rf$U - rf$L + 1), labels=rf$L:rf$U )
        abline(v = which.min(rf$objective), col = "purple")
      }

    }
    if(do.plot) {
      par(mfrow = c(1,1))
      pl$refinement <- recordPlot();
    }
  }
  if(!is.null(n.cps)){ # select change point number
    ms.rank <- rank(-mosum[cps])
    cps <- cps[ms.rank %in% 1:n.cps]
    refined.cps <- refined.cps[ms.rank %in% 1:n.cps]
    q <- n.cps
  }
  stopCluster(cl)
  out <- list(mosum = mosum, cps = cps, refined.cps = refined.cps,
              G=G, lambda = lambda, threshold = threshold,
              grid.resolution = grid.resolution, family = family,
              plots = pl, model_list = model_list)
  attr(out, "class") <- "moseg"
  if(do.plot) {
    plot.moseg(out)
    out$plots$mosum <- recordPlot()
  }
  return(out)
}


#' @title Plot the moseg detector
#' @method plot moseg
#' @description Plotting method for S3 objects of class \code{moseg}.
#' @param x \code{moseg} object
#' @param ... unused
#'
#' @return A detector plot
#' @export
plot.moseg <- function(x, ...){
  if(x$grid.resolution == 1/x$G) plot.ts(x$mosum, ylab="Detector", xlab = "Time") # plot test statistic
  else plot(which(x$mosum>0), x$mosum[which(x$mosum>0)], ylab="Detector", xlab = "Time", ylim = c(0, max(x$mosum)), xlim = c(1,length(x$mosum)) )
  abline(h = x$threshold, col = "blue") #add threshold
  if(length(x$cps)>0) { # add estimated cps
    abline(v = x$cps, col = "red")
    abline(v = x$refined.cps , col = "purple")
  }
}




#' Fit a piecewise sparse regression model and evaluates a penalty. Fits a model to each stationary segment using \code{glmnet::cv.glmnet}.
#'
#' @param X design matrix
#' @param y response vector
#' @param cps vector of change point locations
#' @param lambda regularisation parameter; either a numeric, or one of \code{"min","1se"} (see \link[glmnet]{cv.glmnet})
#' @param family response type, one of \code{"gaussian","binomial","poisson"}
#' @param type type of prediction; see \code{?glmnet.predict}
#' @param do.plot Boolean - return coefficent plot
#' @param do.scale Boolean - scale \code{X,y}
#' @param ... optional arguments to \code{glmnet}
#'
#' @return list containing
#' \itemize{
#'   \item{\code{model}}{ list of fitted models}
#'   \item{\code{likelihood}}{ likelihood value}
#'   \item{\code{preds}}{ vector of fitted values for each time step}
#'   \item{\code{coeffs}}{ coefficient matrix}
#'   \item{\code{plot}}{ coefficient heatmap}
#' }
#' @import glmnet
#' @import doParallel
#' @export
#'
#' @examples

#' eqX <- eqdata[,-c(1,9)]
#' eq_mosum <- moseg.fit(as.matrix(eqX), as.matrix(eqdata[,9]), c(500))
moseg.fit <- function(X, y, cps, lambda = c("min","1se"), family =  c("gaussian","binomial","poisson"),
                        type = c("link","response", "coefficients", "nonzero", "class"), do.plot = TRUE, do.scale = TRUE, ...){
  cps <- sort(cps)
  family <- match.arg(family, c("gaussian","binomial","poisson"))
  type <- match.arg(type, c("link","response", "coefficients", "nonzero", "class"))
  n <- nrow(X)
  p <- ncol(X)
  X <- as.matrix(X)
  if(do.scale) {
    X <- scale(X)
    y <- scale(y)
  }
  if(is.character(lambda))  lambda <-  match.arg(lambda, c("min","1se"))
  starts <- c(0, cps); ends <- c(cps, n)
  q <- length(cps)
  out <- as.list(1:(q+1) )
  preds <- c()
  coeffs <- matrix(0, p+1,n)
  if(!is.null(colnames(X))){
    rownames(coeffs)[1:p+1] <- colnames(X)
    rownames(coeffs)[1] <- "intercept"
  }
  if(is.character(lambda)){
    for (ii in 1:(q+1)) {
      out[[ii]] <- cv.glmnet(X[(starts[ii]+1):ends[ii],], y[(starts[ii]+1):ends[ii]],
                                     nfolds = 10,family=family, ... )
      if(lambda == "1se") lambdaii <- out[[ii]]$lambda.1se
      if(lambda == "min") lambdaii <- out[[ii]]$lambda.min
      preds <- c(preds, predict(out[[ii]], X[(starts[ii]+1):ends[ii],], lambdaii, type=type) )
      if(do.plot) coeffs[,(starts[ii]+1):ends[ii]] <- as.matrix(coef(out[[ii]], s=lambdaii))
    }
  }
  if(is.numeric(lambda) & length(lambda) == 1){
    for (ii in 1:(q+1)) {
      out[[ii]] <- glmnet(X[(starts[ii]+1):ends[ii],], y[(starts[ii]+1):ends[ii]],family=family, lambda = lambda, ... )
      preds <- c(preds, predict(out[[ii]], X[(starts[ii]+1):ends[ii],], lambda, type=type) )
      if(do.plot) coeffs[,(starts[ii]+1):ends[ii]] <- as.matrix(coef(out[[ii]], s=lambda))
    }
  }
  RSS <- sum(likelihood(y, preds, family))

  if(do.plot){
    image(t(coeffs), axes=F);
    axis(1, at = (1:n -1) / (n-1), labels = 1:n)#;
    axis(2, at = (0:p ) / p, labels = rownames(coeffs), las = 2, cex.axis = .5 )
    title(xlab="Time")
    pl <- recordPlot()
  } else pl <- c()
  return(list(model = out, likelihood = RSS, preds = preds, coeffs = t(coeffs), plot = pl))
}

#' @title logarithmic factorial of `n`
#' @keywords internal
log.factorial <- function(n)  sum(log(1:max(n,1) ))



#' @title Add residuals to \code{glmnet} object
#' @keywords internal
add_residuals <- function(mod, X, y, s = "lambda.min"){
  mod$residuals <- predict(mod, X, s, type = "response") - y
  return(mod)
}




