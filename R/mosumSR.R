## mosumsr

#' @title Create grid
#' @keywords internal
indexLE <- function(n, G, kap=1){ #define grid
  a <- 1:n
  R <- floor(kap*G)
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
#' @param G integer bandwidth; defaults to \code{round(30 + ncol(X)/100)} 
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
#'   \item{\code{refined_cps}}{ integer vector of refined change points}
#'   \item{\code{threshold}}
#'   \item{\code{lambda}}
#'   \item{\code{model_list}}{ list of fitted \code{glmnet} models at each grid point}
#'   \item{\code{plots}}{ list of detector and refinement plots}
#'   \item{\code{family}}{ input}
#' }
#' @details The default threshold is chosen as a product of exponents of \code{lambda} and \code{G}.
#'  For a more accurate, but slightly slower, procedure, see \link[mosumsr]{mosumsr.cv}.
#' @export
#'
#' @examples
#' eqX <- eqdata[,-c(1,9)]
#' eq_mosum <- mosumsr(as.matrix(eqX), eqdata[,9], 120, grid.resolution = 1/10, ncores = 2)
mosumsr <- function(X, y, G = NULL, lambda = c("min","1se"), family = c("gaussian","binomial","poisson"), threshold = NULL, n.cps = NULL,
                    grid.resolution = 1/G, nu = 0.5, do.refinement = TRUE, do.plot = TRUE, do.scale = TRUE,
                    ncores = NULL, ...){
  n <- dim(X)[1]
  p <- dim(X)[2]
  X <- as.matrix(X)
  if(do.scale) {
    X <- scale(X)
    y <- scale(y)
  }
  if(is.null(G)) G <- 30 + p/100
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
  if(is.null(threshold)) threshold <- prod(c(G,lambda) ^c(0.2382334,-0.6929399 ) ) #
  cps <- get_local_maxima(mosum, threshold, G, nu= max(nu, grid.resolution/G) )

  if( is.null(cps) ){
    q <- 0
    refined_cps <- c()
  } else {
    ## Step 2
    q <- length(cps)
    refined_cps <- cps
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

      rf <- refinement(X, y, cps[k], G, lambda, model_list[[L_ind]],  model_list[[R_ind ]], family=family)
      refined_cps[k] <- rf$cp

      if(do.plot){
        plot.ts(rf$objective, ylab = "Q", xaxt = "n")
        axis(1, at=1:(2*G), labels= (cps[k]-G+1):(cps[k]+G))
        abline(v = which.min(rf$objective), col = "purple")
      }

    }
    if(do.plot) {
      pl$refinement <- recordPlot();
      par(mfrow = c(1,1))
    }
  }
  if(!is.null(n.cps)){ # select change point number
    ms.rank <- rank(-mosum[cps])
    cps <- cps[ms.rank %in% 1:n.cps]
    refined_cps <- refined_cps[ms.rank %in% 1:n.cps]
    q <- n.cps
  }
  ## Plot
  if(do.plot){
    if(grid.resolution == 1/G) plot.ts(mosum, ylab="Detector", xlab = "Time") # plot test statistic
    else plot(which(mosum>0), mosum[which(mosum>0)], ylab="Detector", xlab = "Time", ylim = c(0, max(mosum)), xlim = c(1,n) )
    abline(h = threshold, col = "blue") #add threshold
    if(q>0) { # add estimated cps
      par(mfrow =c(1,1))
      abline(v = cps, col = "red")
      abline(v = refined_cps , col = "purple")
    }
    pl$mosum <- recordPlot()
  }
  stopCluster(cl)
  out <- list(mosum = mosum, cps = cps, refined_cps = refined_cps,
              threshold = threshold, lambda = lambda, model_list = model_list, #
              plots = pl, family = family)
  attr(out, "class") <- "mosumsr"
  return(out)
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
#' @export
#'
#' @examples

#' eqX <- eqdata[,-c(1,9)]
#' eq_mosum <- mosumsr.fit(as.matrix(eqX), as.matrix(eqdata[,9]), c(500))
mosumsr.fit <- function(X, y, cps, lambda = c("min","1se"), family =  c("gaussian","binomial","poisson"),
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
      out[[ii]] <- glmnet::cv.glmnet(X[(starts[ii]+1):ends[ii],], y[(starts[ii]+1):ends[ii]],
                                     nfolds = 10,family=family, ... )
      if(lambda == "1se") lambdaii <- out[[ii]]$lambda.1se
      if(lambda == "min") lambdaii <- out[[ii]]$lambda.min
      preds <- c(preds, predict(out[[ii]], X[(starts[ii]+1):ends[ii],], lambdaii, type=type) )
      if(do.plot) coeffs[,(starts[ii]+1):ends[ii]] <- as.matrix(coef(out[[ii]], s=lambdaii))
    }
  }
  if(is.numeric(lambda) & length(lambda) == 1){
    for (ii in 1:(q+1)) {
      out[[ii]] <- glmnet::glmnet(X[(starts[ii]+1):ends[ii],], y[(starts[ii]+1):ends[ii]],family=family, lambda = lambda, ... )
      preds <- c(preds, predict(out[[ii]], X[(starts[ii]+1):ends[ii],], lambda, type=type) )
      if(do.plot) coeffs[,(starts[ii]+1):ends[ii]] <- as.matrix(coef(out[[ii]], s=lambda))
    }
  }
  RSS <- sum(likelihood(y, preds, family))

  if(do.plot){
    image(t(coeffs), axes=F);
    axis(1, at = (1:n -1) / (n-1), labels = 1:n)#;
    axis(2, at = (0:p ) / p, labels = rownames(coeffs) )
    title(xlab="Time")
    pl <- recordPlot()
  } else pl <- c()
  return(list(model = out, likelihood = RSS, preds = preds, coeffs = t(coeffs), plot = pl))
}

#' @title logarithmic factorial of `n`
#' @keywords internal
log.factorial <- function(n)  sum(log(1:max(n,1) ))





# mosumsr.IC <- function(X, y, G, lambda = NULL, max.cps = NULL, family = c("gaussian","binomial","poisson"),
#                        path.length = 5, penalty = NULL, nu = 0.25, do.plot = TRUE,
#                        ncores = NULL, ...){
#   n <- nrow(X)
#   family <- match.arg(family, c("gaussian","binomial","poisson"))
#   if(is.null(lambda)){
#     lambda.max <- max( abs(t(y - mean(y)*(1-mean(y)) ) %*% X ) )/n #/2
#     lambda <- round(exp(seq(log(lambda.max), log(lambda.max * .0001), length.out = path.length)), digits = 10)
#   }
#   path.length <- length(lambda)
#   if(!is.null(max.cps)) out_IC <- matrix(Inf,path.length, max.cps+1)
#   for (jj in 1:path.length) {
#     ms <- mosumsr(X,y,G,lambda[jj],family,threshold = 0, nu = nu, do.plot = do.plot, ncores = ncores, ...)
#     cps <- ms$refined_cps
#     max.cps.jj <- min(max.cps, length(cps))
#     if(jj==1 & is.null(max.cps)) {max.cps <- max.cps.jj; out_IC <- matrix(Inf,path.length, max.cps+1)}
#     ranks <- rank(ms$mosum[cps])
#     out_list <- list()
#
#     for (ii in 0:max.cps) {
#       ii_cps <- cps[ranks<=ii]
#       if(min(diff(ii_cps)) > 10  ){
#         ii_fit <- mosumsr.fit(X, y, ii_cps, family = family, penalty = penalty, lambda = lambda[jj], ...)
#         out_IC[jj,ii+1] <- ii_fit$IC
#         ii_list <- list(lambda = lambda[jj], cps = ii_cps, threshold = ms$mosum[cps[ranks==ii]], fit = ii_fit)
#         out_list[[ (jj-1)*path.length +  ii+1 ]] <- ii_list
#         if (min(out_IC) == ii_fit$IC) min_list <- ii_list
#       } else warning("Segment too short. Returning Inf")
#
#     }
#   }
#   if(do.plot) {
#     min.point <- which(out_IC==min(out_IC), arr.ind=TRUE)
#     matplot(0:max.cps, t(out_IC), type = 'b', col = 1+1:path.length, pch = 1+1:path.length,
#             xlab = 'q', ylab = 'IC', main = 'IC for change point number estimation')
#     abline(v = (0:max.cps)[min.point[2]])
#     legend('topleft', legend = lambda, col = 1+1:path.length, pch = 1+1:path.length, lty = 1)
#
#     # matplot(lambda, out_IC, type = 'b', col = 2:(max.cps + 2), pch = 2:(max.cps + 2),
#     #         log = 'x', xlab = 'lambda (log scale)', ylab = 'IC', main = 'IC for change point number estimation')
#     # abline(v = lambda[min.point[1]])
#     # legend('topleft', legend = 0:max.cps, col = 2:(max.cps + 2), pch = 2:(max.cps + 2), lty = 1)
#   }
#   min_list$IC <- out_IC
#   min_list$models <- out_list
#   return(min_list)
# }




#' @title Add residuals to \code{glmnet} object
#' @keywords internal
add_residuals <- function(mod, X, y, s = "lambda.min"){
  mod$residuals <- predict(mod, X, s, type = "response") - y
  return(mod)
}




