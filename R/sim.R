## data simulation

#' Simulate from a sparse piecewise regression model
#'
#' @param n sample size
#' @param p number of parameters
#' @param sparsity number of non-zero parameters
#' @param q number of change points
#' @param sigma.noise error standard deviation
#' @param sigma.x covariance structure of \code{X}, one of \code{"id", "band", "ar"}
#' @param kappa change size
#'
#' @return List containing
#' \itemize{
#'   \item{\code{y}}{ response vector}
#'   \item{\code{X}}{ matrix of covariates}
#'   \item{\code{cps}}{ vector of change points}
#'   \item{\code{beta}}{ matrix of parameters corresponding to each regime}
#' }
#' and all inputs
#' @export
#'
#' @examples
#' moseg.sim(100, 50)
moseg.sim <- function(n, p, sparsity = floor(sqrt(p)), q = 1, sigma.noise = 1, sigma.x = c("id","band","ar"), kappa = 1){
  sigma.x <- match.arg(sigma.x, c("id","band","ar"))
  if(sigma.x=="id"){
     Sigma.x <- diag(1,p)
  } else if(sigma.x=="band"){
    Sigma.x <- stats::toeplitz(c(1, .6, .3, rep(0, p-3)))
  } else if(sigma.x=="ar"){
    Sigma.x <- stats::toeplitz(.6^(1:p - 1))
  }
  sparsity <- min(sparsity,p)
  beta0 <- c( rep_len(c(1,-1), sparsity), rep(0, p-sparsity))
  kappa.factor <-  kappa/norm(2*beta0,"2")
  cps <- floor(n/(q+1) * 0:(q+1) )
  X <- MASS::mvrnorm(n, rep(0,p), Sigma.x) #MASS

  y <- rnorm(n, 0, sigma.noise)
  beta <- c()
  for (kk in 1:(q+1)) {
    t.range <- (cps[kk]+1):cps[kk+1]
    y[t.range] <- y[t.range] + X[t.range,] %*% beta0 * kappa.factor * (-1)^kk
    beta <- cbind(beta, beta0 * kappa.factor * (-1)^kk)
  }

  out <- list(y=y, X=X, cps = cps, beta = beta,
    n=n, p=p, sparsity=sparsity, q=q, sigma.noise = sigma.noise, sigma.x = sigma.x, kappa = kappa)
  return(out)
}
