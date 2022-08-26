## refinement

#' Refine a change point location estimate
#'
#' @param X design matrix
#' @param y response vector
#' @param k integer change point location estimate
#' @param G integer bandwidth
#' @param lambda numeric regularisation parameter
#' @param l_mod \code{glmnet} object from left of change point
#' @param r_mod \code{glmnet} object from right of change point
#' @param family response type, one of \code{"gaussian","binomial","poisson"}
#'
#' @return list containing
#'  \itemize{
#'    \item{\code{cp}}{ refined change point estimate}
#'    \item{\code{objective}}{ Loss function}
#' }
#' @keywords internal
#'
refinement <- function(X, y, k, G, lambda, l_mod, r_mod, L_min = NULL, U_max = NULL,
                       family = c("gaussian","binomial","poisson")){
  family <-  match.arg(family, c("gaussian","binomial","poisson"))

  n <- nrow(X)
  L <- max(1,k-G+1)
  if(!is.null(L_min)) L <- max(L_min, L)
  U <- min(n,k+G)
  if(!is.null(U_max)) U <- min(U_max, U)
  yy <- y[L:U]
  XX <- X[L:U,]
  objective <- rep(0, U-L)

    l_preds <- predict(l_mod, XX, lambda, type = "response")
    l_res <-  likelihood(yy,l_preds,family)
    r_preds <- predict(r_mod, XX, lambda, type = "response")
    r_res <- likelihood(yy,r_preds,family)

  for (t in 1:(U-L) ) {
    objective[t] <- (sum(l_res[1:t])) + (sum(r_res[(t+1):(U-L+1)])) #norm(y[1:t] - l_pred[1:t], "2")  + norm(y[(t+1):n] - l_pred[(t+1):n] , "2")
  }
  cp <- L - 1 +  which.min(objective)
  out <- list(cp = cp, objective = objective, L=L, U=U)
  return(out)
}
#' @title Evaluate gaussian likelihood
#' @keywords internal
gauss_likelihood <- function(yy,preds) (yy - preds)^2

#' @title Evaluate binomial likelihood
#' @keywords internal
bin_likelihood <- function(yy,preds) -(yy*log(preds) + (1-yy) *log(1-preds))

#' @title Evaluate poisson likelihood
#' @keywords internal
pois_likelihood <- function(yy,preds) -yy*log(preds) + preds

#' @title Evaluate likelihood
#' @keywords internal
likelihood <- function(yy, preds, family = "gaussian"){
  family <- match.arg(family, c("gaussian","binomial","poisson"))
  if(family == "gaussian") out <- gauss_likelihood(yy,preds)
  if(family == "binomial") out <- bin_likelihood(yy,preds)
  if(family == "poisson") out <- pois_likelihood(yy,preds)
  #else paste0("family must be \"gaussian\", \"binomial\", or \"poisson\" ")
  return(out)
}








