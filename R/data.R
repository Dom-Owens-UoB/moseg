#' Equity premium and macro/financial variables
#'
#' A dataset containing time series of the equity premium and other macro/financial variables
#'
#' @format A data frame with 1032 rows and 16 variables:
#' \describe{
#'   \item{yyyymm}{Date}
#'   \item{b.m}{Book-to-Market Ratio}
#'   \item{tbl}{Treasury Bills}
#'   \item{lty}{Long Term Yield}
#'   \item{ntis}{Net Equity Expansion}
#'   \item{infl}{Inflation}
#'   \item{ltr}{Long Term Rate of Returns}
#'   \item{svar}{Stock Variance}
#'   \item{equity_premium}{Equity Premium}
#'   \item{dp}{Dividend Price Ratio}
#'   \item{dy}{Dividend Yield}
#'   \item{ep}{Earnings Price Ratio}
#'   \item{de}{Dividend Payout Ratio}
#'   \item{tms}{Term Spread}
#'   \item{dfy}{Default Yield Spread}
#'   \item{dfr}{Default Return Spread}
#' }
#' @source \url{https://sites.google.com/view/agoyal145/?redirpath=/}
#' @references Welch, Ivo, and Amit Goyal.
#' "A comprehensive look at the empirical performance of equity premium prediction."
#' The Review of Financial Studies 21.4 (2008): 1455-1508.
#' @export
"eqdata"
