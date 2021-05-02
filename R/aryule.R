# aryule.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Original Octave function:
# Copyright (C) 1999 Paul Kienzle <pkienzle@users.sf.net>
# Copyright (C) 2006 Peter Lanspeary
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#
# Version history
# 20201106  GvB       setup for gsignal v0.1.0
#------------------------------------------------------------------------------

#' Autoregressive model coefficients - Yule-Walker method
#'
#' compute autoregressive all-pole model parameters using the Yule-Walker
#' method.
#'
#' \code{aryule} uses the Levinson-Durbin recursion on the biased estimate of
#' the sample autocorrelation sequence to compute the parameters.
#'
#' @param x input data, specified as a numeric or complex vector or matrix. In
#'   case of a vector it represents a single signal; in case of a matrix each
#'   column is a signal.
#' @param p model order; number of poles in the AR model or limit to the number
#'   of poles if a valid criterion is provided. Must be smaller than the length
#'   of \code{x} minus 1.
#'
#' @return A \code{list} containing the following elements:
#'   \describe{
#'     \item{a}{vector or matrix containing \code{(p + 1)} autoregression
#'     coefficients. If \code{x} is a matrix, then each row of a corresponds to
#'     a column of \code{x}. \code{a} has \code{p + 1} columns.}
#'     \item{e}{white noise input variance, returned as a vector. If \code{x} is
#'     a matrix, then each element of e corresponds to a column of \code{x}.}
#'     \item{k}{Reflection coefficients defining the lattice-filter embodiment
#'     of the model returned as vector or a matrix. If \code{x} is a matrix,
#'     then each column of \code{k} corresponds to a column of \code{x}.
#'     \code{k} has \code{p} rows.}
#'   }
#'
#' @note The power spectrum of the resulting filter can be plotted with
#'   \code{pyulear(x, p)}, or you can plot it directly with
#'   \code{ar_psd(a,v,...)}.
#'
#' @examples
#' a <- Arma(1, c(1, -2.7607, 3.8106, -2.6535, 0.9238))
#' y <- filter(a, rnorm(1024))
#' coefs <- aryule(y, 4)
#'
#' @author Paul Kienzle, \email{pkienzle@@users.sf.net},\cr
#'  Peter V. Lanspeary, \email{pvl@@mecheng.adelaide.edu.au}.\cr
#'  Conversion to R by Geert van Boxtel, \email{gjmvanboxtel@@gmail.com}.
#'
#' @seealso \code{\link{ar_psd}}, \code{\link{arburg}}
#'
#' @export

aryule <- function(x, p)  {

  # check parameters
  if (!(is.vector(x) || is.matrix(x)) || !is.numeric(x)) {
    stop("x must be a numeric or vector or matrix")
  }

  if (is.vector(x)) {
    vec <- TRUE
    x <- as.matrix(x, ncol = 1)
  } else {
    vec <- FALSE
  }
  nr <- nrow(x)
  nc <- ncol(x)

  if (!isScalar(p) || !isWhole(p) || !is.numeric(p) || p <= 0.5) {
    stop("p must be a positive integer")
  }
  if (p >= nr - 1) {
    stop(paste0("p must be less than the length of x (", nr, ") - 1"))
  }
  # end of parameter checking

  # loop over columns
  aggr_a <- aggr_e <- aggr_k <- NULL
  for (icol in seq_len(nc)) {

    xc <- xcorr(x[, icol], maxlag = p + 1, scale = "biased")
    R <- xc$R[-c(1:(p + 1))]    # remove negative autocorrelation lags
    R[1] <- Re(R[1])            # levinson/toeplitz requires
                                # exactly R[1]==Conj(R[1])
    lev <- levinson(R, p)
    aggr_a <- rbind(aggr_a, lev$a)
    aggr_e <- c(aggr_e, lev$e)
    aggr_k <- rbind(aggr_k, lev$k)
  }

  if (vec) {
    rv <- list(a = as.vector(aggr_a), e = aggr_e, k = as.vector(aggr_k))
  } else {
    rv <- list(a = aggr_a, e = aggr_e, k = t(aggr_k))
  }
  rv
}
