# levinson.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Original Octave code:
# Copyright (C) 1999 Paul Kienzle <pkienzle@users.sf.net>
# Copyright (C) 2006 Peter V. Lanspeary <peter.lanspeary@.adelaide.edu.au>
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
# 20201105  GvB       setup for gsignal v0.1.0
#------------------------------------------------------------------------------

#' Durbin-Levinson Recursion
#'
#' Use the Durbin-Levinson algorithm to compute the coefficients of an
#' autoregressive linear process.
#'
#' \code{levinson} uses the Durbin-Levinson algorithm to solve:
#' \deqn{toeplitz(acf(1:p)) * x = -acf(2:p+1)} The solution \code{c(1, x)} is
#' the denominator of an all pole filter approximation to the signal \code{x}
#' which generated the autocorrelation function acf.
#'
#' From ref [2]: Levinson recursion or Levinson–Durbin recursion is a procedure
#' in linear algebra to recursively calculate the solution to an equation
#' involving a Toeplitz matrix. Other methods to process data include Schur
#' decomposition and Cholesky decomposition. In comparison to these, Levinson
#' recursion (particularly split Levinson recursion) tends to be faster
#' computationally, but more sensitive to computational inaccuracies like
#' round-off errors.
#'
#' @param acf autocorrelation function for lags 0 to \code{p}, specified as a
#'   vector or matrix. If r is a matrix, the function finds the coefficients for
#'   each column of \code{acf} and returns them in the rows of \code{a}.
#' @param p model order, specified as a positive integer. Default:
#'   \code{NROW(acf) - 1}.
#'
#' @return A \code{list} containing the following elements:
#'   \describe{
#'     \item{a}{vector or matrix containing \code{(p+1)} autoregression
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
#' @examples
#' ## Estimate the coefficients of an autoregressive process given by
#' ## x(n) = 0.1x(n-1) - 0.8x(n-2) - 0.27x(n-3) + w(n).
#' a <- c(1, 0.1, -0.8, -0.27)
#' v <- 0.4
#' w <- sqrt(v) * rnorm(15000)
#' x <- filter(1, a, w)
#' xc <- xcorr(x, scale = 'biased')
#' acf <- xc$R[-which(xc$lags < 0)]
#' lev <- levinson(acf, length(a) - 1)
#'
#' @author Paul Kienzle, \email{pkienzle@@users.sf.net},\cr
#'  Peter V. Lanspeary, \email{pvl@@mecheng.adelaide.edu.au}.\cr
#'  Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#'
#' @references [1] Steven M. Kay and Stanley Lawrence Marple Jr. (1981).
#'   Spectrum analysis – a modern perspective. Proceedings of the IEEE, Vol 69,
#'   1380-1419.\cr
#'   [2] \url{https://en.wikipedia.org/wiki/Levinson_recursion}
#'
#' @export

levinson <- function(acf, p = NROW(acf)) {

  # check parameters
  if (!(is.vector(acf) || is.matrix(acf))) {
    stop("acf must be a vector or matrix")
  }

  if (is.vector(acf)) {
    vec <- TRUE
    acf <- as.matrix(acf, ncol = 1)
  } else {
    vec <- FALSE
  }
  nr <- nrow(acf)
  nc <- ncol(acf)
  if (nr < 2) {
    stop("acf must be a vector or matrix of length > 1")
  }

  if (!isScalar(p) || !isWhole(p) || !is.numeric(p) || p <= 0.5) {
    stop("p must be a positive integer > 0")
  }
  # end of parameter checking

  aggr_a <- aggr_v <- aggr_k <- NULL
  for (icol in seq_len(nc)) {
    ref <- rep(0L, p)
    g <- -acf[2] / acf[1]
    a <- g
    v <- Re((1 - g * Conj(g)) * acf[1])
    ref[1] <- g
    if (p > 1) {
      for (t in 2:p) {
        g <- as.vector(- (acf[t + 1] + a %*% acf[seq(t, 2, -1)]) / v)
        a <- c((a + g * Conj(a[seq(t - 1, 1, -1)])), g)
        v <- v * (1 - Re(g * Conj(g)))
        ref[t] <- g
      }
    }
    aggr_a <- rbind(aggr_a, c(1, a))
    aggr_v <- c(aggr_v, v)
    aggr_k <- rbind(aggr_k, ref)
  }

  if (vec) {
    rv <- list(a = as.vector(aggr_a), e = aggr_v, k = as.vector(aggr_k))
  } else {
    rv <- list(a = aggr_a, e = aggr_v, k = t(aggr_k))
  }
  rv
}
