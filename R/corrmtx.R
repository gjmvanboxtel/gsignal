# corrmtx.R
# Copyright (C) 2026 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Octave version Copyright (C) 2006 Tang Chonghao <chadholton@qq.com>
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
# 20260903  GvB       setup for gsignal v0.4.0
#------------------------------------------------------------------------------

#' Data matrix for autocorrelation matrix estimation
#'
#' Given a vector \code{x} of length \code{N} and a model order \code{m},
#' compute the rectangular Toeplitz matrix \code{H} such that \code{H’*H} is
#' a biased estimate of the autocorrelation matrix. The size of \code{H}
#' depends on the selected method.
#'
#' @param x Input data, specified as a vector.
#' @param m Prediction model order, specified as a positive real integer.
#' @param method haracter string specifying matrix computation method,
#'   one of:
#'    \describe{
#'     \item{\code{autocorrelation}}{Uses both prewindowed and postwindowed
#'     data (default). The matrix can be used to perform autoregressive 
#'     parameter estimation using the Yule-Walker method. For more details,
#'     see \code{\link{aryule}}}
#'     \item{\code{prewindowed}}{Uses prewindowed data only}
#'     \item{\code{postwindowed}}{Uses postwindowed data only}
#'     \item{\code{covariance}}{Uses nonwindowed data}
#'     \item{\code{modified}}{Uses forward and backward prediction error 
#'     estimates}
#'  }
#'
#' @return A list containing the following elements:
#'   \describe{
#'    \item{H}{Data matrix, returned for autocorrelation matrix estimation. 
#'    The size of H depends on the matrix computation method specified in
#'    \code{method}}
#'    \item{R}{Biased autocorrelation matrix estimate \code{H’*H}}
#'   }
#'
#' @examples
#' x <- 1:5
#' m <- 2
#' HR1 <- corrmtx (x, m)
#' HR2 <- corrmtx (x, m, 'cov')
#'
#' @seealso \code{\link{aryule}}
#'
#' @author Tang Chonghao \email{chadholton@@qq.com}.\cr Conversion to R by Geert
#'   van Boxtel \email{G.J.M.vanBoxtel@@gmail.com}.
#
#' @export

corrmtx <- function(x, m, method = c('autocorrelation',
                                     'prewindowed', 'postwindowed',
                                     'covariance', 'modified')) {

  x <- as.vector(x)
  if (!isPosscal(m) || !isWhole(m)) {
    stop("m must be a positive integer.")
  }
  method <- match.arg(method)

  n <- length(x)
  if (n <= m) {
    stop ("length of x must be greater than m.")
  }

  l <- n - m 
  if (method == 'autocorrelation') {
    H <- pracma::Toeplitz(c(x, rep(0L, m)), c(x[1], rep(0L, m))) / sqrt(n)
  }
  else if (method == 'prewindowed') {
    H <- pracma::Toeplitz(c(x, rep(0L, m)), c(x[1], rep(0L,m)))
    H <- H[1:n,] / sqrt(n)
  }
  else if (method == 'postwindowed') {
    H <- pracma::Toeplitz(c(x[(m + 1):n], rep(0L, m)),
                          c(x[m + 1], rev(x[1:m]))) / sqrt(n)
  }
  else if (method == 'covariance') {
    H <- pracma::Toeplitz(x[(m + 1):n], c(x[m + 1], rev(x[1:m]))) / sqrt(l)
  }
  else if (method == 'modified') {
      H_f <- pracma::Toeplitz(x[(m + 1):n], c(x[m + 1], rev(x[1:m])))
      H_b <- pracma::hankel(x[1:l], x[l:n])
      H <- rbind(H_f, Conj(H_b)) / sqrt(2 * l)
  }

  R <- t(H) %*% H
  list(H = H, R = R)
}
