# xcov.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Octave version:
# Copyright (C) 1999-2001 Paul Kienzle <pkienzle@users.sf.net>
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
# 2020313  GvB       setup for gsignal v0.1.0
#------------------------------------------------------------------------------

#' Cross-covariance
#'
#' Compute covariance at various lags (= correlation(x-mean(x), y-mean(y))).
#'
#' @param x Input, numeric or complex vector or matrix. Must not be missing.
#' @param y Input, numeric or complex vector data.  If \code{x} is a matrix (not
#'   a vector), \code{y} must be omitted. \code{y} may be omitted if \code{x} is
#'   a vector; in this case \code{xcov} estimates the autocovariance of
#'   \code{x}.
#' @param maxlag Integer scalar. Maximum covariance lag. If omitted, the
#'   default value is \code{N-1}, where \code{N} is the greater of the lengths
#'   of \code{x} and \code{y} or, if \code{x} is a matrix, the number of rows in
#'   \code{x}.
#' @param scale Character string. Specifies the type of scaling applied to the
#'   covariation vector (or matrix). matched to one of:
#'   \describe{
#'     \item{"none"}{return the unscaled covariance, C}
#'     \item{"biased"}{return the biased average, C/N}
#'     \item{"unbiased"}{return the unbiased average, C(k)/(N-|k|)}
#'     \item{"coeff"}{return C/(covariance at lag 0)},
#'     where \code{k} is the lag, and \code{N} is the length of \code{x}
#'   }
#'  If omitted, the default value is \code{"none"}. If \code{y} is supplied but
#'  does not have the same length as \code{x}, scale must be \code{"none"}.
#'
#' @return A list containing the following variables:
#' \describe{
#'   \item{C}{array of covariance estimates}
#'   \item{lags}{vector of covariance lags \code{[-maxlag:maxlag]}}
#' }
#' The array of covariance estimates has one of the following forms:
#' \enumerate{
#'   \item Cross-covariance estimate if X and Y are vectors.
#'   \item Autocovariance estimate if is a vector and Y is omitted.
#'   \item If \code{x} is a matrix, \code{C} is a matrix containing the
#'   cross-covariance estimates of each column with every other column. Lag
#'   varies with the first index so that \code{C} has \code{2 * maxlag + 1} rows
#'   and \eqn{P^2} columns where \code{P} is the number of columns in \code{x}.
#' }
#' @seealso \code{\link{xcorr}}.
#'
#' @examples
#' x <- rnorm(1000)
#' cl <- xcov(x, maxlag = 10, scale = 'coeff')
#' plot (cl$lags, cl$C, type = "h", xlab = "", ylab = "")
#' points (cl$lags, cl$C)
#' abline(h = 0)
#'
#' @author Paul Kienzle, \email{pkienzle@@users.sf.net}.\cr
#'   Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#'
#' @export

xcov <- function(x, y = NULL,
                 maxlag = if (is.matrix(x)) nrow(x) - 1
                 else max(length(x), length(y)) - 1,
                 scale = c("none", "biased", "unbiased", "coeff")) {
  if (is.null(y)) {
    ret <- xcorr(x - colMeans(as.matrix(x)), maxlag = maxlag, scale = scale)
  } else {
    ret <- xcorr(x - colMeans(as.matrix(x)), y - colMeans(as.matrix(y)),
                 maxlag = maxlag, scale = scale)
  }
  list(C = ret$R, lags = ret$lags)
}
