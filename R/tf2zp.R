# tf2zp.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
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
# 20200331  GvB       setup for gsignal v0.1.0
# 20200401  GvB       catch k == 0
# 20200403  GvB       compute roots with the "eigen" method, and sort them
# 20200406  GvB       validated
# 20210326  GvB       renamed k to g, return object of class 'Zpg'
# 20210506  GvB       sort output z and p
#------------------------------------------------------------------------------

#' Transfer function to zero-pole-gain form
#'
#' Convert digital filter transfer function parameters to zero-pole-gain form.
#'
#' @param b moving average (MA) polynomial coefficients, specified as a numeric
#'   vector or matrix. In case of a matrix, then each row corresponds to an
#'   output of the system. The number of columns of \code{b} must be less than
#'   or equal to the length of \code{a}.
#' @param a autoregressive (AR) polynomial coefficients, specified as a vector.
#'
#'@return A list of class Zpg with the following list elements:
#' \describe{
#'   \item{z}{complex vector of the zeros of the model (roots of \code{B(z)})}
#'   \item{p}{complex vector of the poles of the model (roots of \code{A(z)})}
#'   \item{g}{overall gain (\code{B(Inf)})}
#' }
#'
#' @seealso \code{\link{filter}}
#'
#' @examples
#' b <- c(2, 3)
#' a <- c(1, 1/sqrt(2), 1/4)
#' zpk <- tf2zp(b, a)
#'
#' @author Geert van Boxtel, \email{gjmvanboxtel@@gmail.com}
#'
#' @export

tf2zp <- function(b, a) {

  if (!(is.vector(b) || is.matrix(b))) {
    stop("b must be a vector or a matrix")
  }
  if (!is.vector(a)) {
    stop("a must be a vector")
  }
  if (NCOL(b) > length(a)) {
    stop("The number of columns of b must be <= length(a)")
  }

  if (length(b) > 0) {
    z <- pracma::roots(b)
  } else {
    z <- NULL
  }

  if (length(a) > 0) {
    p <- pracma::roots(a)
  } else {
    p <- NULL
  }

  if (a[1] != 0) {
    k <- b[1] / a[1]
    if (k <= 0) k <- 1
  } else {
    k <- 1
  }

  Zpg(z = sort(z), p = sort(p), g = k)
}
