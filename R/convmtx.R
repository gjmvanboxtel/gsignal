# convmtx.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Octave version Copyright (C) 2003 David Bateman <adb014@gmail.com>
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
# 20200212  GvB       setup for gsignal v0.1.0
#------------------------------------------------------------------------------

#' Convolution matrix
#'
#' Returns the convolution matrix for a filter kernel.
#'
#' Computing a convolution using \code{conv} when the signals are vectors is
#' generally more efficient than using \code{convmtx}. For multichannel signals,
#' however, when  a large number of vectors are to be convolved with the same
#' filter kernel, \code{convmtx} might be more efficient.
#'
#' The code \code{cm <- convmtx(h, n)} computes the convolution matrix of the
#' filter kernel \code{h} with a vector of length \code{n}. Then, \code{cm %*%
#' x} gives the convolution of \code{h} and \code{x}.
#'
#' @param h Input, coerced to a vector, representing the filter kernel
#' @param n Length of vector(s) that \code{h} is to be convolved with.
#'
#' @return Convolution matrix of input \code{h} for a vector of length \code{n}.
#'   If \code{h} is a vector of length \code{m}, then the convolution matrix has
#'   \code{m + n - 1} rows and \code{n} columns.
#'
#' @examples
#' N <- 1000
#' a <- runif(N)
#' b <- runif(N)
#' cm <- convmtx(b, N)
#' d <- cm %*% a
#'
#' cref = conv(a, b)
#' all.equal(max(d - cref), 0)
#'
#' @seealso \code{\link{conv}}
#'
#' @author David Bateman \email{adb014@@gmail.com}.\cr Conversion to R by Geert
#'   van Boxtel \email{G.J.M.vanBoxtel@@gmail.com}.
#
#' @export

convmtx <- function(h, n) {

  h <- as.vector(h)
  if (!isPosscal(n) || !isWhole(n)) {
    stop("n must be a positive integer.")
  }

  y <- pracma::Toeplitz(c(h, rep(0L, (n - 1))), c(h[1], rep(0L, (n - 1))))
  y
}
