# upsamplefill.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Original Octave code:
# Copyright (C) 2013 Juan Pablo Carbajal <carbajal@ifi.uzh.ch>
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
# 20201128  GvB       setup for gsignal v0.1.0
# 20220328  GvB       copy dimnames of x to output object
#------------------------------------------------------------------------------

#' Upsample and Fill
#'
#' Upsample and fill with given values or copies of the vector elements.
#'
#' @param x input data, specified as a numeric vector or matrix. In case of a
#'   vector it represents a single signal; in case of a matrix each column is a
#'   signal.
#' @param v vector of values to be placed between the elements of \code{x}.
#' @param copy logical. If TRUE then \code{v} should be a scalar
#'   (\code{length(v) == 1)} and each value in \code{x} are repeated \code{v}
#'   times. If FALSE (default), the values in the vector \code{v} are placed
#'   between the elements of \code{x}.
#'
#' @return upsampled vector or matrix
#'
#' @examples
#' u <- upsamplefill(diag(2), 2, TRUE)
#' u <- upsamplefill(diag(2), rep(-1, 3))
#'
#' @seealso \code{\link{upsample}}
#'
#' @author Juan Pablo Carbajal, \email{carbajal@@ifi.uzh.ch}.\cr
#' Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#
#' @export

upsamplefill <- function(x, v, copy = FALSE) {

  if (!is.numeric(x)) {
    stop("x must be a numeric vector or matrix")
  }

  if (is.vector(x)) {
    x <- matrix(x, ncol = 1)
    vec <- TRUE
  } else if (is.matrix(x)) {
    vec <- FALSE
  } else {
    stop("x must be a numeric vector or matrix")
  }
  nc <- ncol(x)
  nr <- nrow(x)

  if (!is.numeric(v) || !is.vector(v)) {
    stop("v must be a numeric vector")
  }
  if (!is.logical(copy)) {
    stop("copy must be a logical value TRUE or FALSE")
  }

  if (copy) {
    v <- v[1]
    if (v < 0) {
      stop("v must be a scalar value >= 0")
    }
    y <- pracma::kron(x, rep(1, v + 1))
  } else {
    n <- length(v) + 1
    N <- n * nr
    if (any(c(nr, nc) == 1)) {
      N        <- N * nc
      idx      <- seq(1, N, n)
      idx_c    <- setdiff(seq(1, N), seq(1, N, n))
      y        <- rep(0, N)
      y[idx]   <- x
      y[idx_c] <- pracma::repmat(v, 1, max(nr, nc))
    } else {
      idx        <- seq(1, N, n)
      idx_c      <- setdiff(seq(1, N), seq(1, N, n))
      y          <- matrix(0, N, nc)
      y[idx, ]   <- x
      y[idx_c, ] <- pracma::repmat(v, nr, nc)
    }
  }

  if (vec) {
    y <- as.vector(y)
  }
  dimnames(y) <- dimnames(x)
  y
}
