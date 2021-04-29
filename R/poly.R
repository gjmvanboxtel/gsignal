# poly.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Original Octave function Copyright (C) 1994-2017 John W. Eaton
# Author: KH <Kurt.Hornik@wu-wien.ac.at>
# Created: 24 December 1993
# Adapted-By: jwe
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
# 20200126  GvB       setup for gsignal v0.1.0
#------------------------------------------------------------------------------

#' Polynomial with specified roots
#'
#' Compute the coefficients of a polynomial when the roots are given, or the
#' characteristic polynomial of a matrix.
#'
#' If a vector is passed as an argument, then \code{poly((x)} is a vector of the
#' coefficients of the polynomial whose roots are the elements of \code{x}.
#'
#' If an \eqn{N x N} square matrix is given, \code{poly((x)}
#' is the row vector of the coefficients of \code{det (z * diag (N) - x)},
#' which is the characteristic polynomial of \code{x}.
#'
#' @param x Real or complex vector, or square matrix.
#'
#' @return A vector of the coefficients of the polynomial in order from highest
#'   to lowest polynomial power.
#'
#' @examples
#' p <- poly(c(1, -1))
#' p <- poly(pracma::roots(1:3))
#' p <- poly(matrix(1:9, 3, 3))
#'
#' @seealso \code{\link{roots}}
#'
#' @author Kurt Hornik.\cr
#'  Conversion to R by Tom Short,\cr
#'   adapted by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}
#
#' @export

poly <- function(x) {

  n <- NROW(x)
  m <- NCOL(x)

  if (is.null(x) || length(x) == 0)
    return(1)

  if (m == 1) {
    v <- x
  } else if (m == n) {
    v <- eigen(x)$values
  } else {
    stop("x must be a vector or a square matrix")
  }

  y <- numeric(n + 1)
  y[1] <- 1
  for (j in seq_len(n)) {
    y[2:(j + 1)] <- y[2:(j + 1)] - v[j] * y[1:j]
  }

  zapIm(y)
}
