# peak2peak.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Octave signal package:
# Copyright (C) 2014 Georgios Ouzounis <ouzounis_georgios@hotmail.com>
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
# 20200110  GvB       setup for gsignal v0.1.0
#------------------------------------------------------------------------------

#' Maximum-to-minimum difference
#'
#' Compute the maximum-to-minimum difference of the input data \code{x}.
#'
#' The input \code{x} can be a vector, a matrix or an array. If the input is a
#' vector, a single value is returned representing the maximum-to-minimum
#' difference of the vector. If the input is a matrix or an array, a vector or
#' an array of values is returned representing the maximum-to-minimum
#' differences of the dimensions of \code{x} indicated by the \code{MARGIN}
#' argument.
#'
#' Support for complex valued input is provided. In this case, the function
#' \code{peak2peak} identifies the maximum and minimum in complex magnitude, and
#' then subtracts the complex number with the minimum modulus from the complex
#' number with the maximum modulus.
#'
#' @param x the data, expected to be a vector, a matrix, an array.
#' @param MARGIN a vector giving the subscripts which the function will be
#'   applied over. E.g., for a matrix 1 indicates rows, 2 indicates columns,
#'   c(1, 2) indicates rows and columns. Where \code{x} has named dimnames, it
#'   can be a character vector selecting dimension names. Default: 2 (columns)
#'
#' @return Vector or array of values containing the maximum-to-minimum
#'   differences of the specified \code{MARGIN} of \code{x}.
#'
#' @examples
#' ## numeric vector
#' x <- c(1:5)
#' pp <- peak2peak(x)
#' 
#' ## numeric matrix
#' x <- matrix(c(1,2,3, 100, 150, 200, 1000, 1500, 2000), 3, 3)
#' pp <- peak2peak(x)
#' pp <- peak2peak(x, 1)
#' 
#' ## numeric array
#' x <- array(c(1, 1.5, 2, 100, 150, 200, 1000, 1500, 2000,
#'              10000, 15000, 20000), c(2,3,2))
#' pp <- peak2peak(x, 1)
#' pp <- peak2peak(x, 2)
#' pp <- peak2peak(x, 3)
#' 
#' ## complex input
#' x <- c(1+1i, 2+3i, 3+5i, 4+7i, 5+9i)
#' pp <- peak2peak(x)
#'
#' @author Georgios Ouzounis, \email{ouzounis_georgios@@hotmail.com}.\cr
#' Conversion to R by Geert van Boxtel \email{G.J.M.vanBoxtel@@gmail.com}.
#
#' @export

peak2peak <- function(x, MARGIN = 2) {

  if (!(is.numeric(x) || is.complex(x)) ||
      !(is.vector(x) || is.matrix(x) || is.array(x))) {
    stop("x must be a numeric or complex vector, matrix or array")
  }
  if (!isPosscal(MARGIN) || !isWhole(MARGIN)) {
    stop("MARGIN must be a positive scalar")
  }

  mm <- function(a) max(a) - min(a)
  cmm <- function(a) a[which.max(abs(a))] - a[which.min(abs(a))]

  if (is.vector(x)) {
    x <- as.matrix(x)
    MARGIN <- 2
  }

  if (is.numeric(x)) {
    y <- apply(x, MARGIN, mm)
  } else {
    y <- apply(x, MARGIN, cmm)
  }
  y
}
