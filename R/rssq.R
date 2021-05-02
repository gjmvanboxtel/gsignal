# rssq.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Octave signal package:
# Copyright (C) 2018-2019 Mike Miller
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
# 20200111  GvB       setup for gsignal v0.1.0
#------------------------------------------------------------------------------

#' Root-sum-of-squares
#'
#' Compute the root-sum-of-squares (SSQ) of the object \code{x}.
#'
#' The input \code{x} can be a vector, a matrix or an array. If the input is a
#' vector, a single value is returned representing the root-sum-of-squares of
#' the vector. If the input is a matrix or an array, a vector or an array of
#' values is returned representing the root-sum-of-squares of the dimensions of
#' \code{x} indicated by the \code{MARGIN} argument.
#'
#' Support for complex valued input is provided. The sum of squares of complex
#' numbers is defined by \code{sum(x * Conj(x))}
#'
#' @param x the data, expected to be a vector, a matrix, an array.
#' @param MARGIN a vector giving the subscripts which the function will be
#'   applied over. E.g., for a matrix 1 indicates rows, 2 indicates columns,
#'   c(1, 2) indicates rows and columns. Where \code{x} has named dimnames, it
#'   can be a character vector selecting dimension names. Default: 2 (usually
#'   columns)
#'
#' @return Vector or array of values containing the root-sum-of-squares of the
#'   specified \code{MARGIN} of \code{x}.
#'
#' @examples
#' ## numeric vector
#' x <- c(1:5)
#' p <- rssq(x)
#' 
#' ## numeric matrix
#' x <- matrix(c(1,2,3, 100, 150, 200, 1000, 1500, 2000), 3, 3)
#' p <- rssq(x)
#' p <- rssq(x, 1)
#' 
#' ## numeric array
#' x <- array(c(1, 1.5, 2, 100, 150, 200, 1000, 1500,
#'             2000, 10000, 15000, 20000), c(2,3,2))
#' p <- rssq(x, 1)
#' p <- rssq(x, 2)
#' p <- rssq(x, 3)
#' 
#' ## complex input
#' x <- c(1+1i, 2+3i, 3+5i, 4+7i, 5+9i)
#' p <- rssq(x)
#'
#' @author Mike Miller.\cr
#' Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#
#' @export

rssq <- function(x, MARGIN = 2) {

  if (!(is.numeric(x) || is.complex(x)) ||
      !(is.vector(x) || is.matrix(x) || is.array(x))) {
    stop("x must be a numeric or complex vector, matrix or array")
  }
  if (!isPosscal(MARGIN) || !isWhole(MARGIN)) {
    stop("MARGIN must be a positive scalar")
  }

  if (is.vector(x)) {
    x <- as.matrix(x)
    MARGIN <- 2
  }

  y <- apply(x, MARGIN, function(x) sqrt(ssq(x)))
  y
}
