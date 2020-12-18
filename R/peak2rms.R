# peak2rms.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Octave signal package:
# Copyright (C) 2015 Andreas Weber <octave@tech-chat.de>
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
#---------------------------------------------------------------------------------------------------------------------

#' Peak-magnitude-to-RMS ratio
#'
#' Compute the ratio of the largest absolute value to the root-mean-square (RMS)
#' value of the object \code{x}.
#'
#' The input \code{x} can be a vector, a matrix or an array. If the input is a
#' vector, a single value is returned representing the peak-magnitude-to-RMS
#' ratio of the vector. If the input is a matrix or an array, a vector or an
#' array of values is returned representing the peak-magnitude-to-RMS ratios of
#' the dimensions of \code{x} indicated by the \code{MARGIN} argument.
#'
#' Support for complex valued input is provided.
#'
#' @note This R implementation of the \code{peak2rms} function differs from the
#'   implementation of the \code{signal} package in Octave in that he parameter
#'   \code{dim} in the Octave function is replaced here by the \code{MARGIN}
#'   argument, which is interpreted as in the \code{apply} family of functions.
#'   That is, given a \code{r} by \code{c} matrix as input, specifying
#'   \code{MARGIN = 1} will result in a vector of length \code{r} containing the
#'   peak-magnitude-to-RMS ratios of the  \code{r} rows, and specifying
#'   \code{MARGIN = 2} will result in a vector of length \code{c} containing the
#'   peak-magnitude-to-RMS ratios of the \code{c} columns. In Octave, by
#'   contrast, \code{dim = 1} will give a row vector containing the ratios over
#'   columns, and  \code{dim = 2} will give a column vector containing the
#'   ratios over rows.
#'
#' @param x the data, expected to be a vector, a matrix, an array.
#' @param MARGIN a vector giving the subscripts which the function will be
#'   applied over. E.g., for a matrix 1 indicates rows, 2 indicates columns,
#'   c(1, 2) indicates rows and columns. Where \code{x} has named dimnames, it
#'   can be a character vector selecting dimension names. Default: 2 (usually
#'   columns)
#'
#' @return Vector or array of values containing the peak-magnitude-to-RMS ratios
#'   of the specified \code{MARGIN} of \code{x}.
#'
#' @examples
#' ## numeric vector
#' x <- c(1:5)
#' peak2rms(x)
#' ## numeric matrix
#' x <- matrix(c(1,2,3, 100, 150, 200, 1000, 1500, 2000), 3, 3)
#' peak2rms(x)
#' peak2rms(x, 1)
#' ## numeric array
#' x <- array(c(1, 1.5, 2, 100, 150, 200, 1000, 1500, 2000, 10000, 15000, 20000), c(2,3,2))
#' peak2rms(x, 1)
#' peak2rms(x, 2)
#' peak2rms(x, 3)
#' ## complex input
#' x <- c(1+1i, 2+3i, 3+5i, 4+7i, 5+9i)
#' peak2rms(x)
#'
#' @author Andreas Weber, \email{octave@@tech-chat.de}.\cr
#' Conversion to R by Geert van Boxtel \email{G.J.M.vanBoxtel@@gmail.com}.
#
#' @export

peak2rms <- function (x, MARGIN = 2) {

  if (!(is.numeric(x) || is.complex(x)) || !(is.vector(x) || is.matrix(x) || is.array(x))) {
    stop ('x must be a vector, matrix or array containing numeric or complex values.')
  }
  if(!isPosscal(MARGIN) || !isWhole(MARGIN)) {
    stop ('MARGIN must be a positive scalar')
  }

  if (is.vector(x)) {
    x <- as.matrix(x)
    MARGIN <- 2
  }

  p2r <- function (a) max(abs(a)) / rmsq(a)

  y <- apply(x, MARGIN, p2r)
  y
}

