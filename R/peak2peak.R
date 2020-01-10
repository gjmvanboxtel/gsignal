# peak2peak.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Octave signal package:
# Copyright (C) 2014 Georgios Ouzounis <ouzounis_georgios@hotmail.com>
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
# See also: http://www.gnu.org/licenses/gpl-2.0.txt
#
# Version history
# 20200110  GvB       setup for gsignal v0.1.0
#---------------------------------------------------------------------------------------------------------------------

#' Maximum-to-minimum difference
#' 
#' Compute the maximum-to-minimum difference of the input data \code{x}.
#' 
#' The input \code{x} can be a vector, a matrix or an array. If the input is a vector, a single value is returned representing
#' the maximum-to-minimum difference of the vector. If the input is a matrix or an array, a vector or an array of values is returned
#' representing the maximum-to-minimum differences of the dimensions of \code{x} indicated by the \code{MARGIN} argument.
#' 
#' @note This R implementation of the \code{peak2peak} function differs from the implementation of the \code{signal} package in Octave
#' in the following ways:
#' \enumerate{
#'   \item The parameter \code{dim} in the Octave function is replaced here by the \code{MARGIN} argument, which is interpreted
#'     as in the \code{apply} family of functions. That is, given a \code{r} by \code{c} matrix as input, 
#'     specifying \code{MARGIN = 1} will result in a vector of length \code{r} containing the maximum-to-minimum differences of the
#'     \code{r} rows, and specifying \code{MARGIN = 2} will result in a vector of length \code{c} containing the maximum-to-minimum
#'     differences of the \code{c} columns. In Octave, by contrast, \code{dim = 1} will give a row vector containing the differences
#'     over columns, and  \code{dim = 2} will give a column vector containing the differences over rows.
#'   \item Unlike the Octave function, the present implementation provides support for complex valued input. In this case, the 
#'     function \code{peak2peak} identifies the maximum and minimum in complex magnitude, and then subtracts the complex number
#'     with the minimum modulus from the complex number with the maximum modulus.
#' }
#' 
#' @param x the data, expected to be a vector, a matrix, an array.
#' @param MARGIN a vector giving the subscripts which the function will be applied over. E.g., for a matrix 1 indicates rows,
#'   2 indicates columns, c(1, 2) indicates rows and columns. Where \code{x} has named dimnames, it can be a character vector
#'   selecting dimension names. Default: 2 (usually columns)
#' 
#' @return vector or array of values containing the maximum-to-minimum differences of the specified \code{MARGIN} of \code{x}.
#'   See 'Details'.
#' 
#' @examples
#' ## numeric vector
#' x <- c(1:5)
#' peak2peak(x)
#' ## numeric matrix
#' x <- matrix(c(1,2,3, 100, 150, 200, 1000, 1500, 2000), 3, 3)
#' peak2peak(x)
#' peak2peak(x, 1)
#' ## numeric array
#' x <- array(c(1, 1.5, 2, 100, 150, 200, 1000, 1500, 2000, 10000, 15000, 20000), c(2,3,2))
#' peak2peak(x, 1)
#' peak2peak(x, 2)
#' peak2peak(x, 3)
#' ## complex input
#' x <- c(1+1i, 2+3i, 3+5i, 4+7i, 5+9i)
#' peak2peak(x)
#'
#' @author Original Octave code Copyright (C) 2014 Georgios Ouzounis \email{ouzounis_georgios@@hotmail.com}.
#' Port to R by Geert van Boxtel \email{G.J.M.vanBoxtel@@gmail.com}.
#
#' @export

peak2peak <- function (x, MARGIN = 2) {
  
  if (!(is.numeric(x) || is.complex(x)) || !(is.vector(x) || is.matrix(x) || is.array(x))) {
    stop ('x must be a vector, matrix or array containing numeric or complex values.')
  }
  if(!isPosscal(MARGIN) || !isWhole(MARGIN)) {
    stop ('MARGIN must be a positive scalar')
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

