# wconv.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Original Octave function Copyright (C) 2013 Lukas F. Reichlin
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
# 20200228  GvB       setup for gsignal v0.1.0
#------------------------------------------------------------------------------

#' 1-D or 2-D convolution
#'
#' Compute the one- or two-dimensional convolution of two vectors or matrices.
#'
#' @param type Numeric or character, specifies the type of convolution to
#'   perform:
#' \describe{
#'   \item{"1d"}{For \code{a} and \code{b} as (coerced to) vectors,
#'   perform 1-D convolution of \code{a} and \code{b};}
#'   \item{"2d}{For \code{a} and \code{b} as (coerced to)
#'   matrices, perform 2-D convolution of \code{a} and \code{b};}
#'   \item{"row"}{For \code{a} as (coerced to) a matrix, and \code{b}
#'   (coerced to) a vector, perform the 1-D convolution of the rows of \code{a}
#'   and \code{b};}
#'   \item{"column"}{For \code{a} as (coerced to) a matrix, and \code{b}
#'   (coerced to) a vector, perform the 1-D convolution of the colums of
#'   \code{a} and \code{b};}
#' }
#' @param a,b Input vectors or matrices, coerced to numeric.
#' @param shape Subsection of convolution, partially matched to:
#' \describe{
#'   \item{"full"}{Return the full convolution (default)}
#'   \item{"same"}{Return the central part of the convolution with the same size
#'   as A. The central part of the convolution begins at the indices
#'   \code{floor(c(nrow(b), ncol(b)) / 2 + 1)}}
#'   \item{"valid"}{Return only the parts which do not include zero-padded edges.
#'   The size of the result is \code{max(c(nrow(a), ncol(b)) - c(nrow(b),
#'   ncol(b)) + 1, 0)}}
#' }
#'
#' @return Convolution of input matrices, returned as a matrix or a vector.
#'
#' @examples
#' a <- matrix(1:16, 4, 4)
#' b <- matrix(1:9, 3,3)
#' w <- wconv('2', a, b)
#' w <- wconv('1', a, b, 'same')
#' w <- wconv('r', a, b)
#' w <- wconv('r', a, c(0,1), 'same')
#' w <- wconv('c', a, c(0,1), 'valid')
#'
#' @seealso \code{\link{conv}}
#'
#' @author Lukas Reichlin, \email{lukas.reichlin@@gmail.com}.\cr
#' Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#'
#' @export

wconv <- function(type = c("1d", "2d", "row", "column"),
                  a, b, shape = c("full", "same", "valid")) {

  type <- match.arg(type)
  shape <- match.arg(shape)

  y <- switch(type,
              "1d" = conv(as.vector(a), as.vector(b), shape),
              "2d" = conv2(as.matrix(a), as.matrix(b), shape),
              "row" = conv2(as.matrix(a), t(as.vector(b)), shape),
              "column" = t(conv2(t(as.matrix(a)), t(as.vector(b)), shape))
  )
  y
}
