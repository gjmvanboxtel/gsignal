# unwrap.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Original Octave function:
# Copyright (C) 2000-2017 Bill Lash
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
# 20200413  GvB       setup for gsignal v0.1.0
#---------------------------------------------------------------------------------------------------------------------

#' Unwrap phase angles
#' 
#' Unwrap radian phases by adding or subtracting multiples of \code{2 * pi}. 
#' 
#' @param x Input array, specified as a vector, matrix, or multidimensional array.
#' @param tol Jump threshold to apply phase shift, specified as a scalar. A jump
#'   threshold less than \eqn{pi} has the same effect as the threshold
#'   \eqn{pi}. Default: \deqn{pi}.
#' @param MARGIN A vector giving the subscripts which the function will be
#'   applied over. E.g., for a matrix 1 indicates rows, 2 indicates columns,
#'   c(1, 2) indicates rows and columns. Where x has named dimnames, it can be a
#'   character vector selecting dimension names. If MARGIN is larger than the
#'   dimensions of x, the result will have MARGIN dimensions. Default: 2
#'   (columns).
#' 
#' @return Unwraped phase angle, returned as a vector, matrix, or
#'   multidimensional array. The size of the output is always the same as the
#'   size of the input.
#' 
#' @examples
#' ## Define spiral shape.
#' t <- seq(0, 6 * pi, length.out = 201)
#' x <- t / pi * cos(t)
#' y <- t / pi * sin(t)
#' plot(x, y, type = "l")
#' ## find phase angle
#' p = atan2(y, x)
#' plot(t, p, type="l")
#' ## unwrap it
#' q = unwrap(p)
#' plot(t, q, type ="l")
#' 
#' @author Bill Lash, port to R by Geert van Boxtel
#'   \email{gjmvanboxtel@@gmail.com}
#' 
#' @export

unwrap <- function(x, tol = pi, MARGIN = 2) {
  
  if (is.vector(x)) {
    vec <- TRUE
  } else {
    vec <- FALSE
    dims <- dim(x)
    ld <- length(dims)
    if (is.character(MARGIN)) {
      if (is.null(dnn <- names(dimnames(x)))) 
        stop("'x' must have named dimnames")
      MARGIN <- match(MARGIN, dnn)
      if (anyNA(MARGIN)) 
        stop("not all elements of 'MARGIN' are names of dimensions")
    } else if(!isPosscal(MARGIN) || MARGIN > ld) {
      stop("'MARGIN' must be a positive integer and a valid margin")
    }
  }
  if (!is.numeric (x)) {
    stop("x must be a numeric matrix or vector");
  }
  tol = abs (tol)
  
  unw <- function(v, tol) {
    
    if (length(v) == 1) {
      y <- v
    } else {
      rng <- 2 * pi
      d <- diff(v)
      p <- round(abs(d) / rng) * rng * (((d > tol) > 0) - ((d < -tol) > 0))
      r <- cumsum(p)
      y <- v - c(0, r)
    }
    y
  }
  
  if (vec) {
    y <- unw(x, tol)
  } else {
    y <- unlist(apply(x, MARGIN, unw, tol = tol))
  }
  y
}
