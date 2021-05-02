# unwrap.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Original Octave function:
# Copyright (C) 2000-2017 Bill Lash
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
# 20200413  GvB       setup for gsignal v0.1.0
# 20210425  GvB       bugfix: handle NA, Nan, -Inf, Inf
#                     unwrap vector or matrix along columns
#------------------------------------------------------------------------------

#' Unwrap phase angles
#'
#' Unwrap radian phases by adding or subtracting multiples of \code{2 * pi}.
#'
#' @param x Input array, specified as a vector or a matrix. If \code{x} is a
#'   matrix, unwrapping along the columns of \code{x} is applied.
#' @param tol Jump threshold to apply phase shift, specified as a scalar. A jump
#'   threshold less than \eqn{pi} has the same effect as the threshold
#'   \eqn{pi}. Default: \deqn{pi}.
#'
#' @return Unwrapped phase angle, returned as a vector, matrix, or
#'   multidimensional array.
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
#' @author Bill Lash.\cr
#' Conversion to R by Geert van Boxtel, \email{gjmvanboxtel@@gmail.com}
#'
#' @export

unwrap <- function(x, tol = pi) {

  if (is.vector(x)) {
    x <- as.matrix(x, ncol = 1)
    vec <- TRUE
  } else {
    vec <- FALSE
  }
  nr <- nrow(x)
  nc <- ncol(x)
  if (!is.numeric (x)) {
    stop("x must be a numeric matrix or vector")
  }
  tol <- abs (tol)

  y <- x
  if (nr > 1) {
    rng <- 2 * pi
    for (col in seq_len(nc)) {
      valid <- which(is.finite(x[, col]))
      d <- diff(x[valid, col])
      p <- round(abs(d) / rng) * rng * (((d > tol) > 0) - ((d < -tol) > 0))
      r <- cumsum(p)
      y[valid, col] <- x[valid, col] - c(0, r)
    }
  }

  if (vec) {
    y <- as.vector(y)
  }
  y
}
