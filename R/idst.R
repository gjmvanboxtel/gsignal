# idst.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Original Octave code:
# Author: Paul Kienzle <pkienzle@users.sf.net> (2006)
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
# 20201016  GvB       setup for gsignal v0.1.0
#------------------------------------------------------------------------------

#' Inverse Discrete Sine Transform
#'
#' Compute the inverse discrete sine transform of a signal.
#'
#' The discrete sine transform (DST) is closely related to the discrete Fourier
#' transform. but using a purely real matrix. It is equivalent to the imaginary
#' parts of a DFT of roughly twice the length.
#'
#' @param x input discrete cosine transform, specified as a numeric vector or
#'   matrix. In case of a vector it represents a single signal; in case of a
#'   matrix each column is a signal.
#' @param n transform length, specified as a positive integer scalar. Default:
#'   \code{NROW(x)}.
#'
#' @return Inverse discrete sine transform, returned as a vector or matrix.
#'
#' @examples
#' x <- seq_len(100) + 50 * cos(seq_len(100) * 2 * pi / 40)
#' X <- dst(x)
#' xx <- idst(X)
#' all.equal(x, xx)
#'
#' @author Paul Kienzle, \email{pkienzle@@users.sf.net}.\cr
#' Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#'
#' @seealso \code{\link{dst}}
#'
#' @export

idst <- function(x, n = NROW(x)) {

  if (!(is.vector(x) || is.matrix(x)) || !(is.numeric(x) || is.complex(x))) {
    stop("x must be a numeric or complex vector or matrix")
  } else {
    realx <- is.numeric(x)
  }

  if (is.vector(x)) {
    vec <- TRUE
    x <- as.matrix(x, ncol = 1)
  } else {
    vec <- FALSE
  }
  nr <- nrow(x)

  if (!isPosscal(n) || !isWhole(n)) {
    stop("n must be a positive integer")
  }

  if (n != nr) {
    x <- postpad(x, n)
  }

  y <- dst(x, n) * 2 / (n + 1)

  if (realx) {
    y <- Re(y)
  }

  if (vec) {
    y <- as.vector(y)
  }
  y
}
