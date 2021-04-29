# idct2.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Original Octave code:
# Copyright (C) 2001 Paul Kienzle <pkienzle@users.sf.net>
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
# 20201015  GvB       setup for gsignal v0.1.0
#------------------------------------------------------------------------------

#' Inverse 2-D Discrete Cosine Transform
#'
#' Compute the inverse two-dimensional discrete cosine transform of a matrix.
#'
#' The discrete cosine transform (DCT) is closely related to the discrete
#' Fourier transform. It is a separable linear transformation; that is, the
#' two-dimensional transform is equivalent to a one-dimensional DCT performed
#' along a single dimension followed by a one-dimensional DCT in the other
#' dimension.
#'
#' @param x 2-D numeric matrix
#' @param m Number of rows, specified as a positive integer. \code{dct2} pads or
#'   truncates \code{x} so that it has \code{m} rows. Default: \code{NROW(x)}.
#' @param n Number of columns, specified as a positive integer. \code{dct2} pads
#'   or truncates \code{x} so that it has \code{n} columns. Default:
#'   \code{NCOL(x)}.
#'
#' @return \code{m}-by-\code{n} numeric discrete cosine transformed matrix.
#'
#' @examples
#' A <- matrix(50 * runif(100), 10, 10)
#' B <- dct2(A)
#' B[which(B < 1)] <- 0
#' AA <- idct2(B)
#'
#' @author Paul Kienzle, \email{pkienzle@@users.sf.net}.\cr
#' Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#'
#' @seealso \code{\link{dct2}}
#'
#' @export

 idct2 <- function(x, m = NROW(x), n = NCOL(x)) {

  # check parameters
  if (!is.matrix(x) || !is.numeric(x)) {
    stop("x must be a numeric matrix")
  }
  nr <- nrow(x)
  nc <- ncol(x)

  if (!isPosscal(m) || !isWhole(m)) {
    stop("m must be a positive integer")
  }
  if (!isPosscal(n) || !isWhole(n)) {
    stop("n must be a positive integer")
  }

  if (m != nr) {
    x <- postpad(x, n, MARGIN = 2)
  }
  if (n != nc) {
    x <- postpad(x, n, MARGIN = 1)
  }

  if (m == 1) {
    y <- t(idct(t(x), n))
  } else if (n == 1) {
    y <- idct(x, m)
  } else {
    y <- t(idct(t(idct(x, m)), n))
  }
  y
}
