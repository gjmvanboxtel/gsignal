# dct.R
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
# 20210506  GvB       use matrix() instead of as.matrix()
#------------------------------------------------------------------------------

#' Discrete Sine Transform
#'
#' Compute the discrete sine transform of a signal.
#'
#' The discrete sine transform (DST) is closely related to the discrete Fourier
#' transform. but using a purely real matrix. It is equivalent to the imaginary
#' parts of a DFT of roughly twice the length.
#'
#' The DST has four standard variants. This function implements the DCT-I
#' according to the definition in [1], which is the most common variant, and
#' the original variant first proposed for image processing.
#'
#' The 'Matlab' documentation for the DST warns that the use of the function is
#' not recommended. They do not state the reason why, but it is likely that use
#' of the discrete cosine transform (DCT)is preferred for image processing.
#' Because cos(0) is 1, the first coefficient of the DCT (II) is the mean of the
#' values being transformed. This makes the first coefficient of each 8x8 block
#' represent the average tone of its constituent pixels, which is obviously a
#' good start. Subsequent coefficients add increasing levels of detail, starting
#' with sweeping gradients and continuing into increasingly fiddly patterns, and
#' it just so happens that the first few coefficients capture most of the signal
#' in photographic images. Sin(0) is 0, so the DSTs start with an offset of 0.5
#' or 1, and the first coefficient is a gentle mound rather than a flat plain.
#' That is unlikely to suit ordinary images, and the result is that DSTs require
#' more coefficients than DCTs to encode most blocks. This explanation was
#' provided by Douglas Bagnall on Stackoverflow.
#'
#' @param x input data, specified as a numeric vector or matrix. In case of a
#'   vector it represents a single signal; in case of a matrix each column is a
#'   signal.
#' @param n transform length, specified as a positive integer scalar. Default:
#'   \code{NROW(x)}.
#'
#' @return Discrete sine transform, returned as a vector or matrix.
#'
#' @examples
#' x <- matrix(seq_len(100) + 50 * cos(seq_len(100) * 2 * pi / 40))
#' ct <- dct(x)
#' st <- dst(x)
#'
#' @author Paul Kienzle, \email{pkienzle@@users.sf.net}.\cr
#' Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#'
#' @references [1] \url{https://en.wikipedia.org/wiki/Discrete_sine_transform}
#'
#' @seealso \code{\link{idst}}
#'
#' @export

 dst <- function(x, n = NROW(x)) {

  # check parameters
  if (!(is.vector(x) || is.matrix(x)) || !(is.numeric(x) || is.complex(x))) {
    stop("x must be a numeric or complex vector or matrix")
  } else {
    realx <- is.numeric(x)
  }

  if (is.vector(x)) {
    vec <- TRUE
    x <- matrix(x, ncol = 1)
  } else {
    vec <- FALSE
  }
  nr <- nrow(x)
  nc <- ncol(x)

  if (!isPosscal(n) || !isWhole(n)) {
    stop("n must be a positive integer")
  }

  if (n != nr) {
    x <- postpad(x, n)
  }

  y <- stats::mvfft(rbind(rep(0, nc), x, rep(0, nc),
                          -matrix(pracma::flipud(x), ncol = nc))) / -2i
  y <- y[2:(nr + 1), , drop = FALSE]

  if (realx) {
    y <- Re(y)
  }
  if (vec) {
    y <- as.vector(y)
  }
  y
}
