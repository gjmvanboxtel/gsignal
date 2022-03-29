# upsample.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Original Octave code:
# Copyright (C) 2007 Paul Kienzle <pkienzle@users.sf.net>
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
# 20201121  GvB       setup for gsignal v0.1.0
# 20220328  GvB       copy dimnaes of x to output object
#------------------------------------------------------------------------------

#' Increase sample rate
#'
#' Upsample a signal by an integer factor.
#'
#' @param x input data, specified as a numeric vector or matrix. In case of a
#'   vector it represents a single signal; in case of a matrix each column is a
#'   signal.
#' @param n upsampling factor, specified as a positive integer. The signal is
#'   upsampled by inserting \code{n - 1} zeros between samples.
#' @param phase offset, specified as a positive integer from \code{0} to \code{n
#'   - 1}. Default: 0.
#'
#' @return Upsampled signal, returned as a vector or matrix.
#'
#' @examples
#' x <- seq_len(4)
#' u <- upsample(x, 3)
#' u <- upsample(x, 3, 2)
#'
#' x <- matrix(seq_len(6), 3, byrow = TRUE)
#' u <- upsample(x, 3)
#'
#' @seealso \code{\link{downsample}}, \code{\link{interp}},
#'   \code{\link{decimate}}, \code{\link{resample}}
#'
#' @author Paul Kienzle, \email{pkienzle@@users.sf.net}.\cr
#' Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#
#' @export

upsample <- function(x, n, phase = 0) {

  if (is.vector(x)) {
    ns <- 1
    lx <- length(x)
    x <- matrix(x, ncol = 1)
    vec <- TRUE
  } else if (is.matrix(x)) {
    ns <- ncol(x)
    lx <- nrow(x)
    vec <- FALSE
  } else {
    stop("x must be a numeric vector or matrix")
  }

  if (!(isPosscal(n) && isWhole(n))) {
    stop("n must be a positive integer")
  }
  if (!(isPosscal(phase) && isWhole(phase)) || phase > n - 1) {
    stop("phase must be a positive integer between 0 and n - 1")
  }

  y <- matrix(0, n * lx, ns)
  y[seq(phase + 1, length(y), n)] <- x

  if (vec) {
    y <- as.vector(y)
  }
  dimnames(y) <- dimnames(x)
  y
}
