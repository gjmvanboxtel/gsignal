# decimate.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Original Octave code:
# Copyright (C) 2000 Paul Kienzle <pkienzle@users.sf.net>
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
# 20201129  GvB       setup for gsignal v0.1.0
# 20220328  GvB       copy dimnames of x to output object
# 20220513  GvB       Github Issue #7
# 20240821  GvB       Github Issue #17: replace fftfilt() with filtfilt()
#                     when ftype == 'fir'; changed order of arguments
#------------------------------------------------------------------------------

#' Decrease sample rate
#'
#' Downsample a signal by an integer factor.
#'
#' @param x input data, specified as a numeric vector or matrix. In case of a
#'   vector it represents a single signal; in case of a matrix each column is a
#'   signal.
#' @param q decimation factor, specified as a positive integer.
#' @param n Order of the filter used prior to the downsampling, specified as a
#'   positive integer. Default: 8 if \code{ftype} equals \code{"iir"}; 30 of
#'   \code{ftype} equals \code{"fir"}.
#' @param ftype filter type; either \code{"fir"}, specifying a FIR filter of
#'   length \code{n} designed with the function \code{\link{fir1}}, or
#'   \code{"iir"} (default), specifying an IIR Chebyshev filter of order 8 using
#'   the function \code{\link{cheby1}}.
#'
#' @return downsampled signal, returned as a vector or matrix.
#'
#' @examples
#' t <- seq(0, 2, 0.01)
#' x <- chirp(t, 2, .5, 10, 'quadratic') + sin(2 * pi * t * 0.4)
#' w <- 1:121
#' plot(t[w] * 1000, x[w], type = "h", col = "green")
#' points(t[w] * 1000, x[w], col = "green")
#' y = decimate(x, 4)
#' lines(t[seq(1, 121, 4)] * 1000, y[1:31], type = "h", col = "red")
#' points(t[seq(1, 121, 4)] * 1000, y[1:31], col = "red")
#'
#' @seealso \code{\link{cheby1}}, \code{\link{fir1}}
#'
#' @author Paul Kienzle, \email{pkienzle@@users.sf.net}.\cr
#'   Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#
#' @export

decimate <- function(x, q, ftype = c("iir", "fir"),
                     n = ifelse(ftype == "iir", 8, 30)) {

  if (!is.numeric(x)) {
    stop("x must be numeric")
  }

  if (is.vector(x)) {
    x <- matrix(x, ncol = 1)
    l <- length(x)
    vec <- TRUE
  } else if (is.matrix(x)) {
    vec <- FALSE
    l <- nrow(x)
  } else {
    stop("x must be a numeric vector or matrix")
  }
  ftype <- match.arg(ftype)
  if (!(isPosscal(q) && isWhole(q))) {
    stop("q must be a positive integer")
  }
  if (!(isPosscal(n) && isWhole(n))) {
    stop("n must be a positive integer")
  }

  if (ftype == "fir") {
    b <- fir1(n, 1 / q)
    y <- filtfilt(b, x)
  } else {
    ba <- cheby1(n, 0.05, 0.8 / q)
    y <- filtfilt(ba, x)
  }
  y <- y[seq(1, l, q), ]
  if (vec) {
    y <- as.vector(y)
  }
  dimnames(y) <- dimnames(x)
  y
}
