# interp.R
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
# 20201121  GvB       setup for gsignal v0.1.0
#------------------------------------------------------------------------------

#' Interpolation
#'
#' Increase sample rate by integer factor.
#'
#' @param x input data, specified as a numeric vector.
#' @param q interpolation factor, specified as a positive integer.
#' @param n Half the number of input samples used for interpolation, specified
#'   as a positive integer. For best results, use \code{n} no larger than 10.
#'   The low-pass interpolation filter has length \code{2 × n × q + 1}. Default:
#'   4.
#' @param Wc Normalized cutoff frequency of the input signal, specified as a
#'   positive real scalar not greater than 1 that represents a fraction of the
#'   Nyquist frequency. A value of 1 means that the signal occupies the full
#'   Nyquist interval. Default: 0.5.
#'
#' @return interpolated signal, returned as a vector.
#'
#' @examples
#' # Generate a signal
#' t <- seq(0, 2, 0.01)
#' x <- chirp(t, 2, .5, 10,'quadratic') + sin(2 * pi * t * 0.4)
#' w <- seq(1, 121, 4)
#' plot(t[w] * 1000, x[w], type = "h", xlab = "", ylab = "")
#' points(t[w] * 1000, x[w])
#' abline (h = 0)
#' y <- interp(x[seq(1, length(x), 4)], 4, 4, 1)
#' lines(t[1:121] * 1000, y[1:121], type = "l", col = "red")
#' points(t[1:121] * 1000, y[1:121], col = "red", pch = '+')
#' legend("topleft", legend = c("original", "interpolated"),
#'   lty = 1, pch = c(1, 3), col = c(1, 2))
#'
#' @seealso \code{\link{decimate}}, \code{\link{resample}}
#'
#' @author Paul Kienzle, \email{pkienzle@@users.sf.net}.\cr
#' Conversion to R by Geert van Boxtel \email{G.J.M.vanBoxtel@@gmail.com}.
#
#' @export

interp <- function(x, q, n = 4, Wc = 0.5) {

  if (!is.numeric(x) || !is.vector(x)) {
    stop("x must be a numeric vector")
  }
  if (!(isPosscal(q) && isWhole(q))) {
    stop("q must be a positive integer")
  }
  if (!(isPosscal(n) && isWhole(n))) {
    stop("n must be a positive integer")
  }
  if (!isPosscal(Wc) || Wc > 1) {
    stop("n must be a numeric value between 0 and 1")
  }

  y <- rep(0, length(x) * q + q * n + 1)
  y[seq(1, length(x) * q, q)] <- x
  b <- fir1(2 * q * n + 1, Wc / q)
  y <- q * fftfilt(b, y)
  y[- (1:(q * n + 1))]  # adjust for zero filter delay
}
