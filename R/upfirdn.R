# upfirdn.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
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
# 20200929  GvB       setup for gsignal v0.1.0
#------------------------------------------------------------------------------

#' Upsample, apply FIR filter, downsample
#'
#' Filter and resample a signal using polyphase interpolation.
#'
#' upfirdn performs a cascade of three operations:
#' \enumerate{
#'   \item Upsample the input data in the matrix \code{x} by a factor of the
#'     integer \code{p} (inserting zeros)
#'   \item FIR filter the upsampled signal data with the impulse response
#'     sequence given in the vector or matrix \code{h}
#'   \item Downsample the result by a factor of the integer \code{q} (throwing
#'     away samples)
#' }
#'
#' The FIR filter is usually a lowpass filter, which you must design using
#' another function such as \code{fir1}.
#'
#' @param x input data, specified as a numeric vector or matrix. In case of a
#'   vector it represents a single signal; in case of a matrix each column is a
#'   signal.
#' @param h Impulse response of the FIR filter specified as a numeric vector or
#'   matrix. If it is a vector, then it represents one FIR filter to may be
#'   applied to multiple signals in \code{x}; if it is a matrix, then each
#'   column is a separate FIR impulse response.
#' @param p Upsampling factor, specified as a positive integer (default: 1).
#' @param q downsampling factor, specified as a positive integer (default: 1).
#'
#' @return output signal, returned as a vector or matrix. Each column has length
#'   \code{ceiling(((length(x) - 1) * p + length(h)) / q)}.
#'
#' @note This function uses a polyphase implementation, which is generally
#'   faster than using \code{filter} by a factor equal to the downsampling
#'   factor, since it only calculates the needed outputs.
#'
#' @seealso \code{\link{fir1}}
#'
#' @examples
#' ## Change the sample rate of a signal by a rational conversion factor
#' ## from the DAT rate of 48 kHz to the CD sample rate of 44.1 kHz.
#' Fdat <- 48e3
#' Fcd <- 44.1e3
#' LM <- scan(text = capture.output(pracma::rats(Fcd / Fdat)), sep = '/',
#'  quiet = TRUE)
#' L <- LM[1]; M <- LM[2]
#'
#' ## Generate a 1.5 kHz sinusoid sampled at fDAT for 0.25 seconds. Plot the
#' ## first millisecond.
#' t <- seq(0, 0.25 - 1 / Fdat, 1 / Fdat)
#' x <- sin(2 * pi * 1.5e3 * t)
#' plot(t, x, type = "h", xlim = c(0, 0.001))
#' points(t, x)
#'
#' ## Design an antialiasing lowpass filter using a Kaiser window.
#' ## band edges 90% and 110%
#' f <- (Fdat / 2) * min(1 / L, 1 / M)
#' ko <- kaiserord(c(f - 0.1*f, f + 0.1*f), c(1, 0), c(0.1, 0.1), Fdat)
#' h <- L * fir1(ko$n, ko$Wc, ko$type, kaiser(ko$n + 1, ko$beta))
#'
#' y <- upfirdn(x, h, L, M)
#'
#' delay <- floor(((ko$n - 1) / 2 - (L - 1)) / L) + 1
#' y <- y[(delay + 1):length(y)]
#' ty <- seq(0, (length(y) - 1)) / Fcd
#' points(ty, y, col = "red", pch = 2)
#'
#' @author Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#
#' @export

upfirdn <- function(x, h, p = 1, q = 1) {

  # x and h must be numeric
  if (!is.numeric(x) || ! is.numeric(h)) {
    stop("x and h must be numeric")
  }

  if (is.vector(x)) {
    # if x is a vector then h must be a vector too
    if (!(is.vector(h) || "Ma" %in% class(h))) {
      stop("h must be a numeric vector")
    }
    ns <- 1
    x <- matrix(x, ncol = 1)
    h <- matrix(h, ncol = 1)
    vec <- TRUE
  } else if (is.matrix(x)) {
    ns <- ncol(x)
    if (is.vector(h)) {
      h <- matrix(rep(h, ns), ncol = ns, byrow = FALSE)
    } else if (!is.matrix(h)) {
      stop("h must be a numeric matrix")
    }
    vec <- FALSE
  } else {
    stop("x and h must be numeric vectors or matrices")
  }

  if (!(isPosscal(p) && isWhole(p)) ||
      !(isPosscal(q) && isWhole(q))) {
    stop("p and q must be positive integers")
  }

  y <- .Call("_gsignal_upfirdn", PACKAGE = "gsignal", x, h, p, q)

  if (vec) {
    y <- as.vector(y)
  }
  y
}
