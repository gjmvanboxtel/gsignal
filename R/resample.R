# resample.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Original Octave code:
# Copyright (C) 2008 Eric Chassande-Mottin, CNRS (France)
# <ecm@apc.univ-paris7.fr>
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

#' Change sampling rate
#'
#' Resample using a polyphase algorithm.
#'
#' If \code{h} is not specified, this function will design an optimal FIR filter
#' using a Kaiser-Bessel window. The filter length and the parameter \eqn{\beta}
#' are computed based on ref [2], Chapter 7, Eq. 7.63 (p. 476), and Eq. 7.62 (p.
#' 474), respectively.
#'
#' @param x input data, specified as a numeric vector or matrix. In case of a
#'   vector it represents a single signal; in case of a matrix each column is a
#'   signal.
#' @param p,q resampling factors, specified as positive integers. \code{p / q}
#'   is the resampling factor.
#' @param h Impulse response of the FIR filter specified as a numeric vector or
#'   matrix. If it is a vector, then it represents one FIR filter to may be
#'   applied to multiple signals in \code{x}; if it is a matrix, then each
#'   column is a separate FIR impulse response. If not specified, a FIR filter
#'   based on a Kaiser window is designed.
#'
#' @return output signal, returned as a vector or matrix. Each column has length
#'   \code{ceiling(((length(x) - 1) * p + length(h)) / q)}..
#'
#' @examples
#' lx <- 60
#' tx <- seq(0, 360, length.out = lx)
#' x <- sin(2 * pi * tx / 120)
#'
#' # upsample
#' p <- 3; q <- 2
#' ty <- seq(0, 360, length.out = lx * p / q)
#' y <- resample(x, p, q)
#'
#' # downsample
#' p <- 2; q <- 3
#' tz <- seq(0, 360, length.out = lx * p / q)
#' z <- resample(x, p, q)
#'
#' # plot
#' plot(tx, x, type = "b", col = 1, pch = 1,
#'  xlab = "", ylab = "")
#' points(ty, y, col = 2, pch = 2)
#' points(tz, z, col = 3, pch = 3)
#' legend("bottomleft", legend = c("original", "upsampled", "downsampled"),
#'   lty = 1, pch = 1:3, col = 1:3)
#'
#' @references [1] Proakis, J.G., and Manolakis, D.G. (2007).
#' Digital Signal Processing: Principles, Algorithms, and Applications,
#' 4th ed., Prentice Hall, Chap. 6.\cr
#' [2] Oppenheim, A.V., Schafer, R.W., and Buck, J.R. (1999).
#' Discrete-time signal processing, Signal processing series,
#' Prentice-Hall.
#'
#' @seealso \code{\link{kaiser}}
#'
#' @author Eric Chassande-Mottin, \email{ecm@@apc.univ-paris7.fr}.\cr
#' Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#
#' @export

resample <- function(x, p, q, h) {

  if (!is.numeric(x)) {
    stop("x must be numeric")
  }

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

  if (!(isPosscal(p) && isWhole(p)) ||
      !(isPosscal(q) && isWhole(q))) {
    stop("p and q must be positive integers")
  }

  # simplify decimation and interpolation factors
  great_common_divisor <- pracma::gcd(p, q)
  if (great_common_divisor > 1) {
    p <- as.double(p) / as.double(great_common_divisor)
    q <- as.double(q) / as.double(great_common_divisor)
  } else {
    p <- as.double(p)
    q <- as.double(q)
  }

  # filter design if required
  if (missing(h)) {

    # properties of the antialiasing filter
    log10_rejection <- -3.0
    stopband_cutoff_f <- 1 / (2 * max(p, q))
    roll_off_width <- stopband_cutoff_f / 10.0

    # determine filter length
    # use empirical formula from ref [2], Chap 7, Eq. (7.63) p 476
    rejection_dB <- -20.0 * log10_rejection
    L <- ceiling((rejection_dB - 8.0) / (28.714 * roll_off_width))

    ## ideal sinc filter
    t <- -L:L
    ideal_filter <- 2 * p * stopband_cutoff_f * sinc(2 * stopband_cutoff_f * t)

    # determine parameter of Kaiser window
    # use empirical formula from [2] Chap 7, Eq. (7.62) p 474
    if ((rejection_dB >= 21) && (rejection_dB <= 50)) {
      beta <- 0.5842 * (rejection_dB - 21.0)^0.4 +
        0.07886 * (rejection_dB - 21.0)
    } else if (rejection_dB > 50) {
      beta <- 0.1102 * (rejection_dB - 8.7)
    } else {
      beta <- 0.0
    }

    # apodize ideal filter response
    h <- kaiser(2 * L + 1, beta) * ideal_filter
  }

  if (is.vector(h)) {
    lh <- length(h)
    h <- matrix(rep(h, ns), ncol = ns, byrow = FALSE)
  } else if (!is.matrix(h)) {
    stop("h must be a numeric matrix")
  }
  lh <- nrow(h)

  L <- (lh - 1) / 2.0
  ly <- ceiling(lx * p / q)

  # pre and postpad filter response
  nz_pre <- floor(q - L %% q)
  hpad <- prepad(h, lh + nz_pre)
  offset <- floor((L + nz_pre) / q)
  nz_post <- 0
  while (ceiling(((lx - 1) * p + nz_pre + lh + nz_post) / q) - offset < ly) {
    nz_post <- nz_post + 1
  }
  hpad <- postpad(hpad, lh + nz_pre + nz_post)

  ##filtering
  xfilt <- upfirdn(x, hpad, p, q)
  y <- xfilt[(offset + 1):(offset + ly), ]
  if (vec) {
    y <- as.vector(y)
  }
  y
}
