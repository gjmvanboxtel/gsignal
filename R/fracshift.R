# fracshift.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Original Octave code:
# Copyright (C) 2008 Eric Chassande-Mottin, CNRS (France)
# Copyright (C) 2018 Juan Pablo Carbajal
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
# 20201122  GvB       setup for gsignal v0.1.0
#------------------------------------------------------------------------------

#' Fractional shift
#'
#' Shift a signal by a (possibly fractional) number of samples.
#'
#' The function calculates the initial index and end index of the sequences of
#' 1’s in the rows of \code{x}. The clusters are sought in the rows of the array
#' \code{x}. The function works by finding the indexes of jumps between
#' consecutive values in the rows of \code{x}.
#'
#' @param x input data, specified as a numeric vector.
#' @param d number of samples to shift \code{x} by, specified as a numeric value
#' @param h interpolator impulse response, specified as a numeric vector. If
#'   NULL (default), the interpolator is designed by a Kaiser-windowed sinecard.
#'
#' @return A list of matrices size \code{nr}, where \code{nr} is the number of
#'   rows in \code{x}. Each element of the list contains a matrix with two rows.
#'   The first row is the initial index of a sequence of 1’s and the second row
#'   is the end index of that sequence. If \code{nr == 1} the output is a matrix
#'   with two rows.
#'
#' @examples
#' N = 1024
#' t <- seq(0, 1, length.out = N)
#' x <- exp(-t^2 / 2 / 0.25^2) * sin(2 * pi * 10 * t)
#' dt <- 0.25
#' d  <- dt / (t[2] - t[1])
#' y <- fracshift(x, d)
#' plot(t, x, type = "l", xlab = "Time", ylab = "Sigfnal")
#' lines (t, y, col = "red")
#' legend("topright", legend = c("original", "shifted"), lty = 1, col = 1:2)
#'
#' @author Eric Chassande-Mottin, \email{ecm@@apc.univ-paris7.fr},\cr
#'  Juan Pablo Carbajal, \email{carbajal@@ifi.uzh.ch},\cr
#'  Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#'
#' @references [1] A. V. Oppenheim, R. W. Schafer and J. R. Buck,
#' Discrete-time signal processing, Signal processing series,
#' Prentice-Hall, 1999.\cr
#' [2] T.I. Laakso, V. Valimaki, M. Karjalainen and U.K. Laine
#' Splitting the unit delay, IEEE Signal Processing Magazine,
#' vol. 13, no. 1, pp 30--59 Jan 1996.
#
#' @export

fracshift <- function(x, d, h = NULL) {

  if (!is.vector(x) || !is.numeric(x)) {
    stop("x must be a numeric vector")
  }
  if (!isScalar(d)) {
    stop("d must be a scalar")
  }
  if (!is.null(h)) {
    if (!is.numeric(h) || !is.vector(h)) {
      stop("h must be a numeric vector")
    }
  } else {
    h <- design_filter(d)

    Lx <- length(x)
    Lh <- length(h)
    L  <- (Lh - 1) / 2.0
    Ly <- Lx

    ## pre and postpad filter response
    hpad   <- prepad(h, Lh)
    offset <- floor(L)
    hpad   <- postpad(hpad, Ly + offset)

    ## filtering
    xfilt <- upfirdn(x, hpad, 1, 1)
    x     <- xfilt[(offset + 1):(offset + Ly)]
  }

  y <- pracma::circshift(x, trunc(d))


  y
}

design_filter <- function(d) {

  ## properties of the interpolation filter
  log10_rejection <- -3.0
  ## use empirical formula from [1] Chap 7, Eq. (7.63) p 476
  rejection_dB <- -20.0 * log10_rejection
  ## determine parameter of Kaiser window
  ## use empirical formula from [1] Chap 7, Eq. (7.62) p 474
  ## FIXME since the parameters are fix the conditional below is not needed
  if (rejection_dB >= 21 && rejection_dB <= 50) {
    beta <- 0.5842 * (rejection_dB - 21.0)^0.4 + 0.07886 * (rejection_dB - 21.0)
  } else if (rejection_dB > 50) {
    beta <- 0.1102 * (rejection_dB - 8.7)
  } else {
    beta <- 0.0
  }
  ## properties of the interpolation filter
  stopband_cutoff_f <- 0.5
  roll_off_width    <- stopband_cutoff_f / 10

  ## ideal sinc filter
  ## determine filter length
  L <- ceiling((rejection_dB - 8.0) / (28.714 * roll_off_width))
  t <- (-L:L)
  ideal_filter <- 2 * stopband_cutoff_f *
    sinc(2 * stopband_cutoff_f * (t - (d - trunc(d))))

  ## apodize ideal (sincard) filter response
  m <- 2 * L
  t <- (0:m) - (d - trunc(d))
  # kludge to prevent sqrt of negative number
  # besselI does not take complex input
  qq <- t * (m - t)
  qq[which(qq < 0)] <- 0
  t <- 2 * beta / m * sqrt(qq)
  w <- besselI(t, 0) / besselI(beta, 0)
  h <- w * ideal_filter
  h
}
