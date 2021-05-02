# kaiserord.R
# Copyright (C) 2019 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Octave signal package:
# Copyright (C) 2000 Paul Kienzle <pkienzle@users.sf.net>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; see the file COPYING. If not, see
# <https://www.gnu.org/licenses/>.
#
# 20200708 Geert van Boxtel          First version for v0.1.0
#------------------------------------------------------------------------------

#' Kaiser filter order and cutoff frequency
#'
#' Return the parameters needed to produce a FIR filter of the desired
#' specification from a Kaiser window.
#'
#' Given a set of specifications in the frequency domain, \code{kaiserord}
#' estimates the minimum FIR filter order that will approximately meet the
#' specifications. \code{kaiserord} converts the given filter specifications
#' into passband and stopband ripples and converts cutoff frequencies into the
#' form needed for windowed FIR filter design.
#'
#' \code{kaiserord} uses empirically derived formulas for estimating the orders
#' of lowpass filters, as well as differentiators and Hilbert transformers.
#' Estimates for multiband filters (such as band-pass filters) are derived from
#' the low-pass design formulas.
#'
#' The design formulas that underlie the Kaiser window and its application to
#' FIR filter design are
#' \deqn{\beta =}
#' \deqn{0.1102(\alpha - 8.7), \alpha > 50}
#' \deqn{0.5842(\alpha -21)^{0.4} + 0.07886(\alpha - 21), 21 \le \alpha \le 50}
#' \deqn{0, \alpha < 21}
#'
#' where \eqn{\alpha = -20log_{10}(\delta)} is the stopband attenuation
#' expressed in decibels, \eqn{n=(\alpha - 8) / 2.285(\Delta\omega)}, where
#' \eqn{n} is the filter order and \eqn{\Delta\omega} is the width of the
#' smallest transition region.
#'
#' @param f frequency bands, given as pairs, with the first half of the first
#'   pair assumed to start at 0 and the last half of the last pair assumed to
#'   end at 1. It is important to separate the band edges, since narrow
#'   transition regions require large order filters.
#' @param m magnitude within each band. Should be non-zero for pass band and
#'   zero for stop band. All passbands must have the same magnitude, or you will
#'   get the error that pass and stop bands must be strictly alternating.
#' @param dev deviation within each band. Since all bands in the resulting
#'   filter have the same deviation, only the minimum deviation is used. In this
#'   version, a single scalar will work just as well.
#' @param fs sampling rate. Used to convert the frequency specification into the
#'   c(0, 1) range, where 1 corresponds to the Nyquist frequency, \code{fs / 2}.
#'
#' @return A list of class \code{\link{FilterSpecs}} with the following list
#'   elements:
#' \describe{
#'   \item{n}{filter order}
#'   \item{Wc}{cutoff frequency}
#'   \item{type}{filter type, one of "low", "high", "stop", "pass", "DC-0", or
#'   "DC-1".}
#'   \item{beta}{shape parameter}
#' }
#'
#' @examples
#' fs <- 11025
#' op <- par(mfrow = c(2, 2), mar = c(3, 3, 1, 1))
#' for (i in 1:4) {
#'   if (i == 1) {
#'     bands <- c(1200, 1500)
#'     mag <- c(1, 0)
#'     dev <- c(0.1, 0.1)
#'   }
#'   if (i == 2) {
#'     bands <- c(1000, 1500)
#'     mag <- c(0, 1)
#'     dev <- c(0.1, 0.1)
#'   }
#'   if (i == 3) {
#'     bands <- c(1000, 1200, 3000, 3500)
#'     mag <- c(0, 1, 0)
#'     dev <- 0.1
#'   }
#'   if (i == 4) {
#'     bands <- 100 * c(10, 13, 15, 20, 30, 33, 35, 40)
#'     mag <- c(1, 0, 1, 0, 1)
#'     dev <- 0.05
#'   }
#'   kaisprm <- kaiserord(bands, mag, dev, fs)
#'   d <- max(1, trunc(kaisprm$n / 10))
#'   if (mag[length(mag)] == 1 && (d %% 2) == 1) {
#'      d <- d + 1
#'   }
#'   f1 <- freqz(fir1(kaisprm$n, kaisprm$Wc, kaisprm$type,
#'                    kaiser(kaisprm$n + 1, kaisprm$beta),
#'                    scale = FALSE),
#'               fs = fs)
#'   f2 <- freqz(fir1(kaisprm$n - d, kaisprm$Wc, kaisprm$type,
#'                    kaiser(kaisprm$n - d + 1, kaisprm$beta),
#'                    scale = FALSE),
#'               fs = fs)
#'   plot(f1$w, abs(f1$h), col = "blue", type = "l",  xlab = "", ylab = "")
#'   lines(f2$w, abs(f2$h), col = "red")
#'   legend("right", paste("order", c(kaisprm$n-d, kaisprm$n)),
#'          col = c("red", "blue"), lty = 1, bty = "n")
#'   b <- c(0, bands, fs/2)
#'   for (i in seq(2, length(b), by=2)) {
#'     hi <- mag[i/2] + dev[1]
#'     lo <- max(mag[i/2] - dev[1], 0)
#'     lines(c(b[i-1], b[i], b[i], b[i-1], b[i-1]), c(hi, hi, lo, lo, hi))
#'   }
#' }
#' par(op)
#'
#' @author Paul Kienzle.\cr
#'  Conversion to R by Tom Short,\cr
#'   adapted by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#'
#' @seealso \code{\link{hamming}}, \code{\link{kaiser}}
#'
#' @export

kaiserord <- function(f, m, dev, fs = 2) {

  ## parameter checking
  if (length(f) != 2 * length(m) - 2) {
    stop("One magnitude for each frequency band is required")
  }
  if (length(m) > 2 && any(m[1:(length(m) - 2)] != m[3:length(m)])) {
    stop("Ppass and stop bands must be strictly alternating")
  }
  if (length(dev) != length(m) && length(dev) != 1) {
    stop("One deviation for each frequency band is required")
  }
  dev <- min(dev)
  if (dev <= 0) {
    stop("dev must be  0")
  }
  if (!isPosscal(fs) || fs == 0) {
    stop("Sampling frequency fs must be a positive scalar")
  }

  ## use midpoints of the transition region for band edges
  w <- (f[seq(1, length(f), by = 2)] + f[seq(2, length(f), by = 2)]) / fs

  ## determine ftype
  if (length(w) == 1) {
    if (m[1] > m[2]) {
      ftype <- "low"
    } else {
      ftype <- "high"
    }
  } else if (length(w) == 2) {
    if (m[1] > m[2]) {
      ftype <- "stop"
    } else {
      ftype <- "pass"
    }
  } else {
    if (m[1] > m[2]) {
      ftype <- "DC-1"
    } else {
      ftype <- "DC-0"
    }
  }

  ## compute beta from dev
  A <- -20 * log10(dev)
  if (A > 50) {
    beta <- 0.1102 * (A - 8.7)
  } else if (A >= 21) {
    beta <- 0.5842 * (A - 21)^0.4 + 0.07886 * (A - 21)
  } else {
    beta <- 0.0
  }

  ## compute n from beta and dev
  dw <- 2 * pi * min(f[seq(2, length(f), by = 2)] -
                      f[seq(1, length(f), by = 2)]) / fs
  n <- max(1, ceiling((A - 8) / (2.285 * dw)))

  ## if last band is high, make sure the order of the filter is even.
  if ((m[1] > m[2]) == (length(w) %% 2 == 0) && n %% 2 == 1) {
    n <- n + 1
  }

  FilterSpecs(n = n, Wc = w, type = ftype, beta = beta)
}
