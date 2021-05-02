# ar_psd.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Original Octave function:
# Copyright (C) 2006 Peter V. Lanspeary <pvl@mecheng.adelaide.edu.au>
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
# 20201101  GvB       setup for gsignal v0.1.0
#------------------------------------------------------------------------------

#' Power spectrum of AR model
#'
#' Compute the power spectral density of an autoregressive model.
#'
#' This function calculates the power spectrum of the autoregressive model
#' \if{latex}{
#'   \deqn{x(n) = \sqrt{v} \cdot e(n) + \sum_{k=1}^{M} a(k) \cdot x(n-k)}
#' }
#' \if{html}{\preformatted{
#'                        M
#' x(n) = sqrt(v).e(n) + SUM a(k).x(n-k)
#'                       k=1
#' }}
#' where \code{x(n)} is the output of the model and \code{e(n)} is white noise.
#'
#' @param a numeric vector of autoregressive model coefficients. The first
#'   element is the zero-lag coefficient, which always has a value of 1.
#' @param v square of the moving average coefficient, specified as a positive
#'   scalar Default: 1
#' @param freq vector of frequencies at which power spectral density is
#'   calculated, or a scalar indicating the number of uniformly distributed
#'   frequency values at which spectral density is calculated. Default: 256.
#' @param fs sampling frequency (Hz). Default: 1
#' @param range character string. one of:
#' \describe{
#'   \item{\code{"half"} or \code{"onesided"}}{frequency range of the spectrum
#'   is from zero up to but not including \code{fs / 2}. Power from negative
#'   frequencies is added to the positive side of the spectrum.}
#'   \item{\code{"whole"} or \code{"twosided"}}{frequency range of the spectrum
#'   is \code{-fs / 2} to \code{fs / 2}, with negative frequencies stored in
#'   "wrap around order" after the positive frequencies; e.g. frequencies for a
#'   10-point \code{"twosided"} spectrum are 0 0.1 0.2 0.3 0.4 0.5 -0.4 -0.3
#'   -0.2. -0.1.}
#'   \item{\code{"shift"} or \code{"centerdc"}}{same as \code{"whole"} but with
#'   the first half of the spectrum swapped with second half to put the
#'   zero-frequency value in the middle. If \code{freq} is a vector,
#'   \code{"shift"} is ignored.}
#' }
#'   Default: If model coefficients \code{a} are real, the default range is
#'   \code{"half"}, otherwise the default range is \code{"whole"}.
#' @param method method used to calculate the power spectral density, either
#'   \code{"fft"} (use the Fast Fourier Transform) or \code{"poly"} (calculate
#'   the power spectrum as a polynomial). This argument is ignored if the
#'   \code{freq} argument is a vector. The default is \code{"poly"} unless the
#'   \code{freq} argument is an integer power of 2.
#' @param x object to plot.
#' @param yscale character string specifying scaling of Y-axis; one of
#'   \code{"linear"}, \code{"log"}, \code{"dB"}
#' @param xlab,ylab,main labels passed to plotting function. Default: NULL
#' @param ... additional arguments passed to functions
#'
#' @return An object of class \code{"ar_psd"} , which is a list containing two
#'   elements, \code{freq} and \code{psd} containing the frequency values and
#'   the estimates of power-spectral density, respectively.
#'
#' @examples
#' a <- c(1, -2.7607, 3.8106, -2.6535, 0.9238)
#' psd <- ar_psd(a)
#'
#' @author Peter V. Lanspeary, \email{pvl@@mecheng.adelaide.edu.au}.\cr
#'  Conversion to R by Geert van Boxtel, \email{gjmvanboxtel@@gmail.com}
#'
#' @rdname ar_psd
#' @export

ar_psd <- function(a, v = 1, freq = 256, fs = 1,
                   range = ifelse(is.numeric(a), "half", "whole"),
                   method = ifelse(length(freq) == 1 &&
                                     bitwAnd(freq, freq - 1) == 0,
                                   "fft", "poly"))  {

  # parameter checking
  if (!is.vector(a) || length(a) < 2 || a[1] != 1) {
    stop("a must be a vector of length >= 2 with the first element equal to 1")
  } else {
    real_model <- ifelse(is.numeric(a), 1L, 0L)
  }

  if (!isPosscal(v) || v <= 0) {
    stop("v must be a positive scalar > 0")
  }
  if (!is.vector(freq)) {
    stop("freq must be a scalar or a vector")
  } else {
    freq_len <- length(freq)
    user_freqs <- freq_len > 1
    if (!user_freqs && (!is.numeric(freq) || !isWhole(freq) || freq < 2)) {
      stop("freq must be an integer >= 2")
    } else if (user_freqs && !is.numeric(freq)) {
      stop("freq vector must be numeric")
    }
  }

  if (!isPosscal(fs) || fs <= 0) {
    stop("fs must be a positive scalar > 0")
  }

  range <- match.arg(range, c("half", "onesided", "whole",
                              "twosided", "shift", "centerdc"))
  if (range == "half" || range == "onesided") {
    pad_fact <- 2L    # FT zero-padding factor (pad FFT to double length)
    do_shift <- FALSE
  } else if (range == "whole" || range == "twosided") {
    pad_fact <- 1L    # FFT zero-padding factor (do not pad)
    do_shift <- FALSE
  } else if (range == "shift" || range == "centerdc") {
    pad_fact <- 1L
    do_shift <- TRUE
  }

  method <- match.arg(method, c("fft", "poly"))
  if (method == "fft") {
    force_fft  <- TRUE
    force_poly <- FALSE
  } else if (method == "poly") {
    force_fft  <- FALSE
    force_poly <- TRUE
  }
  # end of parameter checking

  # frequencies at which to determine psd
  if (user_freqs) {
    # user provides vector of frequencies
    if (any(abs(freq) > fs / 2)) {
      stop("freq cannot exceed half of the sampling frequency")
    } else if (pad_fact == 2L && any(freq < 0)) {
        stop("freq must be positive in a onesided spectrum")
    }
    freq_len <- length(freq)
    fft_len  <- freq_len
    use_fft  <- FALSE
    do_shift <- FALSE
  } else {
    # internally generated frequencies
    freq_len <- freq
    freq <- (fs / pad_fact / freq_len) * seq(0, freq_len - 1)
    # decide which method to use (poly or FFT)
    is_power_of_2 <- length(freq_len) == 1 &&
      bitwAnd(freq_len, freq_len - 1) == 0
    use_fft <- (!force_poly && is_power_of_2) || force_fft
    fft_len <- freq_len * pad_fact
  }

  # calculate denominator of Equation 2.28, Kay and Marple Jr, "Spectrum
  # analysis -- a modern perspective", Proceedings of the IEEE, Vol 69, pp
  # 1380-1419, Nov., 1981
  len_coeffs <- length(a)
  if (use_fft) {
    # FFT method
    fft_out <- stats::fft(postpad(a, fft_len))
  } else {
    # polynomial method
    # complex data on "half" frequency range needs -ve frequency values
    if (pad_fact == 2L && !real_model) {
      freq <- c(freq, -freq[seq(freq_len, 1, -1)])
      fft_len <- 2 * freq_len
    }
    fft_out <- pracma::polyval(a[seq(len_coeffs, 1, -1)],
                               exp((-1i * 2 * pi / fs) * freq))
  }

  # The power spectrum (PSD) is the scaled squared reciprocal of amplitude
  # of the FFT/polynomial. This is NOT the reciprocal of the periodogram.
  # The PSD is a continuous function of frequency.  For uniformly
  # distributed frequency values, the FFT algorithm might be the most
  # efficient way of calculating it.
  psd <- (v / fs) / (fft_out * Conj(fft_out))

  # range='half' or 'onesided',
  #   add PSD at -ve frequencies to PSD at +ve frequencies
  # N.B. unlike periodogram, PSD at zero frequency _is_ doubled.
  if (pad_fact == 2L) {
    freq <- freq[1:freq_len]
    if (real_model) {
      # real data, double the psd
      psd <- 2 * psd[1:freq_len]
    } else if (use_fft) {
      # complex data, FFT method, internally-generated frequencies
      psd <- psd[1:freq_len] + c(psd[1], psd[seq(fft_len, freq_len + 2, -1)])
    } else {
      # complex data, polynomial method
      # user-defined and internally-generated frequencies
      psd <- psd[1:freq_len] + psd[seq(fft_len, freq_len + 1, -1)]
    }

  # range equals 'shift'
  #   disabled for user-supplied frequencies
  #   Shift zero-frequency to the middle (pad_fact==1)
  } else if (do_shift) {
    len2 <- trunc((fft_len + 1) / 2)
    psd  <- c(psd[(len2 + 1):fft_len], psd[1:len2])
    freq <- c(freq[(len2 + 1):fft_len] - fs, freq[1:len2])
  }

  if (real_model == 1L) {
    psd <- Re(psd)
  }

  structure(list(freq = freq, psd = psd, fs = fs), class = "ar_psd")
}

#' @rdname ar_psd
#' @export

plot.ar_psd <- function(x, yscale = c("linear", "log", "dB"),
                        xlab = NULL, ylab = NULL, main = NULL, ...) {

  if (!("ar_psd" %in% class(x))) {
    stop("invalid object type")
  }
  yscale <- match.arg(yscale)

  if (is.null(xlab)) {
    if (x$fs == 1) {
      xlab <- expression(
        paste("Normalized frequency (\u00D7 ",
              pi, " rad/sample)"))
    } else if (x$fs == pi) {
      xlab <- "Frequency (rad/sample)"
    } else {
      xlab <- "Frequency (Hz)"
    }
  }
  sub <- paste("Resolution:", format(x$fs / length(x$freq),
                                     digits = 6, nsmall = 6))

  if (is.null(ylab)) {
    ylab <- switch(yscale,
                   "linear" = "PSD/Frequency",
                   "log" = expression(paste("log"[10], "PSD/Frequency")),
                   "dB" = "PSD (dB/Hz)")
  }
  plt <- switch(yscale,
                "linear" = x$psd,
                "log" = log10(x$psd),
                "dB" = 10 * log10(x$psd))
  graphics::plot(x$freq, plt, type = "l", xlab = xlab, ylab = "", ...)
  graphics::title(main, sub = sub)
  graphics::title(ylab = ylab, line = 2)
}

print.ar_psd <- function(x, yscale = c("linear", "log", "dB"),
                         xlab = NULL, ylab = NULL, main = NULL, ...) {
  plot.ar_psd(x, yscale, xlab, ylab, main, ...)
}