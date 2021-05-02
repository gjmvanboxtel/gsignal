# specgram.R
# Copyright (C) 2019 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Octave signal package:
# Copyright (C) 1999-2001 Paul Kienzle <pkienzle@users.sf.net>
# Original conversion to R by Tom Short
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
# 20191204  Geert van Boxtel          First version for v0.1.0
# 20210418  GvB                       v0.3.0 plotting via S3 method
#------------------------------------------------------------------------------

#' Spectrogram
#'
#' Spectrogram using short-time Fourier transform.
#'
#' Generate a spectrogram for the signal \code{x}. The signal is chopped into
#' overlapping segments of length \code{n}, and each segment is windowed and
#' transformed into the frequency domain using the FFT. The default segment size
#' is 256. If \code{fs} is given, it specifies the sampling rate of the input
#' signal. The argument \code{window} specifies an alternate window to apply
#' rather than the default of \code{hanning(n)}. The argument overlap specifies
#' the number of samples overlap between successive segments of the input
#' signal. The default overlap is \code{length (window)/2}.
#'
#' When results of \code{specgram} are printed, a spectrogram will be plotted.
#' As with \code{lattice} plots, automatic printing does not work inside loops
#' and function calls, so explicit calls to \code{print} or \code{plot} are
#' needed there.
#'
#' The choice of window defines the time-frequency resolution. In speech for
#' example, a wide window shows more harmonic detail while a narrow window
#' averages over the harmonic detail and shows more formant structure. The shape
#' of the window is not so critical so long as it goes gradually to zero on the
#' ends.
#'
#' Step size (which is window length minus overlap) controls the horizontal
#' scale of the spectrogram. Decrease it to stretch, or increase it to compress.
#' Increasing step size will reduce time resolution, but decreasing it will not
#' improve it much beyond the limits imposed by the window size (you do gain a
#' little bit, depending on the shape of your window, as the peak of the window
#' slides over peaks in the signal energy). The range 1-5 msec is good for
#' speech.
#'
#' FFT length controls the vertical scale. Selecting an FFT length greater than
#' the window length does not add any information to the spectrum, but it is a
#' good way to interpolate between frequency points which can make for prettier
#' spectrograms.
#'
#' AFTER you have generated the spectral slices, there are a number of decisions
#' for displaying them. First the phase information is discarded and the energy
#' normalized:
#'
#' \code{S <- abs(S); S <- S / max(S)}
#'
#' Then the dynamic range of the signal is chosen. Since information in speech
#' is well above the noise floor, it makes sense to eliminate any dynamic range
#' at the bottom end. This is done by taking the max of the magnitude and some
#' minimum energy such as minE = -40dB. Similarly, there is not much information
#' in the very top of the range, so clipping to a maximum energy such as maxE =
#' -3dB makes sense:
#'
#' \code{S <- max(S, 10^(minE / 10)); S <- min(S, 10^(maxE / 10))}
#'
#' The frequency range of the FFT is from 0 to the Nyquist frequency of one half
#' the sampling rate. If the signal of interest is band limited, you do not need
#' to display the entire frequency range. In speech for example, most of the
#' signal is below 4 kHz, so there is no reason to display up to the Nyquist
#' frequency of 10 kHz for a 20 kHz sampling rate. In this case you will want to
#' keep only the first 40% of the rows of the returned \code{S} and \code{f}.
#' More generally, to display the frequency range from minF to maxF, you could
#' use the following row index:
#'
#' \code{idx <- (f >= minF & f <= maxF)}
#'
#' Then there is the choice of colormap. A brightness varying colormap such as
#' copper or bone gives good shape to the ridges and valleys. A hue varying
#' colormap such as jet or hsv gives an indication of the steepness of the
#' slopes. In the field that I am working in (neuroscience / electrophysiology)
#' rainbow color palettes such as jet are very often used. This is an
#' unfortunate choice mainly because (a) colors do not have a natural order, and
#' (b) rainbow palettes are not perceptually linear. It would be better to use a
#' grayscale palette or the 'cool-to-warm' scheme. The examples show how to do
#' this in R.
#'
#' The final spectrogram is displayed in log energy scale and by convention has
#' low frequencies on the bottom of the image.
#'
#' @param x Input signal, specified as a vector.
#' @param n Size of the FFT window. Default: 256 (or less if \code{x} is
#'   shorter).
#' @param fs Sample rate in Hz. Default: 2
#' @param window Either an integer indicating the length of a Hanning window, or
#'   a vector of values representing the shape of the FFT tapering window.
#'   Default: hanning(n)
#' @param overlap Overlap with previous window. Default: half the window length
#' @param col Colormap to use for plotting. Default: \code{grDevices::gray(0:512
#'   / 512)}
#' @param xlab Label for x-axis of plot. Default: \code{"Time"}
#' @param ylab Label for y-axis of plot. Default: \code{"Frequency"}
#' @param ... Additional arguments passed to the \code{image} plotting function
#'
#' @return A list of class \code{specgram} consisting of the following elements:
#' \describe{
#'   \item{S}{the complex output of the FFT, one row per slice}
#'   \item{f}{the frequency indices corresponding to the rows of S}
#'   \item{t}{the time indices corresponding to the columns of S}
#' }
#'
#' @examples
#'
#' sp <- specgram(chirp(seq(-2, 15, by = 0.001), 400, 10, 100, 'quadratic'))
#' specgram(chirp(seq(0, 5, by = 1/8000), 200, 2, 500,
#'         "logarithmic"), fs = 8000)
#'
#' # use other color palettes than grayscale
#' jet <- grDevices::colorRampPalette(
#'          c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F",
#'            "yellow", "#FF7F00", "red", "#7F0000"))
#' plot(specgram(chirp(seq(0, 5, by = 1/8000), 200, 2, 500, "logarithmic"),
#'          fs = 8000), col = jet(20))
#' c2w <- grDevices::colorRampPalette(colors = c("red", "white", "blue"))
#' plot(specgram(chirp(seq(0, 5, by = 1/8000), 200, 2, 500, "logarithmic"),
#'          fs = 8000), col = c2w(50))
#'
#' @author Paul Kienzle, \email{pkienzle@@users.sf.net}.\cr
#' Conversion to R by Tom Short\cr
#' This conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#'
#' @rdname specgram
#' @export

specgram <- function(x, n = min(256, length(x)), fs = 2, window = hanning(n),
                     overlap = ceiling(n / 2)) {

  if (!is.numeric(x) || !is.vector(x)) stop("x must be a numeric vector")
  if (!isPosscal(n) || !isWhole(n)) stop("n must be a positive integer")
  if (n > length(x)) {
    n <- length(x)
    warning(paste("FFT segment size adjusted to", n))
  }

  ## if only the window length is given, generate hanning window
  if (isScalar(window)) {
    if (!isPosscal(window) || !isWhole(window))
      stop("if window is a scalar, it must be a positive integer")
    window <- hanning(window)
  }

  if (!isPosscal(overlap) || !isWhole(overlap))
    stop("if window is a scalar, it must be a positive integer")

  ## compute window offsets
  win_size <- length(window)
  if (win_size > n) {
    n <- win_size
    warning(paste("FFT segment size adjusted to", n))
  }
  if (overlap >= n) stop("overlap should be smaller than n")
  step <- win_size - overlap

  ## build matrix of windowed data slices
  if (length(x) > win_size) {
    offset <- seq(1, length(x) - win_size, step)
  } else {
    offset <- 1
  }
  S <- matrix(0L, n, length(offset))
  for (i in seq_along(offset)) {
    S[1:win_size, i] <- x[offset[i]:(offset[i] + win_size - 1)] * window
  }

  ## compute Fourier transform
  S <- stats::mvfft(S)

  ## extract the positive frequency components
  if (n %% 2 == 1) {
    ret_n <- (n + 1) / 2
  } else {
    ret_n <- n / 2
  }
  S <- S[1:ret_n, ]

  f <- (0:(ret_n - 1)) * fs / n
  t <- offset / fs

  ret <- list(S = S, f = f, t = t)
  class(ret) <- "specgram"
  ret
}

#' @rdname specgram
#' @export

print.specgram <- function(x, col = grDevices::gray(0:512 / 512),
                           xlab = "Time", ylab = "Frequency", ...) {
  graphics::image(x$t, x$f, 20 * log10(t(abs(x$S))),
                  col = col, xlab = xlab, ylab = ylab, ...)
}

#' @rdname specgram
#' @export

plot.specgram <- function(x, col = grDevices::gray(0:512 / 512),
                           xlab = "Time", ylab = "Frequency", ...) {
  graphics::image(x$t, x$f, 20 * log10(t(abs(x$S))),
                  col = col, xlab = xlab, ylab = ylab, ...)
}
