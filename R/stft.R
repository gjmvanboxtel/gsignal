# stft.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Original Octave code:
# Copyright (C) 1995-2019 Andreas Weingessel
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
# 20201206  GvB       setup for gsignal v0.1.0
# 20210411  GvB       v0.3.0 bugfix in output time points
#------------------------------------------------------------------------------

#' Short-Term Fourier Transform
#'
#' Compute the short-term Fourier transform of a vector or matrix.
#'
#' @param x input data, specified as a numeric or complex vector or matrix. In
#'   case of a vector it represents a single signal; in case of a matrix each
#'   column is a signal.
#' @param window If \code{window} is a vector, each segment has the same length
#'    as \code{window} and is multiplied by \code{window} before (optional)
#'    zero-padding and calculation of its periodogram. If \code{window} is a
#'    scalar, each segment has a length of \code{window} and a Hamming window is
#'    used. Default: \code{nextpow2(sqrt(length(x)))} (the square root of the
#'    length of \code{x} rounded up to the next power of two). The window length
#'    must be larger than 3.
#' @param overlap segment overlap, specified as a numeric value expressed as a
#'   multiple of window or segment length. 0 <= overlap < 1. Default: 0.5.
#' @param nfft Length of FFT, specified as an integer scalar. The default is the
#'   length of the \code{window} vector or has the same value as the scalar
#'   \code{window} argument.  If \code{nfft} is larger than the segment length,
#'   (seg_len), the data segment is padded \code{nfft - seg_len} zeros. The
#'   default is no padding. Nfft values smaller than the length of the data
#'   segment (or window) are ignored. Note that the use of padding to increase
#'   the frequency resolution of the spectral estimate is controversial.
#' @param fs sampling frequency (Hertz), specified as a positive scalar.
#'   Default: 1.
#'
#' @return A list containing the following elements:
#'   \describe{
#'     \item{\code{f}}{vector of frequencies at which the STFT is estimated.
#'     If \code{x} is numeric, power from negative frequencies is added to the
#'     positive side of the spectrum, but not at zero or Nyquist (fs/2)
#'     frequencies. This keeps power equal in time and spectral domains. If
#'     \code{x} is complex, then the whole frequency range is returned.}
#'     \item{\code{t}}{vector of time points at which the STFT is estimated.}
#'     \item{\code{s}}{Short-time Fourier transform, returned as a matrix or
#'     a 3-D array. Time increases across the columns of \code{s} and frequency
#'     increases down the rows. The third dimension, if present, corresponds to
#'     the input channels.}
#'   }
#'
#' @examples
#' fs <- 8000
#' y <- chirp(seq(0, 5 - 1/fs, by = 1/fs), 200, 2, 500, "logarithmic")
#' ft <- stft (y, fs = fs)
#' filled.contour(ft$t, ft$f, t(ft$s), xlab = "Time (s)",
#'                ylab = "Frequency (Hz)")
#'
#' @author Andreas Weingessel, \email{Andreas.Weingessel@@ci.tuwien.ac.at}.\cr
#' Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#'
#' @export

stft <- function(x, window = nextpow2(sqrt(NROW(x))), overlap = 0.75,
                    nfft = ifelse(isScalar(window), window, length(window)),
                    fs = 1) {

  # check parameters
  if (!(is.vector(x) || is.matrix(x)) && !(is.numeric(x) || is.complex(x))) {
    stop("x must be a numeric or complex vector or matrix")
  }

  if (is.vector(x)) {
    x <- as.matrix(x, ncol = 1)
  }
  nr <- nrow(x)
  nc <- ncol(x)

  if (!is.vector(window) || !is.numeric(window)) {
    stop("window must be a numeric vector or scalar")
  } else {
    if (isPosscal(window)) {
      if (window <= 3) {
        stop("window must be a scalar > 3 or a vector with length > 3")
      } else {
        window <- hamming(window)
      }
    } else if (length(window) <= 3) {
      stop("window must be a scalar > 3 or a vector with length > 3")
    }
  }

  if (!isScalar(overlap) || !(overlap >= 0 && overlap < 1)) {
    stop("overlap must be a numeric value >= 0 and < 1")
  }

  if (!isPosscal(nfft) || !isWhole(nfft)) {
    stop("nfft must be a positive integer")
  }

  if (!isPosscal(fs) || fs <= 0) {
    stop("fs must be a numeric value > 0")
  }

  # initialize variables
  seg_len <- length(window)
  overlap <- trunc(seg_len * overlap)
  nfft <- nextpow2(max(nfft, seg_len))
  win_meansq <- as.vector(window %*% window / seg_len)
  num_win <- floor((nr - overlap) / (seg_len - overlap))
  if (overlap >= seg_len) {
    stop("overlap must be smaller than window length")
  }
  # Pad data with zeros if shorter than segment. This should not happen.
  if (nr < seg_len) {
    x <- c(x, rep(0, seg_len - nr))
    nr <- seg_len
  }

  # MAIN CALCULATIONS

  # Calculate and accumulate periodograms
  Pxx <- array(0, dim = c(nfft, num_win, nc))
  t <- rep(0, num_win)
  n_ffts <- 0
  for (start_seg in seq(1, nr - seg_len + 1, seg_len - overlap)) {
    end_seg <- start_seg + seg_len - 1
    xx <- window * x[start_seg:end_seg, ]
    if (nc == 1) {
      fft_x <- stats::fft(xx)
    } else {
      fft_x <- stats::mvfft(xx)
    }
    # accumulate periodogram
    n_ffts <- n_ffts + 1
    Pxx[1:nfft, n_ffts, 1:nc] <- Re(fft_x * Conj(fft_x))
    t[n_ffts] <- 0.5 * (end_seg - start_seg + 1) / fs
  }

  # Convert two-sided spectra to one-sided spectra if the input is numeric
  # For one-sided spectra, contributions from negative frequencies are added
  # to the positive side of the spectrum -- but not at zero or Nyquist
  # (half sampling) frequencies.  This keeps power equal in time and spectral
  # domains, as required by Parseval theorem.
  if (is.numeric(x)) {
    if (nfft %% 2 == 0) {    # one-sided, nfft is even
      psd_len <- nfft / 2 + 1
      Pxx <- apply(Pxx, c(2, 3),
                   function(x) x[1:psd_len] +
                     c(0, x[seq(nfft, psd_len + 1, -1)], 0))
    } else {                    # one-sided, nfft is odd
      psd_len <- (nfft + 1) / 2
      Pxx <- apply(Pxx, c(2, 3),
                   function(x) x[1:psd_len] +
                     c(0, x[seq(nfft, psd_len + 1, -1)]))
    }
  } else {
    psd_len <- nfft
  }
  # end MAIN CALCULATIONS

  ## SCALING AND OUTPUT
  scale <- n_ffts * seg_len * fs * win_meansq
  s <- Pxx / scale
  f <- seq(0, psd_len - 1) * (fs / nfft)
  t <- seq(0, (nr / fs) - 1 / fs, length.out = num_win)

  if (nc == 1) {
    s <- s[, , 1]
  }

  list(f = f, t = t, s = s)
}
