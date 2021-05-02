# fftfilt.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Octave function:
# Copyright (C) 1994-2017 John W. Eaton
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
# 20200417  GvB       setup for gsignal v0.1.0
# 20200420  GvB       adapted slightly (rounding in case of whole numbers)
# 20201121  GvB       done away with FFTfilt, only method for Ma()
#------------------------------------------------------------------------------

#' FFT-based FIR filtering
#'
#' FFT-based FIR filtering using the overlap-add method.
#'
#' This function combines two important techniques to speed up filtering of long
#' signals, the overlap-add method, and FFT convolution. The overlap-add method
#' is used to break long signals into smaller segments for easier processing or
#' preventing memory problems. FFT convolution uses the overlap-add method
#' together with the Fast Fourier Transform, allowing signals to be convolved by
#' multiplying their frequency spectra. For filter kernels longer than about 64
#' points, FFT convolution is faster than standard convolution, while producing
#' exactly the same result.
#'
#' The overlap-add technique works as follows. When an \code{N} length signal is
#' convolved with a filter kernel of length \code{M}, the output signal is
#' \code{N + M - 1} samples long, i.e., the signal is expanded 'to the right'.
#' The signal is then broken into \code{k} smaller segments, and the convolution
#' of each segment with the f kernel will have a result of length \code{N / k +
#' M -1}. The individual segments are then added together. The rightmost \code{M
#' - 1} samples overlap with the leftmost \code{M - 1} samples of the next
#' segment. The overlap-add method produces exactly the same output signal as
#' direct convolution.
#'
#' FFT convolution uses the principle that multiplication in the frequency
#' domain corresponds to convolution in the time domain. The input signal is
#' transformed into the frequency domain using the FFT, multiplied by the
#' frequency response of the filter, and then transformed back into the time
#' domain using the inverse FFT. With FFT convolution, the filter kernel can be
#' made very long, with very little penalty in execution time.
#'
#' @param b moving average (Ma) coefficients of a FIR filter, specified as a
#'   vector.
#' @param x the input signal to be filtered. If x is a matrix, its columns are
#'   filtered.
#' @param n FFT length, specified as a positive integer. The FFT size must be an
#'   even power of 2 and must be greater than or equal to the length of
#'   \code{filt}. If the specified \code{n} does not meet these criteria, it is
#'   automatically adjusted to the nearest value that does. If \code{n = NULL}
#'   (default), then the overlap-add method is not used.
#'
#' @return The filtered signal, returned as a vector or matrix with the same
#'   dimensions as \code{x}.
#'
#' @examples
#' t <- seq(0, 1, len = 10000)                          # 1 second sample
#' x <- sin(2* pi * t * 2.3) + 0.25 * rnorm(length(t))  # 2.3 Hz sinusoid+noise
#' filt <- rep(0.1, 10)                                 # filter kernel
#' y1 <- filter(filt, 1, x)                             # use normal convolution
#' y2 <- fftfilt(filt, x)                               # FFT convolution
#' plot(t, x, type = "l")
#' lines(t, y1, col = "red")
#' lines(t, y2, col = "blue")
#'
#' ## use 'filter' with different classes
#' t <- seq(0, 1, len = 10000)                          # 1 second sample
#' x <- sin(2* pi * t * 2.3) + 0.25 * rnorm(length(t))  # 2.3 Hz sinusoid+noise
#' ma <- Ma(rep(0.1, 10))                               # filter kernel
#' y1 <- filter(ma, x)                                  # convulution filter
#' y2 <- fftfilt(ma, x)                                 # FFT filter
#' all.equal(y1, y2)                                    # same result
#'
#' @seealso \code{\link{filter}}
#'
#' @references \url{https://en.wikipedia.org/wiki/Overlap-add_method}.
#'
#' @author Kurt Hornik, \email{Kurt.Hornik@@wu-wien.ac.at},\cr
#'  adapted by John W. Eaton.\cr
#'  Conversion to R by Tom Short,\cr
#'  adapted by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#'
#' @rdname fftfilt
#' @export

fftfilt <- function(b, x, n = NULL) UseMethod("fftfilt")

#' @rdname fftfilt
#' @method fftfilt default
#' @export

fftfilt.default <- function(b, x, n = NULL) {

  if (!is.vector(b)) {
    stop("'b' must be a vector")
  } else {
    lb <- length(b)
  }

  if (is.vector(x)) {
    lx <- length(x)
  } else if (is.matrix(x)) {
    nrx <- nrow(x)
    ncx <- ncol(x)
  } else {
    stop("'x' must be a vector, a matrix, or a 2-D array")
  }

  if (is.null(n)) {
    ## Use FFT with the smallest power of 2 which is >= length (x) +
    ## length (b) - 1 as number of points ...
    if (is.vector(x)) {
      n <- nextpow2(lx + lb - 1)
      B <- stats::fft(postpad(b, n))
      y <- ifft(stats::fft(postpad(x, n)) * B)
    } else {
      n <- nextpow2(nrx + lb - 1)
      B <- stats::fft(postpad(b, n))
      y <- imvfft(stats::mvfft(postpad(x, n)) * replicate(ncx, B))
    }
  } else {
    ## Use overlap-add method ...
    if (!isPosscal(n) || !isWhole(n)) {
      stop("'n' must be a positive integer")
    }
    n <- nextpow2(max(n, lb))
    L <- n - lb + 1

    if (is.vector(x)) {
      R <- ceiling(lx / L)
      y <- rep(0L, lx)
      for (r in seq_len(R)) {
        lo <- (r - 1) * L + 1
        hi <- min(r * L, nrx)
        tmp <- rep(0L, n)
        tmp[1:(hi - lo + 1)] <- x[lo:hi]
        tmp <- ifft(stats::fft(postpad(tmp, n)) * B)
        hi <- min(lo + n - 1, lx)
        y[lo:hi] <- y[lo:hi] + tmp[1:(hi - lo + 1)]
      }
    } else {
      R <- ceiling(nrx / L)
      y <- matrix(0L, nrx, ncx)
      for (r in seq_len(R)) {
        lo <- (r - 1) * L + 1
        hi <- min(r * L, nrx)
        tmp <- matrix(0L, n, ncx)
        tmp[1:(hi - lo + 1), ] <- x[lo:hi, ]
        tmp <- imvfft(stats::mvfft(postpad(tmp, n)) * replicate(ncx, B))
        hi <- min(lo + n - 1, nrx)
        y[lo:hi, ] <- y[lo:hi, ] + tmp[1:(hi - lo + 1), ]
      }
    }
  }

  if (is.vector(x)) {
    y <- y[1:lx]
  } else {
    y <- y[1:nrx, ]
  }

  ## Final cleanups: if both x and b are real respectively integer, y
  ## should also be
  if (is.numeric(b) && is.numeric(x))
    y <- Re(y)
  if (!any(as.logical(b - round(b)))) {
    idx <- !any(as.logical(x - round(x)))
    y[idx] <- round(y[idx])
  }
  y
}

#' @rdname fftfilt
#' @method fftfilt Ma
#' @export

fftfilt.Ma <- function(b, x, n = NULL) {
  fftfilt.default(unclass(b), x, n)
}
