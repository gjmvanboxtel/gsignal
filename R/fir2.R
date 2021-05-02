# fir2.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
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
# 20200703 Geert van Boxtel           First version for v0.1.0
#------------------------------------------------------------------------------

#' Frequency sampling-based FIR filter design
#'
#' Produce a FIR filter with arbitrary frequency response over frequency bands.
#'
#' The function linearly interpolates the desired frequency response onto a
#' dense grid and then uses the inverse Fourier transform and a Hamming window
#' to obtain the filter coefficients.
#'
#' @param n filter order (1 less than the length of the filter).
#' @param f vector of frequency points in the range from 0 to 1, where 1
#'   corresponds to the Nyquist frequency. The first point of \code{f} must be 0
#'   and the last point must be 1. \code{f} must be sorted in increasing order.
#'   Duplicate frequency points are allowed and are treated as steps in the
#'   frequency response.
#' @param m vector of the same length as \code{f} containing the desired
#'   magnitude response at each of the points specified in \code{f}.
#' @param grid_n length of ideal frequency response function. \code{grid_n}
#'   defaults to 512, and should be a power of 2 bigger than \code{n}.
#' @param ramp_n transition width for jumps in filter response (defaults to
#'   \code{grid_n / 20}). A wider ramp gives wider transitions but has better
#'   stopband characteristics.
#' @param window smoothing window. The returned filter is the same shape as the
#'   smoothing window. Default: \code{hamming(n + 1)}.
#'
#' @return The FIR filter coefficients, a vector of length \code{n + 1}, of
#'   class \code{Ma}.
#'
#' @examples
#' f <- c(0, 0.3, 0.3, 0.6, 0.6, 1)
#' m <- c(0, 0, 1, 1/2, 0, 0)
#' fh <- freqz(fir2(100, f, m))
#' op <- par(mfrow = c(1, 2))
#' plot(f, m, type = "b", ylab = "magnitude", xlab = "Frequency")
#' lines(fh$w / pi, abs(fh$h), col = "blue")
#' # plot in dB:
#' plot(f, 20*log10(m+1e-5), type = "b", ylab = "dB", xlab = "Frequency")
#' lines(fh$w / pi, 20*log10(abs(fh$h)), col = "blue")
#' par(op)
#'
#' @seealso \code{\link{Ma}}, \code{\link{filter}}, \code{\link{fftfilt}},
#'   \code{\link{fir1}}
#'
#' @author Paul Kienzle, \email{pkienzle@@users.sf.net}.\cr
#'  Conversion to R Tom Short,\cr
#'  adapted by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#'
#' @export

fir2 <- function(n, f, m, grid_n = 512, ramp_n = NULL,
                 window = hamming(n + 1)) {

  # filter length must be a scalar > 0
  if (!isPosscal(n) || !isWhole(n) || n <= 0) {
    stop("n must be an integer > 0")
  }

  ## verify frequency and magnitude vectors are reasonable
  t <- length(f)
  if (!is.vector(f) || t < 2 || f[1] != 0 ||
      f[t] != 1 || any(diff(f) < 0)) {
    stop("frequency vector f must be nondecreasing between 0 and 1")
  }
  if (t != length(m)) {
    stop("frequency vector f and magnitude vector m must be the same length")
  }

  ## find the grid spacing and ramp width
  if (!isPosscal(grid_n)) {
    stop("grid_n must be a positive scalar")
  }

  ## find the window parameter, or default to hamming
  if (length(window) != n + 1) {
    stop("window must be of length n+1")
  }

  ## ML behavior appears to always round the grid size up to a power of 2
  grid_n <- nextpow2(grid_n)

  ## make sure grid is big enough for the window
  if (2 * grid_n < n + 1) {
    grid_n <- nextpow2(n + 1)
  }

  # do this AFTER adapting grid_n (hence also the clumsy NULL argument)
  if (is.null(ramp_n)) {
    ramp_n <- floor(grid_n / 20)
  }

  ## Apply ramps to discontinuities
  if (ramp_n > 0) {
    ## remember original frequency points prior to applying ramps
    basef <- f
    basem <- m

    ## separate identical frequencies, but keep the midpoint
    idx <- which(diff(f) == 0)
    f[idx] <- f[idx] - ramp_n / grid_n / 2
    f[(idx + 1)] <- f[(idx + 1)] + ramp_n / grid_n / 2
    f <- c(f, basef[idx])

    ## make sure the grid points stay monotonic in [0,1]
    f[f < 0] <- 0
    f[f > 1] <- 1
    f <- sort(unique(c(f, basef[idx])))

    ## preserve window shape even though f may have changed
    m <- stats::approx(basef, basem, f, method = "linear", ties = "ordered")$y
  }

  ## interpolate between grid points
  grid <- stats::approx(f, m, seq(0, 1, length = grid_n + 1),
                        method = "linear", ties = "ordered")$y

  ## Transform frequency response into time response and
  ## center the response about n/2, truncating the excess
  if ((n %% 2) == 0) {
    b <- ifft(c(grid, grid[seq(grid_n, 2, by = -1)]))
    mid <- (n + 1) / 2
    b <- Re(c(b[(2 * grid_n - floor(mid) + 1):(2 * grid_n)], b[1:ceiling(mid)]))
  } else {
    ## Add zeros to interpolate by 2, then pick the odd values below.
    b <- ifft(c(grid, rep(0L, grid_n * 2), grid[seq(grid_n, 2, by = -1)]))
    b <- 2 * Re(c(b[seq(length(b) - n + 1, length(b), by = 2)],
                  b[seq(2, n + 1, by = 2)]))
  }

  ## Multiplication in the time domain is convolution in frequency,
  ## so multiply by our window now to smooth the frequency response.
  b <- b * window
  Ma(b)

}
