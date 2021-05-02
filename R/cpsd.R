# cpsd.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Original Octave code:
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
# 20201104  GvB       setup for gsignal v0.1.0
#------------------------------------------------------------------------------

#' Cross power spectral density
#'
#' Estimates the cross power spectral density (CPSD) of discrete-time signals.
#'
#' \code{cpsd} estimates the cross power spectral density function using
#' Welch’s overlapped averaged periodogram method [1].
#'
#' @param x input data, specified as a numeric vector or matrix. In case of a
#'   vector it represents a single signal; in case of a matrix each column is a
#'   signal.
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
#' @param detrend character string specifying detrending option; one of:
#'   \describe{
#'     \item{\code{"long-mean"}}{remove the mean from the data before
#'     splitting into segments (default)}
#'     \item{\code{"short-mean"}}{remove the mean value of each segment}
#'     \item{\code{"long-linear"}}{remove linear trend from the data before
#'     splitting into segments}
#'     \item{\code{"short-linear"}}{remove linear trend from each segment}
#'     \item{\code{"none"}}{no detrending}
#'  }
#'
#' @return A list containing the following elements:
#'   \describe{
#'     \item{\code{freq}}{vector of frequencies at which the spectral variables
#'     are estimated. If \code{x} is numeric, power from negative frequencies is
#'     added to the positive side of the spectrum, but not at zero or Nyquist
#'     (fs/2) frequencies. This keeps power equal in time and spectral domains.
#'     If \code{x} is complex, then the whole frequency range is returned.}
#'     \item{\code{cross}}{NULL for univariate series. For multivariate series,
#'     a matrix containing the squared coherence between different series.
#'     Column \eqn{i + (j - 1) * (j - 2)/2 } of \code{coh} contains the
#'     cross-spectral estimates between columns \eqn{i} and \eqn{j} of \eqn{x},
#'     where \eqn{i < j}.}
#'   }
#'
#' @examples
#' fs <- 1000
#' f <- 250
#' t <- seq(0, 1 - 1/fs, 1/fs)
#' s1 <- sin(2 * pi * f * t) + runif(length(t))
#' s2 <- sin(2 * pi * f * t - pi / 3) + runif(length(t))
#' rv <- cpsd(cbind(s1, s2), fs = fs)
#' plot(rv$freq, 10 * log10(rv$cross), type="l", xlab = "Frequency",
#'      ylab = "Cross Spectral Density (dB)")
#'
#' @note The function \code{cpsd} (and its deprecated alias \code{csd})
#'   is a wrapper for the function \code{pwelch}, which is more complete and
#'   more flexible.
#'
#' @author Peter V. Lanspeary, \email{pvl@@mecheng.adelaide.edu.au}.\cr
#'  Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#'
#' @references [1] Welch, P.D. (1967). The use of Fast Fourier Transform for
#'   the estimation of power spectra: A method based on time averaging over
#'   short, modified periodograms. IEEE Transactions on Audio and
#'   Electroacoustics, AU-15 (2): 70–73.\cr
#'
#' @rdname cpsd
#' @export

cpsd <- function(x, window = nextpow2(sqrt(NROW(x))), overlap = 0.5,
                 nfft = ifelse(isScalar(window), window, length(window)),
                 fs = 1,
                 detrend = c("long-mean", "short-mean",
                             "long-linear", "short-linear",
                             "none")) {

  pw <- pwelch(x, window, overlap, nfft, fs, detrend)
  rv <- list(freq = pw$freq, cross = pw$cross)
  rv
}

#' @rdname cpsd
#' @export

csd <- cpsd
