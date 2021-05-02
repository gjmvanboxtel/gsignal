# cceps.R
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
# 20200828  GvB       setup for gsignal v0.1.0
#------------------------------------------------------------------------------

#' Complex cepstral analysis
#'
#' Return the complex cepstrum of the input vector.
#'
#' Cepstral analysis is a nonlinear signal processing technique that is applied
#' most commonly in speech and image processing, or as a tool to investigate
#' periodic structures within frequency spectra, for instance resulting from
#' echos/reflections in the signal or to the occurrence of harmonic frequencies
#' (partials, overtones).
#'
#' The cepstrum is used in many variants. Most important are the power cepstrum,
#' the complex cepstrum, and real cepstrum. The function \code{cceps} implements
#' the complex cepstrum by computing the inverse of the log-transformed FFT,
#' i.e.,
#'
#' \deqn{cceps(x) <- ifft(log(fft(x)))}
#'
#' However, because taking the logarithm of a complex number can lead to
#' unexpected results, the phase of \code{fft(x)} needs to be unwrapped before
#' taking the log.
#'
#' @note This function returns slightly different results in comparison with the
#'   'Matlab' and 'Octave' equivalents. The 'Octave' version does not apply phase
#'   unwrapping, but has an optional correction procedure in case of zero phase
#'   at \eqn{\pi} radians. The present implementation does apply phase
#'   unwrapping so that the correction procedure is unnecessary. The 'Matlab'
#'   implementation also applies phase unwrapping, and a circular shift if
#'   necessary to avoid zero phase at \eqn{\pi} radians. The circular shift is
#'   not done here. In addition, the 'Octave' version shifts the zero frequency to
#'   the center of the series, which neither the 'Matlab' nor the present
#'   implementation do.
#'
#' @param x input data, specified as a real vector.
#'
#' @return Complex cepstrum, returned as a vector.
#'
#' @examples
#' ## Generate a sine of frequency 45 Hz, sampled at 100 Hz.
#' fs <- 100
#' t <- seq(0, 1.27, 1/fs)
#' s1 <- sin(2 * pi * 45 * t)
#' ## Add an echo with half the amplitude and 0.2 s later.
#' s2 <- s1 + 0.5 * c(rep(0L, 20), s1[1:108])
#' ## Compute the complex cepstrum of the signal. Notice the echo at 0.2 s.
#' cep <- cceps(s2)
#' plot(t, cep, type="l")
#'
#' @references \url{https://en.wikipedia.org/wiki/Cepstrum}
#'
#' @seealso \code{\link{rceps}}
#'
#' @author Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#
#' @export

cceps <- function(x) {

  if (!is.vector(x) || !is.numeric(x)) {
    stop("x must be a numeric vector")
  }

  X <- stats::fft(x)
  if (min(abs(X)) == 0) {
    stop("signal has Fourier coefficients equal to 0")
  }
  uw <- unwrap(Arg(X))
  logX <- complex(real = log(Mod(X)), imaginary = uw)
  y <- Re(ifft(logX))
  y
}
