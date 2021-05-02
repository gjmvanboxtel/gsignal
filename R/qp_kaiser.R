# qp_kaiser.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Octave signal package:
# Copyright (C) 2002 André Carezia <andre@carezia.eng.br>
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
# 20200710 Geert van Boxtel           First version for v0.1.0
#------------------------------------------------------------------------------

#' Kaiser FIR filter design
#'
#' Compute FIR filter for use with a quasi-perfect reconstruction
#' polyphase-network filter bank.
#'
#' @param nb number of frequency bands, specified as a scalar
#' @param at attenuation (in dB) in the stop band.
#' @param linear logical, indicating linear scaling. If FALSE (default), the
#'   Kaiser window is multiplied by the ideal impulse response \eqn{h(n) = a
#'   sinc(an)} and converted to its minimum-phase version by means of a Hilbert
#'   transform.
#'
#' @return The FIR filter coefficients, of class \code{Ma}.
#'
#' @examples
#'\donttest{
#' freqz(qp_kaiser(1, 20))
#' freqz(qp_kaiser(1, 40))
#'}
#' 
#' @seealso \code{\link{Ma}}, \code{\link{filter}}, \code{\link{fftfilt}},
#'   \code{\link{fir2}}
#'
#' @author André Carezia, \email{andre@@carezia.eng.br}.\cr
#'   Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#'
#' @export

qp_kaiser <- function(nb, at, linear = FALSE) {

  if (!(isPosscal(nb) && isWhole(nb) && nb > 0)) {
    stop("nb must be a positive integer")
  }
  if (!(isPosscal(at) && at > 0)) {
    stop("at must be a positive scalar")
  }
  if (!is.logical(linear)) {
    stop("linear must be logical")
  }

  ## Bandwidth
  bandwidth <- pi / nb

  ## Attenuation correction (empirically
  ## determined by M. Gerken
  ## <mgk@lcs.poli.usp.br>)
  corr <- (1.4 + 0.6 * (at - 20) / 80) ^ (20 / at)
  at <- corr * at

  ## size of window (rounded to next odd integer)
  N <- (at - 8) / (2.285 * bandwidth)
  M <- trunc(N / 2)
  N <- 2 * M + 1

  ## Kaiser window
  if (at > 50) {
    beta <- 0.1102 * (at - 8.7)
  } else if (at > 21) {
    beta <- 0.5842 * (at - 21)^0.4 + 0.07886 * (at - 21)
  } else {
    beta <- 0
  }
  w <- kaiser(N, beta)
  ## squared in freq. domain
  wsquared <- conv(w, w)

  ## multiplied by ideal lowpass filter
  n <- - (N - 1):(N - 1)
  hideal <- 1 / nb * sinc(n / nb)
  hcomp <- wsquared %o% hideal

  ## extract square-root of response and
  ## compute minimum-phase version
  Ndft <- 2^15
  Hsqr <- sqrt(abs(stats::mvfft(postpad(hcomp, Ndft))))
  if (linear) {
    h <- Re(ifft(Hsqr))
    h <- h[2:N]
    h <- c(rev(h), h[1], h)
  } else {
    # prevent log of 0 or negative number
    Hsqr[which(Hsqr <= 0)] <- min(Hsqr[which(Hsqr > 0)])
    Hmin <- Hsqr * exp(-1i * Im(hilbert(log(Hsqr))))
    h <- Re(imvfft(Hmin))
    h <- h[1:N]
  }

  ## truncate and fix amplitude scale (H(0)=1)
  h <- h / sum(h)

  Ma(h)
}
