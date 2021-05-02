# cheby1.R
# Copyright (C) 2019 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Octave signal package:
# Copyright (C) 1999 Paul Kienzle <pkienzle@users.sf.net>
# Copyright (C) 2003 Doug Stewart <dastew@sympatico.ca>
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
# 20200519 Geert van Boxtel          First version for v0.1.0
# 20200708 GvB                       renamed IIRfspec to FilterSpecs
# 20210308 GvB                       added output parameter
#------------------------------------------------------------------------------

#' Chebyshev Type I filter design
#'
#' Compute the transfer function coefficients of a Chebyshev Type I filter.
#'
#' Chebyshev filters are analog or digital filters having a steeper roll-off
#' than Butterworth filters, and have passband ripple (type I) or stopband
#' ripple (type II).
#'
#' Because \code{cheby1} is generic, it can be extended to accept other inputs,
#' using \code{cheb1ord} to generate filter criteria for example.
#'
#' @param n filter order.
#' @param Rp dB of passband ripple.
#' @param w critical frequencies of the filter. \code{w} must be a scalar for
#'   low-pass and high-pass filters, and \code{w} must be a two-element vector
#'   c(low, high) specifying the lower and upper bands in radians/second. For
#'   digital filters, W must be between 0 and 1 where 1 is the Nyquist
#'   frequency.
#' @param type filter type, one of \code{"low"}, \code{"high"}, \code{"stop"},
#'   or \code{"pass"}.
#' @param plane "z" for a digital filter or "s" for an analog filter.
#' @param output Type of output, one of:
#' \describe{
#'   \item{"Arma"}{Autoregressive-Moving average (aka numerator/denominator, aka
#'   b/a)}
#'   \item{"Zpg"}{Zero-pole-gain format}
#'   \item{"Sos"}{Second-order sections}
#' }
#' Default is \code{"Arma"} for compatibility with the 'signal' package and the
#' 'Matlab' and 'Octave' equivalents, but \code{"Sos"} should be preferred for
#' general-purpose filtering because of numeric stability.
#' @param ... additional arguments passed to \code{cheby1}, overriding those
#'   given by \code{n} of class \code{\link{FilterSpecs}}.
#'
#' @return Depending on the value of the \code{output} parameter, a list of
#'   class \code{\link{Arma}}, \code{\link{Zpg}}, or \code{\link{Sos}}
#'   containing the filter coefficients
#'
#' @examples
#' ## compare the frequency responses of 5th-order
#' ## Butterworth and Chebyshev filters.
#' bf <- butter(5, 0.1)
#' cf <- cheby1(5, 3, 0.1)
#' bfr <- freqz(bf)
#' cfr <- freqz(cf)
#' plot(bfr$w / pi, 20 * log10(abs(bfr$h)), type = "l", ylim = c(-40, 0),
#'   xlim = c(0, .5), xlab = "Frequency", ylab = c("dB"))
#' lines(cfr$w / pi, 20 * log10(abs(cfr$h)), col = "red")
#'
#' # compare type I and type II Chebyshev filters.
#' c1fr <- freqz(cheby1(5, .5, 0.5))
#' c2fr <- freqz(cheby2(5, 20, 0.5))
#' plot(c1fr$w / pi, abs(c1fr$h), type = "l", ylim = c(0, 1),
#'   xlab = "Frequency", ylab = c("Magnitude"))
#'   lines(c2fr$w / pi, abs(c2fr$h), col = "red")
#'
#' @references \url{https://en.wikipedia.org/wiki/Chebyshev_filter}
#'
#' @seealso \code{\link{Arma}}, \code{\link{filter}}, \code{\link{butter}},
#'   \code{\link{ellip}}, \code{\link{cheb1ord}}, \code{\link{FilterSpecs}}
#'
#' @author Paul Kienzle, \email{pkienzle@@users.sf.net},\cr
#'   Doug Stewart, \email{dastew@@sympatico.ca}.\cr
#'   Conversion to R Tom Short,\cr
#'   adapted by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#'
#' @rdname cheby1
#' @export

cheby1 <- function(n, ...) UseMethod("cheby1")

#' @rdname cheby1
#' @export

cheby1.FilterSpecs <- function(n, ...)
  cheby1(n$n, n$Rp, n$Wc, n$type, n$plane, ...)

#' @rdname cheby1
#' @export

cheby1.default <- function(n, Rp, w, type = c("low", "high", "stop", "pass"),
                           plane = c("z", "s"),
                           output = c("Arma", "Zpg", "Sos"), ...) {

  # check input arguments
  type <- match.arg(type)
  plane <- match.arg(plane)
  output <- match.arg(output)
  if (!isPosscal(n) || !isWhole(n)) {
    stop("filter order n must be a positive integer")
  }
  if (!isPosscal(Rp) || !is.numeric(Rp)) {
    stop("passband ripple Rp must a non-negative scalar")
  }
  stop <- type == "stop" || type == "high"
  digital <- plane == "z"
  if (!is.vector(w) || (length(w) != 1 && length(w) != 2)) {
    stop(paste("frequency w must be specified as a vector of length 1 or 2",
               "(either w0 or c(w0, w1))"))
  }
  if ((type == "stop" || type == "pass") &&  length(w) != 2) {
    stop("w must be two elements for stop and bandpass filters")
  }
  if (digital && !all(w >= 0 & w <= 1)) {
    stop("critical frequencies w must be in the range [0 1]")
  } else if (!digital && !all(w >= 0)) {
    stop("critical frequencies w must be in the range [0 Inf]")
  }

  ## Prewarp to the band edges to s plane
  if (digital) {
    T <- 2                    # sampling frequency of 2 Hz
    w <- 2 / T * tan(pi * w / T)
  }

  ## Generate splane poles and zeros for the chebyshev type 1 filter
  epsilon <- sqrt(10 ^ (Rp / 10) - 1)
  v0 <- asinh(1 / epsilon) / n
  pole <- exp(1i * pi * seq(- (n - 1), (n - 1), by = 2) / (2 * n))
  pole <- -sinh(v0) * Re(pole) + 1i * cosh(v0) * Im(pole)
  zero <- numeric(0)

  ## compensate for amplitude at s=0
  gain <- prod(-pole)
  ## if n is even, the ripple starts low, but if n is odd the ripple
  ## starts high. We must adjust the s=0 amplitude to compensate.
  if (n %% 2 == 0) {
    gain <- gain / 10 ^ (Rp / 20)
  }
  zpg <- Zpg(z = zero, p = pole, g = gain)

  ## splane frequency transform
  zpg <- sftrans(zpg, w = w, stop = stop)

  ## Use bilinear transform to convert poles to the z plane
  if (digital) {
    zpg <- bilinear(zpg, T = T)
  }

  if (output == "Arma") {
    retval <- as.Arma(zpg)
  } else if (output == "Sos") {
    retval <- as.Sos(zpg)
  } else {
    retval <- zpg
  }

  retval

}
