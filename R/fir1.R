# fir1.R
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
# 20200704 Geert van Boxtel           First version for v0.1.0
# 20240821 Geert van Boxtel           Bugfix determining w_o
#------------------------------------------------------------------------------

#' Window-based FIR filter design
#'
#' FIR filter coefficients for a filter with the given order and frequency
#' cutoff.
#'
#' @param n filter order (1 less than the length of the filter).
#' @param w band edges, strictly increasing vector in the range c(0, 1), where 1
#'   is the Nyquist frequency. A scalar for highpass or lowpass filters, a
#'   vector pair for bandpass or bandstop, or a vector for an alternating
#'   pass/stop filter.
#' @param type character specifying filter type, one of \code{"low"} for a
#'   low-pass filter, \code{"high"} for a high-pass filter, \code{"stop"} for a
#'   stop-band (band-reject) filter, \code{"pass"} for a pass-band filter,
#'   \code{"DC-0"} for a bandpass as the first band of a multiband filter, or
#'   \code{"DC-1"} for a bandstop as the first band of a multiband filter.
#'   Default: \code{"low"}.
#' @param window smoothing window. The returned filter is the same shape as the
#'   smoothing window. Default: \code{hamming(n + 1)}.
#' @param scale whether to normalize or not. Use \code{TRUE} (default) or
#'   \code{"scale"} to set the magnitude of the center of the first passband to
#'   1, and \code{FALSE} or \code{"noscale"} to not normalize.
#'
#' @return The FIR filter coefficients, a vector of length \code{n + 1}, of
#'   class \code{Ma}.
#'
#' @references \url{https://en.wikipedia.org/wiki/Fir_filter}
#'
#' @examples
#' freqz(fir1(40, 0.3))
#' freqz(fir1(10, c(0.3, 0.5), "stop"))
#' freqz(fir1(10, c(0.3, 0.5), "pass"))
#'
#' @seealso \code{\link{Ma}}, \code{\link{filter}}, \code{\link{fftfilt}},
#'   \code{\link{fir2}}
#'
#' @author Paul Kienzle, \email{pkienzle@@users.sf.net},
#'  Conversion to R Tom Short,\cr
#'  adapted by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#'
#' @export

fir1 <- function(n, w, type = c("low", "high", "stop", "pass", "DC-0", "DC-1"),
                 window = hamming(n + 1), scale = TRUE) {

  type <- match.arg(type)
  if (!is.logical(scale)) {
    scale <- match.arg(scale, c("scale", "noscale"))
    scale <- scale == "scale"
  }
  if (is.function(window)) {
    window <- window(n + 1)
  } else if (is.character(window)) {
    window <- do.call(window, list(n + 1))
  }

  ## Assign default window, filter type and scale.
  ## If single band edge, the first band defaults to a pass band to
  ## create a lowpass filter.  If multiple band edges, the first band
  ## defaults to a stop band so that the two band case defaults to a
  ## band pass filter.  Ick.
  ftype <- tolower(type) %in% c("low", "stop", "dc-1")

  ## build response function according to fir2 requirements
  bands <- length(w) + 1
  f <- numeric(2 * bands)
  f[2 * bands] <- 1
  f[seq(2, 2 * bands - 1, by = 2)] <- w
  f[seq(3, 2 * bands - 1, by = 2)] <- w
  m <- numeric(2 * bands)
  m[seq(1, 2 * bands, by = 2)] <- (1:bands - (1 - ftype)) %% 2
  m[seq(2, 2 * bands, by = 2)] <- m[seq(1, 2 * bands, by = 2)]

  ## Increment the order if the final band is a pass band.  Something
  ## about having a nyquist frequency of zero causing problems.
  if (n %% 2 == 1 && m[2 * bands] == 1) {
    warning("n must be even for highpass and bandstop filters. Incrementing.")
    n <- n + 1
    if (is.vector(window) && is.double(window)) {
      ## End the window using interpolation
      M <- length(window)
      if (M == 1)
        window <- c(window, window)
      else
        window <- pracma::interp1(seq(0, 1, length = M),
                                  window, seq(0, 1, length = M + 1),
                                  if (M < 4) "linear" else "spline")
    }
  }

  ## compute the filter
  b <- fir2(n, f, m, 512, 2, window)

  ## normalize filter magnitude
  if (scale) {
    ## find the middle of the first band edge
    ## find the frequency of the normalizing gain
    if (m[1] == 1) {
      ## if the first band is a passband, use DC gain
      w_o <- 0
    } else if (f[4] == 1) {
      ## for a highpass filter,
      ## use the gain at half the sample frequency
      w_o <- 1
    } else{
      ## otherwise, use the gain at the center
      ## frequency of the first passband
      w_o <- f[3] + (f[4] - f[3]) / 2
    }
    ## compute |h(w_o)|^-1
    renorm <- 1 / abs(pracma::polyval(as.vector(b), exp(-1i * pi * w_o)))

    ## normalize the filter
    b <- renorm * b
  }

  Ma(b)
}
