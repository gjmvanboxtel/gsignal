# ellipord.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Octave signal package:
# Copyright (C) 2001 Paulo Neis
# Copyright (C) 2018 Charles Praplan
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
# 20200527 Geert van Boxtel          First version for v0.1.0
# 20200708 GvB                       renamed IIRfspec to FilterSpecs
#------------------------------------------------------------------------------

#' Elliptic Filter Order
#'
#' Compute elliptic filter order and cutoff for the desired
#' response characteristics.
#'
#' @param Wp,Ws pass-band and stop-band edges. For a low-pass or high-pass
#'   filter, \code{Wp} and \code{Ws} are scalars. For a band-pass or
#'   band-rejection filter, both are vectors of length 2. For a low-pass filter,
#'   \code{Wp < Ws}. For a high-pass filter, \code{Ws > Wp}. For a band-pass
#'   \code{(Ws[1] < Wp[1] < Wp[2] < Ws[2])} or band-reject \code{(Wp[1] < Ws[1]
#'   < Ws[2] < Wp[2])} filter design, \code{Wp} gives the edges of the pass
#'   band, and \code{Ws} gives the edges of the stop band. For digital filters,
#'   frequencies are normalized to [0, 1], corresponding to the range [0, fs/2].
#'   In case of an analog filter, all frequencies are specified in radians per
#'   second.
#' @param Rp allowable decibels of ripple in the pass band.
#' @param Rs minimum attenuation in the stop band in dB.
#' @param plane "z" for a digital filter or "s" for an analog filter.
#'
#' @return A list of class \code{\link{FilterSpecs}} with the following list
#'   elements:
#' \describe{
#'   \item{n}{filter order}
#'   \item{Wc}{cutoff frequency}
#'   \item{type}{filter type, one of \code{"low"}, \code{"high"}, \code{"stop"},
#'   or \code{"pass"}.}
#'   \item{Rp}{dB of passband ripple.}
#'   \item{Rs}{dB of stopband ripple.}
#' }

#' @examples
#' fs <- 10000
#' spec <- ellipord(1000/(fs/2), 1200/(fs/2), 0.5, 29)
#' ef <- ellip(spec)
#' hf <- freqz(ef, fs = fs)
#' plot(c(0, 1000, 1000, 0, 0), c(0, 0, -0.5, -0.5, 0),
#'      type = "l", xlab = "Frequency (Hz)", ylab = "Attenuation (dB)",
#'      col = "red", ylim = c(-35,0), xlim = c(0,2000))
#' lines(c(5000, 1200, 1200, 5000, 5000), c(-1000, -1000, -29, -29, -1000),
#'       col = "red")
#' lines(hf$w, 20*log10(abs(hf$h)))
#'
#' @author Paulo Neis, \email{p_neis@@yahoo.com.br},\cr
#'   adapted by Charles Praplan.\cr
#'   Conversion to R by Tom Short,\cr
#'   adapted by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#'
#' @seealso \code{\link{buttord}}, \code{\link{cheb1ord}},
#'   \code{\link{cheb2ord}}, \code{\link{ellip}}
#'
#' @export

ellipord <- function(Wp, Ws, Rp, Rs, plane = c("z", "s")) {

  #input validation
  plane <- match.arg(plane)
  if (! (is.vector(Wp) && is.vector(Ws) && (length(Wp) == length(Ws)))) {
    stop("Wp and Ws must both be scalars or vectors of length 2")
  }
  if (! ((length(Wp) == 1) || (length(Wp) == 2))) {
    stop("Wp and Ws must both be scalars or vectors of length 2")
  }
  if (plane == "z" && !(is.numeric(Wp) && all(Wp >= 0) && all(Wp <= 1))) {
    stop("all elements of Wp must be in the range [0,1]")
  }
  if (plane == "z" && !(is.numeric(Ws) && all(Ws >= 0) && all(Ws <= 1))) {
    stop("all elements of Ws must be in the range [0,1]")
  }
  if (plane == "s" && !(is.numeric(Wp) && all(Wp >= 0))) {
    stop("all elements of Wp must be non-negative")
  }
  if (plane == "s" && !(is.numeric(Ws) && all(Ws >= 0))) {
    stop("all elements of Ws must be non-negative")
  }
  if ((length(Wp) == 2) && (Wp[2] <= Wp[1])) {
    stop("Wp[1] must be smaller than Wp[2]")
  }
  if ((length(Ws) == 2) && (Ws[2] <= Ws[1])) {
    stop("Ws[1] must be smaller than Ws[2]")
  }
  if ((length(Wp) == 2) && (all(Wp > Ws) || all(Ws > Wp))) {
    stop("Wp must be contained by Ws or Ws must be contained by Wp")
  }

  if (plane == "s") {
    # No prewarp in case of analog filter
    Wpw <- Wp
    Wsw <- Ws
  } else {
    ## sampling frequency of 2 Hz
    T <- 2
    Wpw <- (2 / T) * tan(pi * Wp / T)    # prewarp
    Wsw <- (2 / T) * tan(pi * Ws / T)    # prewarp
  }

  ## pass/stop band to low pass filter transform:
  if (length(Wpw) == 2 && length(Wsw) == 2) {

    ## Band-pass filter
    if (Wpw[1] > Wsw[1]) {

      type <- "pass"

      ## Modify band edges if not symmetrical.  For a band-pass filter,
      ## the lower or upper stopband limit is moved, resulting in a smaller
      ## stopband than the caller requested.
      if (Wpw[1] * Wpw[2] < Wsw[1] * Wsw[2]) {
        Wsw[2] <- Wpw[1] * Wpw[2] / Wsw[1]
      } else {
        Wsw[1] <- Wpw[1] * Wpw[2] / Wsw[2]
      }

      wp <- Wpw[2] - Wpw[1]
      ws <- Wsw[2] - Wsw[1]

    ## Band-stop / band-reject / notch filter
    } else {

      type <- "stop"

      ## Modify band edges if not symmetrical.  For a band-stop filter,
      ## the lower or upper passband limit is moved, resulting in a smaller
      ## rejection band than the caller requested.
      if (Wpw[1] * Wpw[2] > Wsw[1] * Wsw[2]) {
        Wpw[2] <- Wsw[1] * Wsw[2] / Wpw[1]
      } else {
        Wpw[1] <- Wsw[1] * Wsw[2] / Wpw[2]
      }

      w02 <- Wpw[1] * Wpw[2]
      wp <- w02 / (Wpw[2] - Wpw[1])
      ws <- w02 / (Wsw[2] - Wsw[1])
    }
    ws <- ws / wp
    wp <- 1

  ## High-pass filter
  } else if (Wpw > Wsw) {
    type <- "high"
    wp <- Wsw
    ws <- Wpw

    ## Low-pass filter
  } else {
    type <- "low"
    wp <- Wpw
    ws <- Wsw
  }

  k <- wp / ws
  k1 <- sqrt(1 - k^2)
  q0 <- (1 / 2) * ((1 - sqrt(k1)) / (1 + sqrt(k1)))
  q <- q0 + 2 * q0^5 + 15 * q0^9 + 150 * q0^13
  D <- (10 ^ (0.1 * Rs) - 1) / (10 ^ (0.1 * Rp) - 1)

  n <- ceiling(log10(16 * D) / log10(1 / q))

  if (plane == "s") {
    # No prewarp in case of analog filter
    Wc <- Wpw
  } else {
    # Inverse frequency warping for discrete-time filter
    Wc <- atan(Wpw * (T / 2)) * (T / pi)
  }

  FilterSpecs(n = n, Wc = Wc, type = type, plane = plane, Rp = Rp, Rs = Rs)
}
