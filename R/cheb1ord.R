# cheb1ord.R
# Copyright (C) 2019 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Octave signal package:
# Copyright (C) 2000 Paul Kienzle
# Copyright (C) 2000 Laurent S. Mazet
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
# 20200517 Geert van Boxtel          First version for v0.1.0
# 20200708 GvB                       renamed IIRfspec to FilterSpecs
#------------------------------------------------------------------------------

#' Chebyshev Type I filter order
#'
#' Compute Chebyshev type-I filter order and cutoff for the desired
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
#' @return A list of class \code{'FilterSpecs'} with the following list
#'   elements:
#' \describe{
#'   \item{n}{filter order}
#'   \item{Wc}{cutoff frequency}
#'   \item{type}{filter type, normally one of \code{"low"}, \code{"high"},
#'     \code{"stop"}, or \code{"pass"}.}
#' }

#' @examples
#' ## low-pass 30 Hz filter
#' fs <- 128
#' spec <- cheb1ord(30/(fs/2), 40/(fs/2), 0.5, 40)
#' cf <- cheby1(spec)
#' freqz(cf, fs = fs)
#'
#' @author Paul Kienzle, Laurent S. Mazet, Charles Praplan.\cr
#'  Conversion to R by Tom Short, adapted by Geert van Boxtel,
#'  \email{G.J.M.vanBoxtel@@gmail.com}.
#'
#' @seealso \code{\link{cheby1}}
#'
#' @export

cheb1ord <- function(Wp, Ws, Rp, Rs, plane = c("z", "s")) {

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

      ## Modify band edges if not symmetrical.  For a band-pass filter,
      ## the lower or upper stopband limit is moved, resulting in a smaller
      ## stopband than the caller requested.
      if (Wpw[1] * Wpw[2] < Wsw[1] * Wsw[2]) {
        Wsw[2] <- Wpw[1] * Wpw[2] / Wsw[1]
      } else {
        Wsw[1] <- Wpw[1] * Wpw[2] / Wsw[2]
      }

      w02 <- Wpw[1] * Wpw[2]
      wp <- Wpw[2] - Wpw[1]
      ws <- Wsw[2] - Wsw[1]

    ## Band-stop / band-reject / notch filter
    } else {

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
    wp <- Wsw
    ws <- Wpw

    ## Low-pass filter
  } else {
    wp <- Wpw
    ws <- Wsw
  }

  Wa <- ws / wp

  ## compute minimum n which satisfies all band edge conditions
  stop_atten <- 10 ^ (abs(Rs) / 10)
  pass_atten <- 10 ^ (abs(Rp) / 10)
  n <- ceiling(acosh(sqrt((stop_atten - 1) / (pass_atten - 1))) / acosh(Wa))

  ## compute stopband frequency limits to make the the filter characteristic
  ## touch either at least one stop band corner or one pass band corner.
  epsilon <- 1 / sqrt(10 ^ (.1 * abs(Rs)) - 1)
  k <- cosh(1 / n * acosh(sqrt(1 / (10 ^ (.1 * abs(Rp)) - 1)) / epsilon))
  # or k = fstop / fpass

  ## compute -3dB cutoff given Wp, Rp and n
  if (length(Wpw) == 2 && length(Wsw) == 2) {

    ## Band-pass filter
    if (Wpw[1] > Wsw[1]) {
      type <- "pass"
      w_prime_p <- Wpw                        #   same formula as for LP
      w_prime_s <- Wsw / k                    #           "

    ## Band-stop / band-reject / notch filter
    } else {
      type <- "stop"
      w_prime_p <- Wpw                        #   same formula as for HP
      w_prime_s <- k * Wsw                    #           "
    }

    ## freq to be returned to match pass band
    w0 <- sqrt(prod(Wpw))
    Q <- w0 / diff(Wpw)                             # BW at -Rp dB not at -3dB
    wc <- Wpw
    W_prime <- w_prime_p[1] / wc[1]                 # same with w_prime(2)/wc(2)
    wa <- abs(W_prime + sqrt(W_prime^2 + 4 * Q ^ 2)) / (2 * Q / w0)
    wb <- abs(W_prime - sqrt(W_prime^2 + 4 * Q ^ 2)) / (2 * Q / w0)
    Wcw_p <- c(wb, wa)

    ## freq to be returned to match stop band
    w0 <- sqrt(prod(Wsw))
    Q <- w0 / diff(Wsw)                             # BW at -Rs dB not at -3dB
    wc <- Wsw
    W_prime <- w_prime_s[1] / wc[1]                 # same with w_prime(2)/wc(2)
    wa <- abs(W_prime + sqrt(W_prime^2 + 4 * Q ^ 2)) / (2 * Q / w0)
    wb <- abs(W_prime - sqrt(W_prime^2 + 4 * Q ^ 2)) / (2 * Q / w0)
    Wcw_s <- c(wb, wa)

  ## High-pass filter
  } else if (Wpw > Wsw) {

    type <- "high"
    Wcw_p <- Wpw                              #   to match pass band
    Wcw_s <- Wsw * k                          #   to match stop band

  ## Low-pass filter
  } else {

    type <- "low"
    Wcw_p <- Wpw                              #   to match pass band
    Wcw_s <- Wsw / k                          #   to match stop band
  }

  if (plane == "s") {
    # No prewarp in case of analog filter
    Wc_p <- Wcw_p
    Wc_s <- Wcw_s
  } else {
    # Inverse frequency warping for discrete-time filter
    Wc_p <- atan(Wcw_p * (T / 2)) * (T / pi)
    Wc_s <- atan(Wcw_s * (T / 2)) * (T / pi)
  }

  FilterSpecs(n = n, Wc = Wc_p, type = type,
              Wc_s = Wc_s, plane = plane, Rp = Rp)
}
