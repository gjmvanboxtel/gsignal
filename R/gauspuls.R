# gauspuls.R
# Copyright (C) 2019 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Matlab/Octave signal package:
# Copyright (C) 2007 Sylvain Pelissier
# Copyright (C) 2018-2019 Mike Miller
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
# 20191124 Geert van Boxtel          First version for v0.1.0
#------------------------------------------------------------------------------

#' Gaussian-modulated sinusoidal RF pulse
#'
#' Generate a Gaussian modulated sinusoidal pulse sampled at times \code{t}.
#'
#' @param t Vector of time values at which the unit-amplitude Gaussian RF pulse
#'   is calculated.
#' @param fc Center frequency of the Gaussian-modulated sinusoidal pulses,
#'   specified as a real positive scalar expressed in Hz. Default: 1000
#' @param bw Fractional bandwidth of the Gaussian-modulated sinusoidal pulses,
#'   specified as a real positive scalar.
#'
#' @return Inphase Gaussian-modulated sinusoidal pulse, returned as a vector of
#'   unit amplitude at the times indicated by the time vector t.
#'
#' @examples
#'
#' fs <- 11025    # arbitrary sample rate
#' t <- seq(-10, 10, 1/fs)
#' yi1 <- gauspuls(t, 0.1, 1)
#' yi2 <- gauspuls(t, 0.1, 2)
#' plot(t, yi1, type="l", xlab = "Time", ylab = "Amplitude")
#' lines(t, yi2, col = "red")
#'
#' fs <- 11025  # arbitrary sample rate
#' f0 <- 100    # pulse train sample rate
#' x <- pulstran (seq(0, 4/f0, 1/fs), seq(0, 4/f0, 1/f0), "gauspuls")
#' plot (0:(length(x)-1) * 1000/fs, x, type="l",
#'       xlab = "Time (ms)", ylab = "Amplitude",
#'       main = "Gaussian pulse train at 10 ms intervals")
#'
#' @author Sylvain Pelissier, Mike Miller.\cr
#'  Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#
#' @export

gauspuls <- function(t, fc = 1e3, bw = 0.5) {

  if (!isPosscal(fc))
    stop("fc must be a non-negative real scalar")
  if (!isPosscal(bw) || bw <= 0)
    stop("bw must be a positive real scalar")

  fv <- - (bw^2 * fc^2) / (8 * log(10 ^ (-6 / 20)))
  tv <- 1 / (4 * pi^2 * fv)
  y <- exp(-t * t / (2 * tv)) * cos(2 * pi * fc * t)
  y
}
