# gmonopuls.R
# Copyright (C) 2019 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Matlab/Octave signal package:
# Copyright (C) 2007 Sylvain Pelissier <sylvain.pelissier@gmail.com>
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
# 20191125 Geert van Boxtel          First version for v0.1.0
#------------------------------------------------------------------------------

#' Gaussian monopulse
#'
#' Returns samples of the unit-amplitude Gaussian monopulse.
#'
#' @param t Vector of time values at which the unit-amplitude Gaussian monopulse
#'   is calculated.
#' @param fc Center frequency of the Gaussian monopulses, specified as a real
#'   positive scalar expressed in Hz. Default: 1000
#'
#' @return Samples of the Gaussian monopulse, returned as a vector of unit
#'   amplitude at the times indicated by the time vector \code{t}.
#'
#' @examples
#' fs <- 11025    # arbitrary sample rate
#' t <- seq(-10, 10, 1/fs)
#' y1 <- gmonopuls(t, 0.1)
#' y2 <- gmonopuls(t, 0.2)
#' plot(t, y1, type="l", xlab = "Time", ylab = "Amplitude")
#' lines(t, y2, col = "red")
#' legend("topright", legend = c("fc = 0.1", "fc = 0.2"),
#'        lty = 1, col = c(1, 2))
#'
#' @author Sylvain Pelissier, \email{sylvain.pelissier@@gmail.com}.\cr
#' Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#
#' @export

gmonopuls <- function(t, fc = 1e3) {

  if (!isPosscal(fc))
    stop("fc must be a non-negative real scalar")

  y <- 2 * sqrt(exp(1)) * pi * t * fc * exp(-2 * (pi * t * fc)^2)
  y
}
