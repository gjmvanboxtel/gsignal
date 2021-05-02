# sawtooth.R
# Copyright (C) 2019 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Octave signal package:
# Copyright (C) 2007 Juan Aguado
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
# 20191127 Geert van Boxtel          First version for v0.1.0
#------------------------------------------------------------------------------

#' Sawtooth or triangle wave
#'
#' Returns samples of the sawtooth function at the times indicated by \code{t}.
#'
#' The code \code{y <- sawtooth(t)} generates a sawtooth wave with period
#' \eqn{2\pi} for the elements of the time array \code{t}. \code{sawtooth()} is
#' similar to the sine function but creates a sawtooth wave with peaks of –1 and
#' 1. The sawtooth wave is defined to be –1 at multiples of \eqn{2\pi} and to
#' increase linearly with time with a slope of \eqn{1/\pi} at all other times.
#'
#' \code{y <- sawtooth(t, width)} generates a modified triangle wave with the
#' maximum location at each period controlled by \code{width}. Set \code{width}
#' to 0.5 to generate a standard triangle wave.
#'
#' @param t Sample times of unit sawtooth wave specified by a vector.
#' @param width Real number between 0 and 1 which specifies the point between 0
#'   and \eqn{2 \pi} where the maximum is. The function increases linearly from
#'   -1 to 1 in the interval from 0 to \eqn{ 2 * \pi * width}, and decreases
#'   linearly from 1 to -1 in the interval from \eqn{2 * \pi * width} to \eqn{2
#'   * \pi}. Default: 1 (standard sawtooth).
#'
#' @return Sawtooth wave, returned as a vector.
#'
#' @examples
#'
#' T <- 10 * (1 / 50)
#' fs <- 1000
#' t <- seq(0, T-1/fs, 1/fs)
#' y <- sawtooth(2 * pi * 50 *t)
#' plot(t, y, type="l", xlab = "", ylab = "", main = "50 Hz sawtooth wave")
#'
#' T <- 10 * (1 / 50)
#' fs <- 1000
#' t <- seq(0, T-1/fs, 1/fs)
#' y <- sawtooth(2 * pi * 50 * t, 1/2)
#' plot(t, y, type="l", xlab = "", ylab = "", main = "50 Hz triangle wave")
#'
#' @author Juan Aguado.\cr
#' Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#
#' @export

sawtooth <- function(t, width = 1) {

  if (length(t) <= 0)
    stop("t must be a vector with length > 0")
  if (!isScalar(width) || width < 0 || width > 1)
    stop("width must be a scalar between 0 and 1")

  t <- (t / (2 * pi)) %% 1
  y <- rep(0L, length(t))

  if (width != 0) {
    y[t < width] <- 2 * t[t < width] / width - 1
  }

  if (width != 1) {
    y[t >= width] <- -2 * (t[t >= width] - width) / (1 - width) + 1
  }
  y
}
