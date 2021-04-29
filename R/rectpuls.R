# rectpuls.R
# Copyright (C) 2019 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Matlab/Octave signal package:
# Copyright (C) 2000 Paul Kienzle, Copyright (C) 2018-2019 Mike Miller
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

#' Rectangular pulse
#'
#' Return samples of the unit-amplitude rectangular pulse at the times
#' indicated by \code{t}.
#'
#' \code{y <- rectpuls(t)} returns a continuous, aperiodic, unit-height
#' rectangular pulse at the sample times indicated in array t, centered about t
#' = 0.
#'
#' \code{y <- rectpuls(t, w)} generates a rectangular pulse over the interval
#' from \code{-w/2} to  \code{w/2}, sampled at times \code{t}. This is useful
#' with the function \code{pulstran} for generating a series of pulses.
#'
#' @param t Sample times of unit rectangular pulse, specified by a vector.
#' @param w Rectangle width, specified by a positive number. Default: 1
#'
#' @return Rectangular pulse of unit amplitude, returned as a vector.
#'
#' @seealso \code{\link{pulstran}}
#'
#' @examples
#'
#' fs <- 10e3
#' t <- seq(-0.1, 0.1, 1/fs)
#' w <- 20e-3
#' y <- rectpuls(t, w)
#' plot(t, y, type="l", xlab = "Time", ylab = "Amplitude")
#'
#' fs <- 11025  # arbitrary sample rate
#' f0 <- 100    # pulse train sample rate
#' w <- 0.3/f0  # pulse width 1/10th the distance between pulses
#' y <- pulstran (seq(0, 4/f0, 1/fs), seq(0, 4/f0, 1/f0), 'rectpuls', w = w)
#' plot (seq(0, length(y)-1) * 1000/fs, y, type ="l", xlab = "Time (ms)",
#'       ylab = "Amplitude",
#'       main = "Rectangular pulse train of 3 ms pulses at 10 ms intervals")
#'
#' @author Paul Kienzle, Mike Miller.\cr
#' Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#
#' @export

rectpuls <- function(t, w = 1) {

  if (length(t) <= 0)
    stop("t must be a vector with length > 0")
  if (!isScalar(w) || w < 0)
    stop("w must be a positive scalar")

  y <- rep(0L, length(t))
  idx <- which((t >= -w / 2) & (t < w / 2))
  y[idx] <- 1
  y
}
