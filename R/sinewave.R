# sinewave.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Octave signal package:
# Copyright (C) 1995-2019 Friedrich Leisch
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
# 20201201 Geert van Boxtel          First version for v0.1.0
#------------------------------------------------------------------------------

#' Sine wave
#'
#' Generate a discrete sine wave.
#'
#' @param m desired length of the generated series, specified as a positive
#'   integer.
#' @param n rate, of the generated series, specified as a positive integer.
#'   Default: \code{m}.
#' @param d delay, specified as a positive integer. Default: 0.
#'
#' @return Sine wave, returned as a vector of length \code{m}.
#'
#' @examples
#' plot(sinewave(100, 10), type = "l")
#'
#' @author Friedrich Leisch.\cr
#' Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#
#' @export

sinewave <- function(m, n = m, d = 0) {

  if (!isPosscal(m) || !isPosscal(n) || !isPosscal(d)) {
    stop("m, n, and d must be positive scalars")
  }
  y <- sin((seq(1, m) + d - 1) * 2 * pi / n)
  y
}
