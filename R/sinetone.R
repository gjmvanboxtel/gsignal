# sinetone.R
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

#' Sine tone
#'
#' Generate discrete sine tone.
#'
#' @param freq frequency of the tone, specified as a vector of positive numeric
#'   values. The length of \code{freq} should equal the length of the
#'   \code{ampl} vector; the shorter of the two is recycled to the longer
#'   vector.
#' @param rate sampling frequency, specified as a positive scalar. Default:
#'   8000.
#' @param sec length of the generated tone in seconds. Default: 1
#' @param ampl amplitude of the tone, specified as a vector of positive numeric
#'   values. The length of \code{ampl} should equal the length of the
#'   \code{freq} vector; the shorter of the two is recycled to the longer
#'   vector. Default: 64.
#'
#' @return Sine tone, returned as a vector of length \code{rate * sec}, or as a
#'   matrix with \code{rate * sec} columns and  \code{max(length(freq),
#'   length(ampl))} columns.
#'
#' @examples
#' fs <- 1000
#' sec <- 2
#' y <- sinetone(10, fs, sec, 1)
#' plot(seq(0, sec, length.out = sec * fs), y, type = "l", xlab = "", ylab = "")
#'
#' y <- sinetone(c(10, 15), fs, sec, c(1, 2))
#' matplot(seq(0, sec, length.out = sec * fs), y, type = "l",
#'         xlab = "", ylab = "")
#'
#' @author Friedrich Leisch.\cr
#' Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#
#' @export

sinetone <- function(freq, rate = 8000, sec = 1, ampl = 64) {

  if (!is.vector(freq) || !is.numeric(freq) || !all(freq >= 0)) {
    stop("freq must be a numeric vector > 0")
  }
  if (!is.numeric(rate) || !isScalar(rate) || rate <= 0) {
    stop("rate must be a numeric value > 0")
  }
  if (!is.numeric(sec) || !isScalar(sec) || sec <= 0) {
    stop("sec must be a numeric value > 0")
  }
  if (!is.vector(ampl) || !is.numeric(ampl)) {
    stop("freq must be a numeric vector")
  }

  lf <- length(freq)
  la <- length(ampl)
  maxl <- max(lf, la)
  if (lf < maxl) {
    freq <- rep(freq, length.out = maxl)
  }
  if (la < maxl) {
    ampl <- rep(ampl, length.out = maxl)
  }

  ns <- round(rate * sec)
  y <- matrix(0, ns, maxl)

  for (k in seq_len(maxl)) {
    y[, k] <- ampl[k] * sin(2 * pi * seq(1, ns) / rate * freq[k])
  }

  if (maxl == 1) {
    y <- as.vector(y)
  }
  y
}
