# square.R
# Copyright (C) 2019 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Octave signal package:
# Copyright (C) 2006 Paul Kienzle
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

#' Square wave
#'
#' Generate a square wave of period \eqn{2\pi} with limits +1 and -1.
#'
#' \code{y <- square(t)} generates a square wave with period \eqn{2\pi} for the
#' elements of the time array \code{t}.
#' \code{square} is similar to the sine function but creates a square wave with
#' values of –1 and 1.
#'
#' \code{y <- square(t, duty)} generates a square wave with specified duty cycle
#' \code{duty}. The duty cycle is the percent of the signal period in which the
#' square wave is positive.
#' \if{latex}{
#'   \deqn{duty cycle = \frac{ontime * 100}{ontime + offtime}}
#' }
#' \if{html}{\preformatted{
#'                      ontime * 100
#'      duty cycle =  ----------------
#'                    ontime + offtime
#' }}
#'
#' @param t Time array, specified as a vector.
#' @param duty Duty cycle, specified as a real scalar from 0 to 100. Default:
#'   50.
#'
#' @return Square wave, returned as a vector.
#'
#' @examples
#'
#' ## Create a vector of 100 equally spaced numbers from 0 to 3pi.
#' ## Generate a square wave with a period of 2pi.
#' t <- seq(0, 3*pi, length.out = 100)
#' y <- square(t)
#' plot(t/pi, y, type="l", xlab = expression(t/pi), ylab = "")
#' lines (t/pi, sin(t), col = "red")
#'
#' ## Generate a 30 Hz square wave sampled at 1 kHz for 70 ms.
#' ## Specify a duty cycle of 37%.
#' ## Add white Gaussian noise with a variance of 1/100.
#' t <- seq(0, 0.07, 1/1e3)
#' y <- square(2 * pi * 30 * t, 37) + rnorm(length(t)) / 10
#' plot(t, y, type="l", xlab = "", ylab = "")
#'
#' @author Paul Kienzle.\cr
#' Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#
#' @export

square <- function(t, duty = 50) {

  if (length(t) <= 0)
    stop("t must be a vector with length > 0")
  if (!isScalar(duty) || duty < 0 || duty > 100)
    stop("width must be a scalar between 0 and 100")

  duty <- duty / 100
  t <- t / (2 * pi)
  y <- rep(1L, length(t))
  y[t - floor(t) >= duty] <- -1
  y
}
