# chirp.R
# Copyright (C) 2019 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Matlab/Octave signal package:
# Copyright (C) 1999-2000 Paul Kienzle <pkienzle@users.sf.net>,
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
# 20191122 Geert van Boxtel          First version for v0.1.0
#------------------------------------------------------------------------------

#' Chirp signal
#'
#' Evaluate a chirp signal (frequency swept cosine wave).
#'
#' A chirp is a signal in which the frequency changes with time, commonly used
#' in sonar, radar, and laser. The name is a reference to the chirping sound
#' made by birds.
#'
#' The chirp can have one of three shapes:
#' \describe{
#'   \item{"linear"}{Specifies an instantaneous frequency sweep \eqn{f_i(t)}
#'   given by \eqn{f_i(t) = f_0 + \beta t}, where \eqn{\beta = (f_1 - f_0) /
#'   t_1} and the default value for \eqn{f_0} is 0. The coefficient \eqn{\beta}
#'   ensures that the desired frequency breakpoint \eqn{f_1} at time \eqn{t_1}
#'   is maintained.}
#'   \item{"quadratic"}{Specifies an instantaneous frequency sweep \eqn{f_i(t)}
#'   given by \eqn{f_i(t) = f_0 + \beta t^2}, where \eqn{\beta = (f_1 - f_0) /
#'   t_1^2} and the default value for \eqn{f_0} is 0. If \eqn{f_0 > f_1}
#'   (downsweep), the default shape is convex. If \eqn{f_0 < f_1} (upsweep), the
#'   default shape is concave.}
#'   \item{"logarithmic"}{Specifies an instantaneous frequency sweep
#'   \eqn{f_i(t)} given by \eqn{f_i(t) = f_0 \times \beta t}, where \eqn{\beta =
#'   \left( \frac {f_1}{f_0} \right) ^ \frac{1}{t1}} and the default value for
#'   \eqn{f_0} is \eqn{10^{-6}}.}
#' }
#'
#' @param t Time array, specified as a vector.
#' @param f0 Initial instantaneous frequency at time 0, specified as a positive
#'   scalar expressed in Hz. Default: 0 Hz for linear and quadratic shapes; 1e-6
#'   for logarithmic shape.
#' @param t1 Reference time, specified as a positive scalar expressed in
#'   seconds. Default: 1 sec.
#' @param f1 Instantaneous frequency at time t1, specified as a positive scalar
#'   expressed in Hz. Default: 100 Hz.
#' @param shape Sweep method, specified as \code{"linear"}, \code{"quadratic"},
#'   or \code{"logarithmic"} (see Details). Default: \code{"linear"}.
#' @param phase Initial phase, specified as a positive scalar expressed in
#'   degrees. Default: 0.

#' @return Chirp signal, returned as an array of the same length as \code{t}.
#'
#' @examples
#' # Shows linear sweep of 100 Hz/sec starting at zero for 5 sec
#' # since the sample rate is 1000 Hz, this should be a diagonal
#' # from bottom left to top right.
#' t <- seq(0, 5, 0.001)
#' y <- chirp (t)
#' specgram (y, 256, 1000)
#'
#' # Shows a quadratic chirp of 400 Hz at t=0 and 100 Hz at t=10
#' # Time goes from -2 to 15 seconds.
#' specgram(chirp(seq(-2, 15, by = 0.001), 400, 10, 100, "quadratic"))
#'
#' # Shows a logarithmic chirp of 200 Hz at t = 0 and 500 Hz at t = 2
#' # Time goes from 0 to 5 seconds at 8000 Hz.
#' specgram(chirp(seq(0, 5, by = 1/8000), 200, 2, 500, "logarithmic"),
#'          fs = 8000)
#'
#' @author Paul Kienzle, \email{pkienzle@@users.sf.net},\cr
#'   Mike Miller.\cr
#'   Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#'
#' @export

chirp <- function(t, f0, t1 = 1, f1 = 100,
                  shape = c("linear", "quadratic", "logarithmic"),
                  phase = 0) {

  shape <- match.arg(shape)

  # The default value for f0 depends on the shape
  if (missing(f0)) {
    if (shape == "logarithmic") {
      f0 <- 1e-6
    } else {
      f0 <- 0
    }
  }

  phase <- 2 * pi * phase / 360
  if (shape == "linear") {
      a <- pi * (f1 - f0) / t1
      b <- 2 * pi * f0
      y <- cos(a * t^2 + b * t + phase)
  } else if (shape == "quadratic") {
      a <- (2 / 3 * pi * (f1 - f0) / t1 / t1)
      b <- 2 * pi * f0
      y <- cos(a * t^3 + b * t + phase)
  } else if (shape == "logarithmic") {
      a <- 2 * pi * f0 * t1 / log(f1 / f0)
      x <- (f1 / f0) ^ (1 / t1)
      y <- cos(a * x^t + phase)
  } else {
    stop(paste("invalid frequency sweep shape", shape))
  }
  y
}
