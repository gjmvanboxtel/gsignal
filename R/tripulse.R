# tripuls.R
# Copyright (C) 2019 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Octave signal package:
# Copyright (C) 2001 Paul Kienzle
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
# 20191127 Geert van Boxtel          First version for v0.1.0
#------------------------------------------------------------------------------

#' Sampled aperiodic triangle
#'
#' Generate a triangular pulse over the interval \code{-w / 2} to \code{w / 2},
#' sampled at times \code{t}.
#'
#' \code{y <- tripuls(t)} returns a continuous, aperiodic, symmetric,
#' unity-height triangular pulse at the times indicated in array \code{t},
#' centered about \code{t = 0} and with a default width of 1.
#'
#' \code{y <- tripuls(t, w)} generates a triangular pulse of width \code{w}.
#'
#' \code{y <- tripuls(t, w, skew)} generates a triangular pulse with skew
#' \code{skew}, where \eqn{-1 \le skew \le 1}. When \code{skew} is 0, a
#' symmetric triangular pulse is generated.
#'
#' @param t Sample times of triangle wave specified by a vector.
#' @param w Width of the triangular pulse to be generated. Default: 1.
#' @param skew Skew, a value between -1 and 1, indicating the relative placement
#'   of the peak within the width. -1 indicates that the peak should be at
#'   \code{-w / 2}, and 1 indicates that the peak should be at \code{w / 2}.
#'   Default: 0 (no skew).
#'
#' @return Triangular pulse, returned as a vector.
#'
#' @examples
#'
#' fs <- 10e3
#' t <- seq(-0.1, 0.1, 1/fs)
#' w <- 40e-3
#' y <- tripuls(t, w)
#' plot(t, y, type="l", xlab = "", ylab = "",
#'      main = "Symmetric triangular pulse")
#'
#' ## displace into paste and future
#' tpast <- -45e-3
#' spast <- -0.45
#' ypast <- tripuls(t-tpast, w, spast)
#' tfutr <- 60e-3
#' sfutr <- 1
#' yfutr <- tripuls(t-tfutr, w/2, sfutr)
#' plot (t, y, type = "l", xlab = "", ylab = "", ylim = c(0, 1))
#' lines(t, ypast, col = "red")
#' lines(t, yfutr, col = "blue")
#'
#' @author Paul Kienzle, Mike Miller.\cr
#' Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#
#' @export

tripuls <- function(t, w = 1, skew = 0) {

  if (length(t) <= 0) stop("t must be a vector with length > 0")
  if (!isScalar(w)) stop("w must be a scalar")
  if (!isScalar(skew) || skew < -1 || skew > 1)
    stop("skew must be a scalar between 0 and 1")

  y <- rep(0L, length(t))
  peak <- skew * w / 2

  idx <- which((t >= -w / 2) & (t <= peak))
  if (length(idx) > 0) {
    y[idx] <- (t[idx] + w / 2) / (peak + w / 2)
  }

  idx <- which((t > peak) & (t < w / 2))
  if (length(idx) > 0) {
    y[idx] <- (t[idx] - w / 2) / (peak - w / 2)
  }

  y
}
