# zerocrossing.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Original Octave code:
# Copyright (C) 2008 Carlo de Falco <carlo.defalco@gmail.com>
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#
# Version history
# 20201128  GvB       setup for gsignal v0.1.0
#------------------------------------------------------------------------------

#' Zero Crossing
#'
#' Estimate zero crossing points of waveform.
#'
#' @param x the x-coordinates of points in the function.
#' @param y the y-coordinates of points in the function.
#'
#' @return Zero-crossing points
#'
#' @examples
#' x <- seq(0, 1, length.out = 100)
#' y <- runif(100) - 0.5
#' x0 <- zerocrossing(x, y)
#' plot(x, y, type ="l", xlab = "", ylab = "")
#' points(x0, rep(0, length(x0)), col = "red")
#'
#' @author Carlo de Falco, \email{carlo.defalco@@gmail.com}.\cr
#' Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#
#' @export

zerocrossing <- function(x, y) {

  if (!is.vector(x) || !is.vector(y) || !is.numeric(x) || !is.numeric(y)) {
    stop("x and y must be numeric vectors")
  }
  len <- length(x)
  if (length(y) != len) {
    stop("x and y must be vectors of the same length")
  }
  crossing_intervals <- (y[1:(len - 1)] * y[2:len] <= 0)
  left_ends          <- (x[1:(len - 1)])[crossing_intervals]
  right_ends         <- (x[2:len])[crossing_intervals]
  left_vals          <- (y[1:(len - 1)])[crossing_intervals]
  right_vals         <- (y[2:len])[crossing_intervals]
  mid_points         <- (left_ends + right_ends) / 2
  zero_intervals     <- which(left_vals == right_vals)
  retval1            <- mid_points[zero_intervals]
  left_ends          <- left_ends[!left_ends %in% zero_intervals]
  right_ends         <- right_ends[!right_ends %in% zero_intervals]
  left_vals          <- left_vals[!left_vals %in% zero_intervals]
  right_vals         <- right_vals[!right_vals %in% zero_intervals]
  retval2            <- left_ends - (right_ends - left_ends) *
    left_vals / (right_vals - left_vals)
  retval             <- union(retval1, retval2)

  retval
}
