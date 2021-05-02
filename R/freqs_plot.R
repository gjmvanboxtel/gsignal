# freqs_plot.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Original Octave function:
# Copyright (C) 2003 Julius O. Smith III <jos@ccrma.stanford.edu>
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
# 20200413  GvB       setup for gsignal v0.1.0
# 20210423  GvB       on.exit(par(op))
#------------------------------------------------------------------------------

#' Plot frequency response
#'
#' Plot the s-plane frequency response of an IIR filter.
#'
#' @param w angular frequencies, specified as a positive real vector expressed
#'   in rad/second.
#' @param h Frequency response, specified as a complex vector.
#' @param ... additional arguments passed to the \code{plot() function}
#' 
#' @return No return value
#'
#' @examples
#' b <- c(1, 2); a <- c(1, 1)
#' w <- seq(0.01, 4, length.out = 128)
#' h <- freqs (b, a, w, plot = FALSE)
#' freqs_plot(w, h)
#' freqs_plot(w, h, log = "x")
#'
#' @author Julius O. Smith III, \email{jos@@ccrma.stanford.edu}.\cr
#' Conversion to R by Geert van Boxtel \email{gjmvanboxtel@@gmail.com}
#'
#' @export

freqs_plot <- function(w, h, ...) {

  mag <- 20 * log10(abs(h))
  phase <- unwrap(Arg(h))

  op <- graphics::par(mfrow = c(2, 1))
  on.exit(graphics::par(op))

  graphics::plot(w, mag, type = "l", xlab = "", ylab = "dB", ...)
  graphics::legend("topright", "Magnitude (dB)", lty = 1)
  graphics::title("Frequency response plot by freqs")

  graphics::plot(w, phase / (2 * pi), type = "l",
                 xlab = "Frequency (rad/s)", ylab = "Phase", ...)
  graphics::legend("topright", "Phase (radians / 2 pi)", lty = 1)
  graphics::title("")

}
