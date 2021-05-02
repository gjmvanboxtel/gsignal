# barthannwin.R
# Copyright (C) 2019 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Octave signal package:
# Copyright (C) 1995-2017 Andreas Weingessel
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
# 20191210 Geert van Boxtel          First version for v0.1.0
#------------------------------------------------------------------------------

#' Modified Bartlett-Hann window
#'
#' Return the filter coefficients of a modified Bartlett-Hann window.
#'
#' Like Bartlett, Hann, and Hamming windows, the Bartlett-Hann window has a
#' mainlobe at the origin and asymptotically decaying sidelobes on both sides.
#' It is a linear combination of weighted Bartlett and Hann windows with near
#' sidelobes lower than both Bartlett and Hann and with far sidelobes lower than
#' both Bartlett and Hamming windows. The mainlobe width of the modified
#' Bartlett-Hann window is not increased relative to either Bartlett or Hann
#' window mainlobes.
#'
#' @param n Window length, specified as a positive integer.
#'
#' @return Modified Bartlett-Hann window, returned as a vector. If you specify a
#'   one-point window \code{(n = 1)}, the value 1 is returned.
#'
#' @examples
#'
#' t <- barthannwin(64)
#' plot (t, type = "l", xlab = "Samples", ylab =" Amplitude")
#'
#' @author Andreas Weingessel, \email{Andreas.Weingessel@@ci.tuwien.ac.at}.
#'   Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#'
#' @seealso \code{\link{bartlett}}, \code{\link{hann}}, \code{\link{hamming}}
#
#' @export

barthannwin <- function(n) {

  if (!isPosscal(n) || ! isWhole(n) || n <= 0) {
    stop("n must be an integer strictly positive")
  }

  if (n == 1) {
    w <- 1
  } else {
    N <- n - 1
    m <- 0:N
    w <- 0.62 - 0.48 * abs(m / (n - 1) - 0.5) +
      0.38 * cos(2 * pi * (m / (n - 1) - 0.5))
  }
  w
}
