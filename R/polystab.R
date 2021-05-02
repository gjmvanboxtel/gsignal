# polystab.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Original Octave function:
# Copyright (C) 2001 Paul Kienzle <pkienzle@users.sf.net>
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
# 20200623  GvB       setup for gsignal v0.1.0
#------------------------------------------------------------------------------

#' Stabilize polynomial
#'
#' Stabilize the polynomial transfer function by replacing all roots outside the
#' unit circle with their reflection inside the unit circle.
#'
#' @param a vector of polynomial coefficients, normally in the z-domain
#'
#' @return Vector of stabilized polynomial coefficients.
#'
#' @examples
#' unstable <- c(-0.5, 1)
#' zplane(unstable, 1)
#' stable <- polystab(unstable)
#' zplane(stable, 1)
#'
#' @author Paul Kienzle, \email{pkienzle@@users.sf.net}.\cr Conversion to R by
#'   Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}
#
#' @export

polystab <- function(a) {

  r <- pracma::roots(a)
  v <- which(abs(r) > 1)
  if (length(v)) {
    r[v] <- 1 / Conj(r[v])
    b <- a[1] * poly(r)
    if (is.numeric(a)) {
      b <- Re(b)
    }
  } else {
    b <- a
  }
  b
}
