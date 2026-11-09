# zp2ss.R
# Copyright (C) 2026 Geert van Boxtel <gjmvanboxtel@gmail.com>
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
# 20261109  GvB       setup for gsignal v0.4.0
#------------------------------------------------------------------------------

#' Zero-pole-gain to state-space form
#'
#' Convert zero-pole-gain filter parameters to state-space form.
#'
#' @param z complex vector of the zeros of the model
#' @param p complex vector of the poles of the model
#' @param g overall gain. Default: 1.
#'
#' @return A list with the following list elements:
#' \describe{
#'   \item{A}{State matrix}
#'   \item{B}{Input-to-state matrix}
#'   \item{C}{State-to-output matrix}
#'   \item{D}{Feedthrough matrix}
#' }
#'
#' @examples
#' g <- 1
#' z <- c(0, 0)
#' p <- pracma::roots(c(1, 0.01, 1))
#' ss <- zp2ss(z, p, g)
#' 
#' @author Geert van Boxtel, \email{gjmvanboxtel@@gmail.com}.
#'
#' @export

zp2ss <- function(z, p, g) {

  ba <- zp2tf(z, p, g)
  tf2ss(ba$b, ba$a)
}
