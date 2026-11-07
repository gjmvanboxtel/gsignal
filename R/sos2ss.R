# sos2ss.R
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
# 20261107  GvB       setup for gsignal v0.4.0
#------------------------------------------------------------------------------

#' Sos to state-space form
#'
#' Convert digital filter second-order section data to state-space form.
#'
#' @param sos Second-order section representation, specified as an nrow-by-6
#'   matrix, whose rows contain the numerator and denominator coefficients of
#'   the second-order sections:\cr \code{sos <- rbind(cbind(B1, A1), cbind(...),
#'   cbind(Bn, An))}, where \code{B1 <- c(b0, b1, b2)}, and \code{A1 <- c(a0,
#'   a1, a2)} for section 1, etc. The b0 entry must be nonzero for each section.
#' @param g Overall gain factor that effectively scales the output \code{b}
#'   vector (or any one of the input \code{Bi} vectors). Default: 1.
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
#' sos <- rbind(c(1, -0.3, 0.09, 1,  0.5562306, 0.81),
#'              c(1,  0.3, 0.00, 1, -1.4562306, 0.81),
#'              c(1,  0.0, 0.00, 1,  0.9000000, 0.00))
#' ss <- sos2ss(sos)
#' ss <- sos2ss(sos, 2)
#'
#' @author Geert van Boxtel, \email{gjmvanboxtel@@gmail.com}.
#'
#' @export

sos2ss <- function(sos, g = 1) {

  ba <- sos2tf(sos, g)
  tf2ss(ba$b, ba$a)
}
