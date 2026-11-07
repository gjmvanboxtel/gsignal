# ss2sos.R
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

#' state-space to second-order sections form
#'
#' Convert digital filter state-space to second-order section form.
#'
#' @param A State matrix (must be square)
#' @param B Input-to-state matrix (must have dimensions nrow(A) x 1)
#' @param C State-to-output matrix (must have nrow(A) columns)
#' @param D Feedthrough matrix (must have dimensions nrow(C) x 1)
#' @param tol tolerance to remove insignificant numerical noise (default: 1e-10)
#'
#' @return A list with the following list elements:
#' \describe{
#'   \item{sos}{Second-order section representation, specified as an nrow-by-6
#'   matrix, whose rows contain the numerator and denominator coefficients of
#'   the second-order sections:\cr \code{sos <- rbind(cbind(B1, A1), cbind(...),
#'   cbind(Bn, An))}, where \code{B1 <- c(b0, b1, b2)}, and \code{A1 <- c(a0,
#'   a1, a2)} for section 1, etc. The b0 entry must be nonzero for each
#'   section.}
#' \item{g}{Overall gain factor that effectively scales the output \code{b}
#'   vector (or any one of the input \code{Bi} vectors).}
#' }
#'
#' @examples
#' A <- rbind(
#'   c(0, 0, 0, 0, -0.9^5),
#'   c(1, 0, 0, 0, 0),
#'   c(0, 1, 0, 0, 0),
#'   c(0, 0, 1, 0, 0),
#'   c(0, 0, 0, 1, 0)
#' )
#' B <- matrix(c(1, 0, 0, 0, 0), ncol = 1)
#' C <- matrix(c(0, 1, 0, 0, 0.3^3), nrow = 1)
#' D <- 0
#' ss2sos(A, B, C, D)
#'
#' @author Geert van Boxtel, \email{gjmvanboxtel@@gmail.com}.
#'
#' @export

ss2sos <- function(A, B, C, D, tol = 1e-10) {

  ba <- ss2tf(A, B, C, D, tol)
  tf2sos(ba$b, ba$a)
}
