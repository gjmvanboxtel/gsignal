# ss2zp.R
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

#' state-space to zero-pole-gain form
#'
#' Convert digital filter state-space filter parameters to zero-pole-gain form.
#'
#' @param A State matrix (must be square)
#' @param B Input-to-state matrix (must have dimensions nrow(A) x 1)
#' @param C State-to-output matrix (must have nrow(A) columns)
#' @param D Feedthrough matrix (must have dimensions nrow(C) x 1)
#' @param tol tolerance to remove insignificant numerical noise (default: 1e-10)
#'
#'@return A list of class Zpg with the following list elements:
#' \describe{
#'   \item{z}{complex vector of the zeros of the model (roots of \code{B(z)})}
#'   \item{p}{complex vector of the poles of the model (roots of \code{A(z)})}
#'   \item{g}{overall gain (\code{B(Inf)})}
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
#' ss2zp(A, B, C, D)
#'
#' @author Geert van Boxtel, \email{gjmvanboxtel@@gmail.com}.
#'
#' @export

ss2zp <- function(A, B, C, D, tol = 1e-10) {

  ba <- ss2tf(A, B, C, D, tol)
  tf2zp(ba$b, ba$a)
}
