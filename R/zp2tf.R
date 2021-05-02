# zp2tf.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
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
# 20200402  GvB       setup for gsignal v0.1.0
# 20200405  GvB       Set default k = 1
# 20200406  GvB       validated
# 20210326  GvB       renamed k to g, return object of class 'Arma'
#------------------------------------------------------------------------------

#' Zero-pole-gain to transfer function
#'
#' Convert digital filter zero-pole-gain data to transfer function form
#'
#' @param z complex vector of the zeros of the model
#' @param p complex vector of the poles of the model
#' @param g overall gain. Default: 1.
#'
#' @return A list of class "Arma" with the following list elements:
#' \describe{
#'   \item{b}{moving average (MA) polynomial coefficients}
#'   \item{a}{autoregressive (AR) polynomial coefficients}
#' }
#'
#' @seealso \code{\link{as.Arma}}, \code{\link{filter}}
#'
#' @examples
#' g <- 1
#' z <- c(0, 0)
#' p <- pracma::roots(c(1, 0.01, 1))
#' ba <- zp2tf(z, p, g)
#'
#' @author Geert van Boxtel, \email{gjmvanboxtel@@gmail.com}
#'
#' @export

zp2tf <- function(z, p, g = 1) {

  b <- Re(g * poly(z))
  a <- Re(poly(p))

  Arma(b = b, a = a)
}
