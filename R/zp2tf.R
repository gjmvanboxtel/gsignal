# zp2tf.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Original Octave version:
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
# See also: http://www.gnu.org/licenses/gpl-2.0.txt
#
# Version history
# 2020402  GvB       setup for gsignal v0.1.0
#---------------------------------------------------------------------------------------------------------------------

#' Zero-pole-gain to transfer function
#' 
#' Convert digital filter zero-pole-gain data to transfer function form
#' 
#' @param z complex vector of the zeros of the model 
#' @param p complex vector of the poles of the model
#' @param k overall gain. Default: 1.
#' 
#' @return A list with the following list elements:
#' \describe{
#'   \item{b}{moving average (MA) polynomial coefficients}
#'   \item{a}{autoregressive (AR) polynomial coefficients}
#' }
#'  
#' @seealso See also \code{\link{filter}}
#' 
#' @examples
#' k <- 1
#' z <- c(0, 0)
#' p <- roots(c(1, 0.01, 1))
#' ba <- zp2tf(z, p, k)
#' 
#' @author Geert van Boxtel \email{gjmvanboxtel@@gmail.com}
#' 
#' @export

zp2tf <- function(z, p, k) {
  
  b <- Re(k * poly(z))
  a <- Re(poly(p))
  
  list(b = b, a = a)
}