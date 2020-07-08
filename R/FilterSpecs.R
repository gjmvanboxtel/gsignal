# FilterSpecs.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
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
# 20200507  GvB       setup for gsignal v0.1.0
# 20200708  GvB       renamed IIRfspec to FilterSpecs
#---------------------------------------------------------------------------------------------------------------------

#' Filter specifications
#' 
#' Filter specifications, including order, frequency cutoff, type, and
#' possibly others.
#' 
#' @param n filter order.
#' @param Wc cutoff frequency.
#' @param type filter type, normally one of "low", "high", "stop", or "pass".
#' @param ...	other filter description characteristics, possibly including Rp
#'   for dB of pass band ripple or Rs for dB of stop band ripple, depending on
#'   filter type (Butterworth, Chebyshev, etc.).
#' 
#' @return A list of class \code{'FilterSpecs'} with the following list elements
#'   (repeats of the input arguments):
#' \describe{
#'   \item{n}{filter order}
#'   \item{Wc}{cutoff frequency}
#'   \item{type}{filter type, normally one of "low", "high", "stop", or "pass".}
#'   \item{...}{other filter description characteristics, possibly including Rp
#'   for dB of pass band ripple or Rs for dB of stop band ripple, depending on
#'   filter type (Butterworth, Chebyshev, etc.).}
#' }
#' 
#' @seealso \code{\link{filter}}, butter and buttord, cheby1 and cheb1ord, ellip and ellipord.
#' 
#' @examples
#' filt <- FilterSpecs(3, 0.1, "low")
#' 
#' @author Tom Short \email{tshort@@eprisolutions.com}, renamed and adapted by
#'   Geert van Boxtel \email{gjmvanboxtel@@gmail.com}
#'
#' @export

FilterSpecs <- function(n, Wc, type, ...) {
  res = list(n = n, Wc = Wc, type = type, ...)
  class(res) = "FilterSpecs"
  res
}
