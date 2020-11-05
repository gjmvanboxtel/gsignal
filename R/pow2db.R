# pow2db.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Original Octave code:
# Copyright (C) 2018 P Sudeepam
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
# 20201105  GvB       setup for gsignal v0.1.0
#---------------------------------------------------------------------------------------------------------------------

#' Power - decibel conversion
#' 
#' Convert power to decibel and decibel to power.
#' 
#' @param x input data, specified as a numeric vector, matrix, or
#'   multidimensional array. Must be non-negative for numeric \code{x}.
#'  
#' @return converted data, same type and dimensions as \code{x}.
#' 
#' @examples
#' pow2db(c(0, 10, 100))
#' db2pow(c(-10, 0, 10))
#' 
#' @author P Sudeepam; port to R by Geert van Boxtel,
#'   \email{G.J.M.vanBoxtel@@gmail.com}.
#' 
#' @rdname pow2db
#' @export

pow2db <- function (x) {

  if (!(is.numeric(x) || is.complex(x))) {
    stop("x must be a numeric or complex vector, matrix, or array")
  }
  if (is.numeric(x) && any(x < 0)) {
    stop("x must be non-negative")
  }
  10 * log10(x)
}

#' @rdname pow2db
#' @export

db2pow <- function (x) {
  10^(x / 10)
}
