# polyreduce.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Original Octave function:
# Copyright (C) 1994-2017 John W. Eaton
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
# 20200606  GvB       setup for gsignal v0.1.0
#------------------------------------------------------------------------------

#' Reduce polynomial
#'
#' Reduce a polynomial coefficient vector to a minimum number of terms by
#' stripping off any leading zeros.
#'
#' @param pc vector of polynomial coefficients
#'
#' @return Vector of reduced polynomial coefficients.
#'
#' @examples
#' p <- polyreduce(c(0, 0, 1, 2, 3))
#'
#' @author Tony Richardson, \email{arichard@@stark.cc.oh.us},\cr
#'  adapted by John W. Eaton.\cr
#'   Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}
#
#' @export

polyreduce <- function(pc) {

  lpc <- length(pc)
  if (is.null(pc) || !is.vector(pc) || lpc == 0) {
    stop("'pc must be a non-empty vector")
  }

  idx <- which(pc != 0)[1]
  if (is.na(idx) || length(idx) == 0) {
    p <- 0
  } else {
    p <- pc[idx:lpc]
  }
  p
}
