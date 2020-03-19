# filtic.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Octave version:
# Copyright (C) 2004 David Billinghurst <David.Billinghurst@riotinto.com>
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
# 20200317    GvB       setup for gsignal v0.1.0
#---------------------------------------------------------------------------------------------------------------------

#' Initial conditions for filter
#' 
#' This function computes teh same values that would be obtained from function
#' \code{filter} given past inputs \code{x} and outputs \code{y}
#' 
#' The vectors \code{x} and \code{y} contain the most recent inputs and outputs
#' respectively, with the newest values first:
#' 
#' \deqn{x = [x(-1) x(-2) ... x(-nb)], nb = length(b)-1}
#' \deqn{y = [y(-1) y(-2) ... y(-na)], na = length(a)-a}
#' 
#' If \code{length(x) < nb} then it is zero padded. If \code{length(y) < na} then it is zero padded.
#' 
#' @param b the moving-average coefficients of an ARMA filter .
#' @param a the autoregressive (recursive) coefficients of an ARMA filter.
#' @param y output vector.
#' @param x input vector. Default: 0
#' 
#' @return Initial conditions for filter with coefficients \code{a} and
#'   \code{b}, input vector \code{x}, and output vector \code{y}, returned as a
#'   vector.
#' 
#' @examples
#' ## Simple low pass filter
#' b <- c(0.25, 0.25)
#' a <- c(1.0, -0.5)
#' ic <- filtic(b, a, 1, 1)
#' 
#' ## Simple high pass filter
#' b <- c(0.25, -0.25)
#' a <- c(1.0, 0.5)
#' ic <- filtic(b, a, 0, 1)
#' 
#' @seealso \code{\link{filter}}, \code{\link{filtfilt}}
#' 
#' @author David Billinghurst \email{David.Billinghurst@@riotinto.com}, port to
#'   R by Geert van Boxtel \email{G.J.M.vanBoxtel@@gmail.com}.
#' 
#' @export

filtic <- function(b, a, y, x = 0) {
  
  nz <- max(length(a) - 1, length(b) - 1)
  zf <- numeric(nz)

  
  # Pad arrays a and b to length nz+1 if required
  if (length(a) < (nz + 1)) {
    a <- postpad(a, nz + 1)
  }
  if (length(b) < (nz + 1)) {
    b <- postpad(b, nz + 1)
  }

  # Pad arrays x and y to length nz if required
  if (length(x) < nz) {
    x <- postpad(x, nz)
  }
  if (length(y) < nz) {
    y <- postpad(y, nz)
  }
  
  for (i in seq(nz, 1, -1)) {
    for (j in i:(nz-1)) {
      zf[j] <- b[j + 1] * x[i] - a[j + 1] * y[i] + zf[j + 1]
    }
    zf[nz] <- b[nz + 1] * x[i] - a[nz + 1] * y[i]
  }

  zf = zf / a[1]
  zf
}