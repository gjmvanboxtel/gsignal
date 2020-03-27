# cplxreal.R
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
# 20200327  GvB       setup for gsignal v0.1.0
#---------------------------------------------------------------------------------------------------------------------

#' Sort complex conjugate pairs and real
#'
#' Sort numbers into into complex-conjugate-valued and real-valued elements.
#' 
#' An error is signalled if some complex numbers could not be paired and if all
#' complex numbers are not exact conjugates (to within tol). Note that here is
#' no defined order for pairs with identical real parts but differing imaginary
#' parts.
#'
#' @param z Vector, matrix, or array of complex numbers.
#' @param tol Weighting factor \code{0 < tol < 1}, which determines the
#'   tolerance of matching. Default: 100 * .Machine$double.eps.
#' @param dim Dimension along which to sort the complex pairs. Default: 2
#'   (columns).
#'
#' @return A list containing two variables:
#' \describe{
#'   \item{zc}{Vector, matrix or array containg ordered complex conjugate pairs
#'   by increasing real parts.}
#'   \item{zr}{Vector, matrix or array containg ordered real numbers.}
#' }
#'
#' @examples
#' cplxreal(c(1, 1 + 3i, 2 - 5i, 1-3i, 2 + 5i, 4, 3))
#'
#' @author Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#'
#' @seealso \code{\link{cplxpair}}
#' 
#' @export

cplxreal <- function (z, tol = 100 * .Machine$double.eps, dim = 2) {
  
  y <- cplxpair (z, tol, dim)
  
  getr <- function (v) {
    lv <- length(v)
    ix <- rep(NA, lv)
    for (i in 1:lv) {
      if (abs(Im(v[i])) / abs(v[i]) <= tol) ix[i] <- i
    }
    Re(v[which(!is.na(ix))])
  }
  getc <- function (v) {
    lv <- length(v)
    ix <- rep(NA, lv)
    for (i in 1:lv) {
      if (abs(Im(v[i])) / abs(v[i]) > tol) ix[i] <- i
    }
    v[which(!is.na(ix))]
  }
  
  if (is.vector(y)) {
    zc <- getc(y)
    zr <- getr(y)
  } else {
    zc <- apply(y, dim, getc)
    zr <- apply(y, dim, getr)
  }
  
  list(zc = zc, zr = zr)
}
