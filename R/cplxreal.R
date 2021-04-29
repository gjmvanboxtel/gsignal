# cplxreal.R
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
# 20200327  GvB       setup for gsignal v0.1.0
# 20200331  GvB       return only positive imaginary numbers in zc
# 20200402  GvB       only test Im(v) against tol to determine whether
#                     complex or real
# 20210405  GvB       changed 'dim' argument to MARGIN
#------------------------------------------------------------------------------

#' Sort complex conjugate pairs and real
#'
#' Sort numbers into into complex-conjugate-valued and real-valued elements.
#'
#' An error is signaled if some complex numbers could not be paired and if all
#' complex numbers are not exact conjugates (to within tol). Note that here is
#' no defined order for pairs with identical real parts but differing imaginary
#' parts.
#'
#' @param z Vector, matrix, or array of complex numbers.
#' @param tol Weighting factor \code{0 < tol < 1}, which determines the
#'   tolerance of matching. Default: \code{100 * .Machine$double.eps}.
#' @param MARGIN Vector giving the subscripts which the function will be applied
#'   over. E.g., for a matrix 1 indicates rows, 2 indicates columns, c(1, 2)
#'   indicates rows and columns. Where X has named dimnames, it can be a
#'   character vector selecting dimension names. Default: 2 (columns).
#'
#' @return A list containing two variables:
#' \describe{
#'   \item{zc}{Vector, matrix or array containing ordered complex conjugate
#'   pairs by increasing real parts. Only the positive imaginary complex numbers
#'   of each complex conjugate pair are returned.}
#'   \item{zr}{Vector, matrix or array containing ordered real numbers.}
#' }
#'
#' @examples
#' r <- cplxreal(c(1, 1 + 3i, 2 - 5i, 1-3i, 2 + 5i, 4, 3))
#'
#' @author Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#'
#' @seealso \code{\link{cplxpair}}
#'
#' @export

cplxreal <- function(z, tol = 100 * .Machine$double.eps, MARGIN = 2) {

  y <- cplxpair(z, tol, MARGIN)

  getr <- function(v) {
    lv <- length(v)
    ix <- rep(NA, lv)
    for (i in 1:lv) {
      if (abs(Im(v[i])) <= tol) ix[i] <- i
    }
    Re(v[which(!is.na(ix))])
  }
  getc <- function(v) {
    lv <- length(v)
    ix <- rep(NA, lv)
    for (i in 1:lv) {
      if (abs(Im(v[i])) > tol) ix[i] <- i
    }
    v <- v[which(!is.na(ix))]
    if (length(v)) {
      v <- v[seq(2, length(v), 2)]  # only pos imag numbers
    } else {
      v <- NULL
    }
    v
  }

  if (is.vector(y)) {
    zc <- getc(y)
    zr <- getr(y)
  } else {
    zc <- apply(y, MARGIN, getc)
    zr <- apply(y, MARGIN, getr)
  }

  list(zc = zc, zr = zr)
}
