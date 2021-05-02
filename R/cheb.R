# cheb.R
# Copyright (C) 2019 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Octave signal package:
# Copyright (C) 2002 André Carezia <acarezia@uol.com.br>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; see the file COPYING. If not, see
# <https://www.gnu.org/licenses/>.
#
# 20191211 GvB          First version for v0.1.0
# 20210405 GvB          v0.3.0 corrected numerical rounding error
#                       when called by chebwin(7)
#------------------------------------------------------------------------------

#' Chebyshev polynomials
#'
#' Return the value of the Chebyshev polynomial at specific points.
#'
#' The Chebyshev polynomials are defined by the equations:
#' \if{latex}{
#'   \deqn{Tn(x) = cos(n \cdot acos(x),    |x|<= 1}
#'   \deqn{Tn(x) = cosh(n \cdot acosh(x),  |x|> 1}
#' }
#' \if{html}{\preformatted{
#'   Tn(x) = cos(n . acos(x),    |x|<= 1
#'   Tn(x) = cosh(n . acosh(x),  |x|> 1
#' }}
#' If \code{x} is a vector, the output is a vector of the same size, where each
#' element is calculated as \eqn{y(i) = Tn(x(i))}.
#'
#' @param n Order of the polynomial, specified as a positive integer.
#' @param x Point or points at which to calculate the Chebyshev polynomial
#'
#' @return Polynomial of order \code{x}, evaluated at point(s) \code{x}.
#'
#' @examples
#'
#' cp <- cheb(5, 1)
#' cp <- cheb(5, c(2,3))
#'
#' @author André Carezia, \email{acarezia@@uol.com.br}.\cr Conversion to R by
#'   Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#
#' @export

cheb <- function(n, x) {

  if (!isPosscal(n) || !isWhole(n) || n <= 0)
    stop("n must be an integer strictly positive")

  T <- rep(0, length(x))

  idx <- abs(x) <= 1
  if (any(idx)) {
    T[idx] <- cos(n * acos(signif(as.complex(x[idx], 15))))
  }

  idx <- abs(x) > 1
  if (any(idx)) {
    T[idx] <- cosh(n * acosh(signif(as.complex(x[idx]))))
  }

  Re(T)
}
