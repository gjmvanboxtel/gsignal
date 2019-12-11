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
# 20191211 Geert van Boxtel          First version for v0.1.0
#---------------------------------------------------------------------------------------------------------------------------------

#' Chebyshev polynomials
#' 
#' Returns the value of the \code{n}th-order Chebyshev polynomial calculated at the point \code{x}.
#' 
#' The Chebyshev polynomials are defined by the equations:
#' \deqn{Tn(x) = cos(n \cdot acos(x),    |x|<= 1}
#' \deqn{Tn(x) = cosh(n \cdot acosh(x),  |x|> 1}
#' 
#' If \code{x} is a vector, the output is a vector of the same size, where each element is calculated 
#' as \eqn{y(i) = Tn(x(i))}.
#' 
#' @param n Order of the polynomial, specified as a positive integer.
#' @param x Point or points at which to calculate the Chebyshev polynomial
#' 
#' @return polynomial of order \code{x}, evaluated at point(s) \code{x}.
#' 
#' @examples
#' 
#' cp <- cheb(5, 1)
#' cp <- cheb(5, c(2,3))
#'
#' @author Original Octave code Copyright (C) 2002 André Carezia \email{acarezia@@uol.com.br}.
#' Port to R by Geert van Boxtel \email{G.J.M.vanBoxtel@@gmail.com}.
#
#' @export

cheb <- function (n, x) {
  
  if (!isPosscal(n) || ! isWhole(n) || n <= 0) stop ("n must be an integer strictly positive")

  T <- rep(0, length(x))
  
  idx <- abs(x) <= 1
  if (any(idx)) {
    T[idx] <- cos(n * acos(as.complex(x[idx])))
  }
  
  idx <- abs(x) > 1
  if (any(idx)) {
    T[idx] <- cosh(n * acosh(as.complex(x[idx])))
  }
  
  Re(T)
}
