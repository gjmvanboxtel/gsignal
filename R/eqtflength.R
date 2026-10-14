# eqtflength.R
# Copyright (C) 2026 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Original Octave version:
# Copyright (C) 2026 John Donoghue
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
# 20261014  GvB       setup for gsignal v0.4.0
#------------------------------------------------------------------------------

#' Equalize numerator and denominator polynomial lengths
#'
#' Pads the shorter polynomial with trailing zeros so that the returned
#' numerator b and denominator a have equal length while representing the same
#' discrete-time transfer function.
#'
#' @param b moving average (MA) polynomial coefficients
#' @param a autoregressive (AR) polynomial coefficients
#'
#' @return A list with the following list elements:
#' \describe{
#'   \item{b}{Numerator polynomial coefficients of transfer function. \code{b}
#'   has the same length as \code{a}.}
#'  \item{a}{Denominator polynomial coefficients of transfer function. \code{a}
#'   has the same length as \code{b}.}
#'  \item{n}{Numerator order. Any trailing zeros in \code{b} are excluded when
#'  computing \code{n}.}
#'  \item{m}{Denominator order. Any trailing zeros in \code{a} are excluded when
#'  computing \code{m}.}
#' }
#'
#' @examples
#' eqtflength (c(1, 2), c(1, 0.5, 0.25))
#'
#' # example with complex coefficients:
#' b <- c(0.2 + 0.1i, 0.3 - 0.2i, 0.1 + 0.05i)
#' a <- c(1, -0.9 - 0.1i)
#' eqtflength(b, a)
#'
#' @author Copyright (C) 2026 John Donoghue.\cr
#' Conversion to R by Geert van Boxtel, \email{gjmvanboxtel@@gmail.com}.
#'
#' @export

eqtflength <- function(b, a) {
  
  if(!(is.numeric(b) || is.complex(b)) ||
     !(is.numeric(a) || is.complex(a))) {
    stop("b and a must be numeric or complex")
  }

  if (length(b) == 0) {
    n <- -1
  }
  idx <- utils::tail(which(b != 0), 1)
  if (length(idx) == 0) {
    ## Zero polynomial
    n <- 0
  } else {
    n <- idx - 1
  }
  
  if (length(a) == 0) {
    m <- -1
  }
  idx <- utils::tail(which(a != 0), 1)
  if (length(idx) == 0) {
    ## Zero polynomial
    m <- 0
  } else {
    m <- idx - 1
  }
  
  ## Equalize lengths
  len <- max(length(b), length(a))
  
  if (length(b) < len) {
    b <- c(b, rep(0, len - length(b)))
  }
  
  if (length(a) < len) {
    a <- c(a, rep(0, len - length(a)))
  }
  
  list(b = b, a = a, n = n, m = m)
}
