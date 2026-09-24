# filtord.R
# Copyright (C) 2026 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Octave signal package:
# Copyright (C) 2023 Leonardo Araujo <leolca@gmail.com>
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
# 20260924 Geert van Boxtel          First version for v0.4.0
#------------------------------------------------------------------------------

#' Filter order
#'
#' Returns the filter order for the specified filter.
#' 
#' The filter order refers to the number of coefficients used in a digital
#' filter, which affects its performance and characteristics. For example, a
#' higher filter order typically results in a steeper roll-off and better
#' frequency response.
#' 
#' The function finds the highest power among numerator polynomial terms with
#' nonzero coefficients, finds the highest power for the denominator, and
#' returns the larger of the two. For FIR filters, the denominator equals 1, so
#' effectively the function returns the highest power among transfer function
#' polynomial terms with nonzero coefficients
#'
#' @param filt For the default case, the moving-average coefficients of an ARMA
#'   filter (normally called \code{b}), specified as a numeric or complex
#'   vector. Generically, \code{filt} specifies an arbitrary model or filter
#'   operation.
#' @param a the autoregressive (recursive) coefficients of an ARMA filter,
#'   specified as a numeric or complex vector. If \code{a[1]} is not equal to 1,
#'   then filter normalizes the filter coefficients by \code{a[1]}. Therefore,
#'   \code{a[1]} must be nonzero.
#' @param ... additional arguments (ignored).
#'
#' @return Filter order, specified as a real positive scalar
#'
#' @examples
#' b <- c(1, 0)
#' a <- c(1, 1)
#' n <- filtord(b, a)
#' 
#' sos <- tf2sos(c(1, 0, 0, 0, 0, 0, 0, 1), c(1, 0, 0, 0, 0, 0, 0, .5))
#' n <- filtord(sos)
#' 
#' @author Leonardo Araujo \email{leolca@@gmail.com},\cr
#'  adapted by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#'
#' @seealso \code{\link{butter}}, \code{\link{FilterSpecs}}
#'
#' @rdname filtord
#' @export
#' 
filtord <- function(filt, ...) UseMethod("filtord")

#' @rdname filtord
#' @method filtord Arma
#' @export
filtord.Arma <- function(filt, ...) # IIR
  filtord(filt$b, filt$a, ...)

#' @rdname filtord
#' @method filtord Ma
#' @export
filtord.Ma <- function(filt, ...) # FIR
  filtord(unclass(filt), ...)

#' @rdname filtord
#' @method filtord Sos
#' @export
filtord.Sos <- function(filt, ...) { # Second-order sections
  filtord(as.Arma(filt), ...)
}

#' @rdname filtord
#' @method filtord Zpg
#' @export
filtord.Zpg <- function(filt, ...) # zero-pole-gain form
  filtord(as.Arma(filt), ...)

#' @rdname filtord
#' @method filtord default
#' @export

filtord.default <- function(filt, a = 1, ...) {
  max(length(filt) - 1, length(a) - 1)
}
