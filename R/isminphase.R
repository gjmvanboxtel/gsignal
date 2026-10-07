# isminphase.R
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
# 20261007 Geert van Boxtel          First version for v0.4.0
#------------------------------------------------------------------------------

#' Determine whether a digital filter is minimum phase
#'
#' A minimum-phase filter is a type of causal filter that has all its poles and
#' zeros inside the unit circle in the z-plane. Note that minimum-phase filters
#' are stable by definition since the poles must be inside the unit circle. In
#' addition, because the zeros must also be inside the unit circle, the inverse
#' filter is also stable when the filter is minimum phase.
#'
#' @param filt For the default case, the moving-average coefficients of an ARMA
#'   filter (normally called \code{b}), specified as a numeric or complex
#'   vector. Generically, \code{filt} specifies an arbitrary model or filter
#'   operation.
#' @param a the autoregressive (recursive) coefficients of an ARMA filter,
#'   specified as a numeric or complex vector. If \code{a[1]} is not equal to 1,
#'   then filter normalizes the filter coefficients by \code{a[1]}. Therefore,
#'   \code{a[1]} must be nonzero.
#' @param tol tolerance to determine when two numbers are close enough to be
#'   considered equal. Default: .Machine$longdouble.eps^(3 / 4)
#' @param ... additional arguments (ignored).
#'
#' @return Logical indicating whether the filter is minimum phase
#'
#' @examples
#' 
#' b <- c(3, 1)
#' a <- c(1, .5)
#' zplane(b, a)
#' isminphase(b, a)
#' 
#' zp <- butter(6, 0.25, output="Zpg")
#' zplane(zp)
#' isminphase(zp)
#'   
#' @author Leonardo Araujo \email{leolca@@gmail.com},\cr
#'  conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#'
#' @references Oppenheim, A., and Schafer, R. Discrete-Time Signal Processing.
#' 3rd edition, Pearson, 2009
#' 
#' @seealso \code{\link{ismaxphase}}
#'
#' @rdname isminphase
#' @export
#' 
isminphase <- function(filt, ...) UseMethod("isminphase")

#' @rdname isminphase
#' @method isminphase Arma
#' @export
isminphase.Arma <- function(filt, ...) # IIR
  isminphase(filt$b, filt$a, ...)

#' @rdname isminphase
#' @method isminphase Ma
#' @export
isminphase.Ma <- function(filt, ...) # FIR
  isminphase(unclass(filt), ...)

#' @rdname isminphase
#' @method isminphase Sos
#' @export
isminphase.Sos <- function(filt, ...) { # Second-order sections
  isminphase(as.Arma(filt), ...)
}

#' @rdname isminphase
#' @method isminphase Zpg
#' @export
isminphase.Zpg <- function(filt, ...) # zero-pole-gain form
  isminphase(as.Arma(filt), ...)

#' @rdname isminphase
#' @method isminphase default
#' @export

isminphase.default <- function(filt, a = 1,
                               tol = .Machine$longdouble.eps^(3 / 4), ...) {
  
  if (missing(filt)) {
    stop("Invalid call to isminphase()")
  }
  
  b <- filt
  if (!(is.numeric(b) || is.complex(b)) ||
      !(is.numeric(a) || is.complex(a))) {
    stop("input must be numeric or complex")
  }
  
  if(!isPosscal(tol)) {
    stop("tol must be a positive scalar")
  }
  
  z <- abs(pracma::roots(b))
  p <- abs(pracma::roots(a))
  # Zeros of a maximum phase filter are constrained to lie inside the unit
  # circle. The filter should be stable (poles inside the unit circle).
  (all(z < 1 - tol) || length(z) == 0) && (all(p < 1 - tol) || length(p) == 0)
}
