# ismaxphase.R
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

#' Determine whether a digital filter is maximum phase
#'
#' A maximum-phase filter is a type of causal filter that has all its zeros
#' outside the unit circle in the z-plane, so its numerator coefficients b are
#' typically larger in magnitude than the denominator coefficients a, resulting
#' in the largest group delay among all causal filters with the same magnitude
#' response.
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
#' @return Logical indicating whether the filter is maximum phase
#'
#' @examples
#' 
#' # maximum-phase FIR filter
#' fir <- Ma(c(0.5, -1))
#' zplane(fir)
#' ismaxphase(fir)
#' 
#' #  maximum-phase IIR filter
#' iir <- Arma(conv(c(1, -1.2), c(1, -1.5)), conv(c(1, -0.6), c(1, -0.7))) 
#' zplane(iir)
#' ismaxphase(iir)
#' 
#' sos <- Sos(c(1, -2.7, 1.8, 1, -1.3, 0.42))
#' zplane(sos)
#' ismaxphase(sos)
#'  
#' @author Leonardo Araujo \email{leolca@@gmail.com},\cr
#'  conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#'
#' @references Oppenheim, A., and Schafer, R. Discrete-Time Signal Processing.
#' 3rd edition, Pearson, 2009
#' 
#' @seealso \code{\link{isminphase}}
#'
#' @rdname ismaxphase
#' @export
#' 
ismaxphase <- function(filt, ...) UseMethod("ismaxphase")

#' @rdname ismaxphase
#' @method ismaxphase Arma
#' @export
ismaxphase.Arma <- function(filt, ...) # IIR
  ismaxphase(filt$b, filt$a, ...)

#' @rdname ismaxphase
#' @method ismaxphase Ma
#' @export
ismaxphase.Ma <- function(filt, ...) # FIR
  ismaxphase(unclass(filt), ...)

#' @rdname ismaxphase
#' @method ismaxphase Sos
#' @export
ismaxphase.Sos <- function(filt, ...) { # Second-order sections
  ismaxphase(as.Arma(filt), ...)
}

#' @rdname ismaxphase
#' @method ismaxphase Zpg
#' @export
ismaxphase.Zpg <- function(filt, ...) # zero-pole-gain form
  ismaxphase(as.Arma(filt), ...)

#' @rdname ismaxphase
#' @method ismaxphase default
#' @export

ismaxphase.default <- function(filt, a = 1,
                               tol = .Machine$longdouble.eps^(3 / 4), ...) {
  
  if (missing(filt)) {
    stop("Invalid call to ismaxphase()")
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
  # Zeros of a maximum phase filter are constrained to lie outside the unit
  # circle. The filter should be stable (poles inside the unit circle).
  (all(z > 1 + tol) || length(z) == 0) && (all(p < 1 - tol) || length(p) == 0)
}
