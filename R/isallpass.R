# isallpass.R
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
# 20261005 Geert van Boxtel          First version for v0.4.0
#------------------------------------------------------------------------------

#' Determine whether a digital filter is allpass
#'
#' An all-pass filter is a filter that passes all frequencies equally in gain,
#' but changes the phase relationship among various frequencies.
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
#' @return Logical indicating whether the filter is allpass
#'
#' @examples
#' a <- c(1, 2, 3)
#' b <- c(3, 2, 1)
#' isallpass (b, a)
#' 
#' # H(z) = (b1 - z^-1) * (b2 - z^-1) / ((1 - b1*z^-1) * (1 - b2*z^-1))
#' b1 <- 0.5 * (1 + 1i)
#' b2 <- 0.7 * (cos(pi / 6) + 1i*sin(pi / 6))
#' b <- conv(c(b1, -1), c(b2, -1))
#' a <- conv(c(1, (-1) * Conj(b1)), c(1, (-1) * Conj(b2)))
#' freqz (b, a)
#' isallpass (b, a)
#'  
#' @author Leonardo Araujo \email{leolca@@gmail.com},\cr
#'  conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#'
#' @references [1] Shyu, Jong-Jy, & Pei, Soo-Chang,
#' A new approach to the design of complex all-pass IIR digital filters,
#' Signal Processing, 40(2–3), 207–215, 1994.
#' https://doi.org/10.1016/0165-1684(94)90068-x
#'
#' @references [2] Vaidyanathan, P. P. Multirate Systems and Filter Banks.
#' 1st edition, Pearson College Div, 1992.
#'
#' @rdname isallpass
#' @export
#' 
isallpass <- function(filt, ...) UseMethod("isallpass")

#' @rdname isallpass
#' @method isallpass Arma
#' @export
isallpass.Arma <- function(filt, ...) # IIR
  isallpass(filt$b, filt$a, ...)

#' @rdname isallpass
#' @method isallpass Ma
#' @export
isallpass.Ma <- function(filt, ...) # FIR
  isallpass(unclass(filt), ...)

#' @rdname isallpass
#' @method isallpass Sos
#' @export
isallpass.Sos <- function(filt, ...) { # Second-order sections
  isallpass(as.Arma(filt), ...)
}

#' @rdname isallpass
#' @method isallpass Zpg
#' @export
isallpass.Zpg <- function(filt, ...) # zero-pole-gain form
  isallpass(as.Arma(filt), ...)

#' @rdname isallpass
#' @method isallpass default
#' @export

isallpass.default <- function(filt, a = 1, ...) {
  
  if (missing(filt)) {
    stop("Invalid call to isallpass()")
  }
  
  b <- filt
  if (!(is.numeric(b) || is.complex(b)) ||
      !(is.numeric(a) || is.complex(a))) {
    stop("input must be numeric or complex")
  }
  
  if (length(b) != length(a)) {
    flag <- FALSE
  } else {
  
    # remove leading and trailing zeros
    b <- b[seq(which(b != 0)[1], which(b != 0)[length(which(b != 0))])]
    a <- a[seq(which(a != 0)[1], which(a != 0)[length(which(a != 0))])]
  
    # normalize
    b <- b / b[length(b)]
    a <- a / a[1]

    flag <- all(b == rev(Conj(a))) || all(b == -rev(Conj(a)))
  }
  
  flag
}
