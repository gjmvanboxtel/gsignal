# isstable.R
# Copyright (C) 2026 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Octave signal package:
# Copyright (C) 2017 Vasilis Lefkopoulos <vlefkopo@gmail.com>
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
# 20261008 Geert van Boxtel          First version for v0.4.0
#------------------------------------------------------------------------------

#' Determine whether a digital filter is stable
#'
#' A stable digital filter is one where the impulse response approaches zero 
#' as time goes to infinity, ensuring that the output remains bounded for any
#' bounded input. Stability is typically guaranteed if all poles of the
#' filter's transfer function are located inside the unit circle in the
#' complex plane.
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
#' @return Logical indicating whether the filter is stable
#'
#' @examples
#' 
#' b <- c(1, 2, 3, 4, 5, 5, 1, 2)
#' a <- c(4, 5, 6, 7, 9, 10, 4, 6)
#' zplane(b, a)
#' isstable (b, a)
#' 
#' zp <- butter(6, 0.7, 'high', output = 'Zpg')
#' zplane(zp)
#' isstable (zp)
#'  
#' @author Vasilis Lefkopoulos \email{vlefkopo@@gmail.com},\cr
#'  conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#'
#' 
#' @rdname isstable
#' @export
#' 
isstable <- function(filt, ...) UseMethod("isstable")

#' @rdname isstable
#' @method isstable Arma
#' @export
isstable.Arma <- function(filt, ...) # IIR
  isstable(filt$b, filt$a, ...)

#' @rdname isstable
#' @method isstable Ma
#' @export
isstable.Ma <- function(filt, ...) # FIR filter is always stable
  TRUE

#' @rdname isstable
#' @method isstable Sos
#' @export
isstable.Sos <- function(filt, ...) { # Second-order sections
  isstable(as.Zpg(filt), ...)
}

#' @rdname isstable
#' @method isstable Zpg
#' @export
isstable.Zpg <- function(filt, ...) { # zero-pole-gain form
  
  if (any(abs(filt$p) >= 1)) {
    flag <- FALSE
  } else {
    flag <- TRUE
  }
  flag
}

#' @rdname isstable
#' @method isstable default
#' @export

isstable.default <- function(filt, a = 1, ...) {
  
  if (missing(filt)) {
    stop("Invalid call to isstable()")
  }
  
  b <- filt
  if (!(is.numeric(b) || is.complex(b)) ||
      !(is.numeric(a) || is.complex(a))) {
    stop("input must be numeric or complex")
  }

  p <- abs(pracma::roots(a))
  if (any(abs(p) >= 1)) {
    flag <- FALSE
  } else {
    flag <- TRUE
  }
  flag
}
