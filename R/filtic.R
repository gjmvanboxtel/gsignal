# filtic.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Octave version:
# Copyright (C) 2004 David Billinghurst <David.Billinghurst@riotinto.com>
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
# 20200317    GvB       setup for gsignal v0.1.0
# 20210322    GvB       adapted to accept missing x and y parameters (all 1's)
#                       defined S3 methods and added method for Sos
#------------------------------------------------------------------------------

#' Filter Initial Conditions
#'
#' Compute the initial conditions for a filter.
#'
#' This function computes the same values that would be obtained from the
#' function \code{filter} given past inputs \code{x} and outputs \code{y}.
#'
#' The vectors \code{x} and \code{y} contain the most recent inputs and outputs
#' respectively, with the newest values first:
#'
#' \code{x = c(x(-1), x(-2), ... x(-nb)); nb = length(b)-1}\cr
#' \code{y = c(y(-1), y(-2), ... y(-na)); na = length(a)-a}
#'
#' If \code{length(x) < nb} then it is zero padded. If \code{length(y) < na}
#' then it is zero padded.
#'
#' @param filt For the default case, the moving-average coefficients of an ARMA
#'   filter (normally called ‘b’), specified as a vector. Generically,
#'   \code{filt} specifies an arbitrary filter operation.
#' @param a the autoregressive (recursive) coefficients of an ARMA filter.
#' @param y output vector, with the most recent values first.
#' @param x input vector, with the most recent values first. Default: 0
#' @param ... additional arguments (ignored).
#'
#' @return Initial conditions for filter specified by \code{filt}, input vector
#'   \code{x}, and output vector \code{y}, returned as a vector.
#'
#' @examples
#' ## Simple low pass filter
#' b <- c(0.25, 0.25)
#' a <- c(1.0, -0.5)
#' ic <- filtic(b, a, 1, 1)
#'
#' ## Simple high pass filter
#' b <- c(0.25, -0.25)
#' a <- c(1.0, 0.5)
#' ic <- filtic(b, a, 0, 1)
#'
#' ## Example from Python scipy.signal.lfilter() documentation
#' t <- seq(-1, 1, length.out =  201)
#' x <- (sin(2 * pi * 0.75 * t * (1 - t) + 2.1)
#'       + 0.1 * sin(2 * pi * 1.25 * t + 1)
#'       + 0.18 * cos(2 * pi * 3.85 * t))
#' h <- butter(3, 0.05)
#' l <- max(length(h$b), length(h$a)) - 1
#' zi <- filtic(h, rep(1, l), rep(1, l))
#' z <- filter(h, x, zi * x[1])
#'
#' @seealso \code{\link{filter}}, \code{\link{sosfilt}}, \code{\link{filtfilt}},
#' \code{\link{filter_zi}}
#'
#' @author David Billinghurst, \email{David.Billinghurst@@riotinto.com}.\cr
#'   Adapted and converted to R by Geert van Boxtel
#'   \email{G.J.M.vanBoxtel@@gmail.com}.
#'
#' @rdname filtic
#' @export

filtic <- function(filt, ...) UseMethod("filtic")

#' @rdname filtic
#' @method filtic default
#' @export

filtic.default <- function(filt, a, y, x = 0, ...) {

  b <- filt
  nz <- max(length(a), length(b)) - 1
  zi <- numeric(nz)

  # Pad arrays a and b to length nz+1 if required
  if (length(a) < (nz + 1)) {
    a <- postpad(a, nz + 1)
  }
  if (length(b) < (nz + 1)) {
    b <- postpad(b, nz + 1)
  }

  # Pad arrays x and y to length nz if required
  if (length(x) < nz) {
    x <- postpad(x, nz)
  }
  if (length(y) < nz) {
    y <- postpad(y, nz)
  }

  for (i in seq(nz, 1, -1)) {
    for (j in i:(nz - 1)) {
      zi[j] <- b[j + 1] * x[i] - a[j + 1] * y[i] + zi[j + 1]
    }
    zi[nz] <- b[nz + 1] * x[i] - a[nz + 1] * y[i]
  }

  zi <- zi / a[1]
  zi
}

#' @rdname filtic
#' @method filtic Arma
#' @export
filtic.Arma <- function(filt, y, x = 0, ...) # IIR
  filtic(filt$b, filt$a, y, x, ...)

#' @rdname filtic
#' @method filtic Ma
#' @export
filtic.Ma <- function(filt, y, x = 0, ...) # FIR
  filtic(unclass(filt), 1, y, x, ...)

#' @rdname filtic
#' @method filtic Sos
#' @export
filtic.Sos <- function(filt, y, x = 0, ...) { # Second-order sections

  if (filt$g != 1) {
    filt$sos[1, 1:3] <- filt$sos[1, 1:3] * filt$g
  }
  L <- NROW(filt$sos)
  zi <- matrix(0, L, 2)
  scale <- 1.0
  for (l in seq_len(L)) {
    b <- filt$sos[l, 1:3]
    a <- filt$sos[l, 4:6]
    zi[l, ] <- scale * filtic.default(b, a, y, x)
    # If H(z) = B(z)/A(z) is this section's transfer function, then
    # b.sum()/a.sum() is H(1), the gain at omega=0.  That's the steady
    # state value of this section's step response.
    scale <- scale *  sum(b) / sum(a)
  }
  zi
}

#' @rdname filtic
#' @method filtic Zpg
#' @export
filtic.Zpg <- function(filt, y, x = 0, ...) # zero-pole-gain form
  filtic(as.Arma(filt), y, x, ...)
