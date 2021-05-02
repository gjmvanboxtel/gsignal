# fftshift.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Original Octave version:
# Author: Vincent Cautaerts <vincent@comf5.comm.eng.osaka-u.ac.jp>
# Created: July 1997
# Adapted-By: jwe
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
# 20200823  GvB       setup for gsignal v0.1.0
#------------------------------------------------------------------------------

#' Zero-frequency shift
#'
#' Perform a shift in order to move the frequency 0 to the center of the input.
#'
#' If \code{x} is a vector of \code{N} elements corresponding to \code{N} time
#' samples spaced by \code{dt}, then \code{fftshift(x)} corresponds to
#' frequencies \code{f = c(-seq(ceiling((N-1)/2), 1, -1), 0, (1:floor((N-1)/2)))
#' * df}, where \code{df = 1 / (N * dt)}. In other words, the left and right
#' halves of \code{x} are swapped.
#'
#' If \code{x} is a matrix, then \code{fftshift} operates on the rows or columns
#' of \code{x}, according to the \code{MARGIN} argument, i.e. it swaps the the
#' upper and lower halves of the matrix \code{(MARGIN = 1)}, or the left and
#' right halves of the matrix \code{(MARGIN = 2)}. Specifying \code{MARGIN =
#' c(1, 2)} swaps along both dimensions, i.e., swaps the first quadrant with the
#' fourth, and the second with the third.
#'
#' @param x input data, specified as a vector or matrix.
#' @param MARGIN dimension to operate along, 1 = row, 2 = columns (default).
#'   Specifying \code{MARGIN = c(1, 2)} centers along both rows and columns.
#'   Ignored when \code{x} is a vector.
#'
#' @return vector or matrix with centered frequency.
#'
#' @examples
#' Xeven <- 1:6
#' ev <- fftshift(Xeven)   # returns 4 5 6 1 2 3
#'
#' Xodd <- 1:7
#' odd <- fftshift(Xodd)   # returns 5 6 7 1 2 3 4
#'
#' fs <- 100                      # sampling frequency
#' t <- seq(0, 10 - 1/fs, 1/fs)   # time vector
#' S <- cos(2 * pi * 15 * t)
#' n <- length(S)
#' X <- fft(S)
#' f <- (0:(n - 1)) * (fs / n);   # frequency range
#' power <- abs(X)^2 / n          # power
#' plot(f, power, type="l")
#' Y <- fftshift(X)
#' fsh <- ((-n/2):(n/2-1)) * (fs / n)  # zero-centered frequency range
#' powersh <- abs(Y)^2 / n             # zero-centered power
#' plot(fsh, powersh, type = "l")
#'
#' @seealso \code{ifftshift}
#'
#' @author Vincent Cautaerts, \email{vincent@@comf5.comm.eng.osaka-u.ac.jp},\cr
#'   adapted by John W. Eaton.\cr
#'   Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#
#' @export

fftshift <- function(x, MARGIN = 2) {

  y <- x
  if (is.vector(y)) {
    xl <- length(y)
    if (xl > 1) {
      xx <- ceiling(xl / 2)
      y <- x[c((xx + 1):xl, 1:xx)]
    }
  } else if (is.matrix(y)) {
    if (! (1 %in% MARGIN ||  2 %in% MARGIN)) {
      stop("MARGIN must be 1, 2, or both")
    }
    if (1 %in% MARGIN) {
      nr <- NROW(y)
      if (nr > 1) {
        xx <- ceiling(nr / 2)
        y <- y[c((xx + 1):nr, 1:xx), ]
      }
    }
    if (2 %in% MARGIN) {
      nc <- NCOL(y)
      if (nc > 1) {
        xx <- ceiling(nc / 2)
        y <- y[, c((xx + 1):nc, 1:xx)]
      }
    }
  } else {
    stop("x must be a vector or matrix")
  }
  y
}
