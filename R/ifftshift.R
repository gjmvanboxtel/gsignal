# ifftshift.R
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

#' Inverse zero-frequency shift
#'
#' Rearranges a zero-frequency-shifted Fourier transform back to the original.
#'
#' Undo the action of the fftshift function. For even length \code{x},
#' \code{fftshift} is its own inverse, but not for odd length input.
#'
#' @param x input data, specified as a vector or matrix.
#' @param MARGIN dimension to operate along, 1 = row, 2 = columns (default).
#'   Specifying \code{MARGIN = c(1, 2)} centers along both rows and columns.
#'   Ignored when \code{x} is a vector.
#'
#' @return back-transformed vector or matrix.
#'
#' @examples
#' Xeven <- 1:6
#' res <- fftshift(fftshift(Xeven))
#'
#' Xodd <- 1:7
#' res <- fftshift(fftshift(Xodd))
#' res <- ifftshift(fftshift(Xodd))
#'
#' @seealso \code{\link{fftshift}}
#'
#' @author Vincent Cautaerts, \email{vincent@@comf5.comm.eng.osaka-u.ac.jp},\cr
#'   adapted by John W. Eaton.\cr
#'   Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#
#' @export

ifftshift <- function(x, MARGIN = 2) {

  y <- x
  if (is.vector(y)) {
    xl <- length(y)
    if (xl > 1) {
      xx <- floor(xl / 2)
      y <- x[c((xx + 1):xl, 1:xx)]
    }
  } else if (is.matrix(y)) {
    if (! (1 %in% MARGIN ||  2 %in% MARGIN)) {
      stop("MARGIN must be 1, 2, or both")
    }
    if (1 %in% MARGIN) {
      nr <- NROW(y)
      if (nr > 1) {
        xx <- floor(nr / 2)
        y <- y[c((xx + 1):nr, 1:xx), ]
      }
    }
    if (2 %in% MARGIN) {
      nc <- NCOL(y)
      if (nc > 1) {
        xx <- floor(nc / 2)
        y <- y[, c((xx + 1):nc, 1:xx)]
      }
    }
  } else {
    stop("x must be a vector or matrix")
  }
  y
}
