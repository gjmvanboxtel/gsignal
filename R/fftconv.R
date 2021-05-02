# fftconv.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Octave function:
# Copyright (C) 1994-2017 John W. Eaton
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
# 20200420  GvB       setup for gsignal v0.1.0
#------------------------------------------------------------------------------

#' FFT-based convolution
#'
#' Convolve two vectors using the FFT for computation.
#'
#' The computation uses the FFT by calling the function \code{fftfilt}. If the
#' optional argument \code{n} is specified, an \code{n}-point overlap-add FFT is
#' used.
#'
#' @param x,y input vectors.
#' @param n FFT length, specified as a positive integer. The FFT size must be an
#'   even power of 2 and must be greater than or equal to the length of
#'   \code{filt}. If the specified \code{n} does not meet these criteria, it is
#'   automatically adjusted to the nearest value that does. If \code{n = NULL}
#'   (default), then the overlap-add method is not used.
#'
#' @return Convoluted signal, specified as a a vector of length equal to
#'   \code{length (x) + length (y) - 1}. If \code{x} and \code{y} are the
#'   coefficient vectors of two polynomials, the returned value is the
#'   coefficient vector of the product polynomial.
#'
#' @examples
#'
#' u <- rep(1L, 3)
#' v <- c(1, 1, 0, 0, 0, 1, 1)
#' w1 <- conv(u, v)              # time-domain convolution
#' w2 <- fftconv(u, v)           # frequency domain convolution
#' all.equal(w1, w2)             # same results
#'
#' @seealso \code{\link{conv}}, \code{\link{conv2}}
#'
#' @author Kurt Hornik, \email{Kurt.Hornik@@wu-wien.ac.at},\cr
#'  adapted by John W. Eaton.\cr
#'  Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#'
#' @export

fftconv <- function(x, y, n = NULL) {

  if (!(is.vector(x) && is.vector(y))) {
    stop("both x and y must be vectors")
  }

  lx <- length(x)
  ly <- length(y)

  if ((lx == 1) || (ly == 1)) {
    z <- x * y
  } else {
    lz <- lx + ly - 1
    x <- postpad(x, lz)
    y <- postpad(y, lz)
    z <- fftfilt(x, y, n)
  }
  z
}
