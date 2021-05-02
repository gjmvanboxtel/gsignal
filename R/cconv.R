# cconv.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Octave version Copyright (C) 2018 Leonardo Araujo
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
# 20200209  GvB       setup for gsignal v0.1.0
# 20200212  GvB       use stats::nextn to ensure that the length is
#                     a highly composite number
#------------------------------------------------------------------------------

#' Circular convolution
#'
#' Compute the modulo-n circular convolution.
#'
#' Linear and circular convolution are fundamentally different operations.
#' Linear convolution of an n-point vector x, and an l-point vector y, has
#' length n + l - 1, and can be computed by the function \code{\link{conv}},
#' which uses \code{\link{filter}}. The circular convolution, by contrast, is
#' equal to the inverse discrete Fourier transform (DFT) of the product of the
#' vectors' DFTs.
#'
#' For the circular convolution of \code{x} and \code{y} to be equivalent to
#' their linear convolution, the vectors must be padded with zeros to length at
#' least \code{n + l - 1} before taking the DFT. After inverting the product of
#' the DFTs, only the first \code{n + l - 1} elements should be retained.
#'
#' For long sequences circular convolution may be more efficient than linear
#' convolution. You can also use \code{cconv} to compute the circular
#' cross-correlation of two sequences.
#'
#' @param a,b Input, coerced to vectors, can be different lengths or data types.
#' @param n Convolution length, specified as a positive integer. Default:
#'   \code{length(a) + length(b) - 1}.
#'
#' @return Circular convolution of input vectors, returned as a vector.
#'
#' @examples
#' a <- c(1, 2, -1, 1)
#' b <- c(1, 1, 2, 1, 2, 2, 1, 1)
#' c <- cconv(a, b)       # Circular convolution
#' cref = conv(a, b)      # Linear convolution
#' all.equal(max(c - cref), 0)
#'
#' cconv(a, b, 6)
#'
#' @seealso \code{\link{conv}}, \code{\link[stats]{convolve}}
#'
#' @author Leonardo Araujo.\cr Conversion to R by Geert van Boxtel,
#'   \email{G.J.M.vanBoxtel@@gmail.com}.
#
#' @export

cconv <- function(a, b, n = length(a) + length(b) - 1) {

  a <- as.vector(a)
  b <- as.vector(b)
  if (!isPosscal(n) || !isWhole(n)) {
    stop("n must be a positive integer")
  }

  la <- length(a)
  lb <- length(b)

  if (la < lb) {
    a <- postpad(a, lb)
  } else if (lb < la) {
    b <- postpad(b, la)
  }

  N <- length(a)
  if (n < N) {
    an <- bn <- rep(0L, n)
    for (i in 0:(N - 1)) {
      modi <- i %% n
      an[modi + 1] <- an[modi + 1] + a[i + 1]
      bn[modi + 1] <- bn[modi + 1] + b[i + 1]
    }
    a <- an
    b <- bn
  } else if (n > N) {
    a <- postpad(a, n)
    b <- postpad(b, n)
  }

  y <- ifft(stats::fft(postpad(a, stats::nextn(length(a)))) *
              stats::fft(postpad(b, stats::nextn(length(b)))))
  y[1:n]
}
