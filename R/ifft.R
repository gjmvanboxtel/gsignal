# ifft.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
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
# 20200113  GvB       setup for gsignal v0.1.0
# 20200123  GvB       adapted rest for Imaginary component == 0 using zapsmall
# 20211031  GvB       replaced zapsmall test by test for conjugate symmetry
#------------------------------------------------------------------------------

#' Inverse Fast Fourier Transform
#'
#' Compute the inverse Fast Fourier Transform compatible with
#' 'Matlab' and 'Octave'.
#'
#' The \code{'fft'} function in the \code{'stats'} package can compute the
#' inverse FFT by specifying \code{inverse = TRUE}. However, that function does
#' \emph{not} divide the result by \code{length(x)}, nor does it return real
#' values when appropriate. The present function does both, and is this
#' compatible with 'Matlab' and 'Octave' (and differs from the \code{'ifft'}
#' function in the \code{'signal'} package, which does not return real values).
#'
#' @param x Real or complex vector, array, or matrix.
#'
#' @return When \code{x} is a vector, the value computed and returned by
#'   \code{ifft} is the univariate inverse discrete Fourier transform of the
#'   sequence of values in \code{x}. Specifically, \code{y <- ifft(x)} is
#'   defined as \code{stats::fft(x, inverse = TRUE) / length(x)}. The
#'   \code{stats::fft} function called with \code{inverse = TRUE} replaces
#'   \code{exp(-2 * pi...)} with \code{exp(2 * pi)} in the definition of the
#'   discrete Fourier transform (see \code{\link[stats]{fft}}).
#'
#' When \code{x} contains an array, \code{ifft} computes and returns the
#' normalized inverse multivariate (spatial) transform. By contrast,
#' \code{imvfft} takes a real or complex matrix as argument, and returns a
#' similar shaped matrix, but with each column replaced by its normalized
#' inverse discrete Fourier transform. This is useful for analyzing
#' vector-valued series.
#'
#' @examples
#' res <- ifft(stats::fft(1:5))
#' res <- ifft(stats::fft(c(1+5i, 2+3i, 3+2i, 4+6i, 5+2i)))
#' res <- imvfft(stats::mvfft(matrix(1:20, 4, 5)))
#'
#' @seealso \code{\link[stats]{fft}}
#'
#' @author Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#
#' @export

ifft <- function(x) {

  y <- stats::fft(x, inverse = TRUE) / length(x)
  if (isConjSymm(x)) Re(y) else y
}

#' @rdname ifft
#' @export

imvfft <- function(x) {

  y <- stats::mvfft(x, inverse = TRUE) / nrow(x)
  if(all(apply(x, 2, isConjSymm))) Re(y) else y
}
