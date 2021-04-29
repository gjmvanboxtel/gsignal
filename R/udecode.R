# udecode.R
# Copyright (C) 2019 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Octave signal package:
# Copyright (C) 2014 Georgios Ouzounis <ouzounis_georgios@hotmail.com>
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
# 20191208 Geert van Boxtel          First version for v0.1.0
#------------------------------------------------------------------------------

#' Uniform decoder
#'
#' Decode \eqn{2^n}-level quantized integer inputs to floating-point outputs.
#'
#' \code{y <- udecode(u, n)} inverts the operation of \code{uencode} and
#' reconstructs quantized floating-point values from an encoded multidimensional
#' array of integers \code{u}. The input argument \code{n} must be an integer
#' between 2 and 32. The integer \code{n} specifies that there are \eqn{2^{n}}
#' quantization levels for the inputs, so that entries in \code{u} must be
#' either:
#' \itemize{
#'   \item Signed integers in the range \eqn{-2^{n}/2} to \eqn{(2^{n}/2) - 1}
#'   \item Unsigned integers in the range 0 to \eqn{2^{n} - 1}
#' }
#'
#' Inputs can be real or complex values of any integer data type. Overflows
#' (entries in u outside of the ranges specified above) are saturated to the
#' endpoints of the range interval. The output has the same dimensions as the
#' input \code{u}. Its entries have values in the range -1 to 1.
#'
#' \code{y <- udecode(u, n, v)} decodes \code{u} such that the output has values
#' in the range \code{-v} to \code{v}, where the default value for \code{v} is
#' 1.
#'
#' \code{y <- udecode(u, n, v, saturate)} decodes \code{u} and treats input
#' overflows (entries in \code{u} outside of the range \code{-v} to \code{v}
#' according to \code{saturate}, which can be set to one of the following:
#' \itemize{
#'   \item TRUE (default). Saturate overflows.
#'     \itemize{
#'       \item Entries in signed inputs \code{u} whose values are outside of the
#'       range \eqn{-2^{n}/2} to \eqn{(2^{n}/2) – 1} are assigned the value
#'       determined by the closest endpoint of this interval.
#'        \item Entries in unsigned inputs \code{u} whose values are outside of
#'        the range 0 to \eqn{2^{n}-1} are assigned the value determined by the
#'        closest endpoint of this interval.
#'     }
#'   \item FALSE Wrap all overflows according to the following:
#'     \itemize{
#'       \item Entries in signed inputs \code{u} whose values are outside of the
#'       range \eqn{-2^{n}/2} to \eqn{(2^{n}/2) – 1} are wrapped back into that
#'       range using modulo \eqn{2^{n}} arithmetic (calculated using \eqn{u =
#'       mod(u+2^{n}/2, 2^{n})-(2^{n}/2))}.
#'       \item Entries in unsigned inputs \code{u} whose values are outside of
#'       the range 0 to \eqn{2^{n}-1} are wrapped back into the required range
#'       before decoding using modulo \eqn{2^{n}} arithmetic (calculated using
#'       \eqn{u = mod(u,2^{n}))}.
#'     }
#' }
#'
#' @param u Input, a multidimensional array of integer numbers (can be complex).
#' @param n Number of levels used in \eqn{2^{n}}-level quantization. \code{n}
#'   must be between 2 and 32
#' @param v Limit on the range of \code{u} to the range from \code{-v} to
#'   \code{v} before saturating them. Default 1.
#' @param saturate Logical indicating to saturate (TRUE, default) or to wrap
#'   (FALSE) overflows. See Details.
#'
#' @return Multidimensional array of the same size as \code{u} containing
#'   floating point numbers.
#'
#' @note The real and imaginary components of complex inputs are decoded
#'   independently.
#'
#' @examples
#'
#' u <- c(-1, 1, 2, -5)
#' ysat <- udecode(u, 3)
#'
#' # Notice the last entry in u saturates to 1, the default peak input
#' # magnitude. Change the peak input magnitude to 6.
#' ysatv <- udecode(u, 3, 6)
#'
#' # The last input entry still saturates. Wrap the overflows.
#' ywrap = udecode(u, 3, 6, FALSE)
#'
#' # Add more quantization levels.
#' yprec <- udecode(u, 5)
#'
#' @author Georgios Ouzounis, \email{ouzounis_georgios@@hotmail.com}.\cr
#' Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#'
#' @export

udecode <- function(u, n, v = 1, saturate = TRUE) {

  if (!isScalar(n) || n < 2 || n > 32 || !isWhole(n))
    stop("n must be an integer in the range 2 to 32")
  if (!isPosscal(v) || !isWhole(v) || v <= 0)
    stop("v must be a positive integer")
  if (!is.logical(saturate))
    stop("signed must be a logical")

  # function to do the actual decoding.
  # needed because it is run twice for complex input (real + imaginary)
  decode_it <- function(x) {

    if (all(x >= 0)) {
      signed <- FALSE
      lowerlevel <- 0
      upperlevel <- (2^n) - 1
    } else {
      signed <- TRUE
      lowerlevel <- - 2 ^ (n - 1)
      upperlevel <- (2 ^ (n - 1)) - 1
    }

    if (saturate) {
      if (signed) {
        x[x < lowerlevel] <- lowerlevel
        x[x > upperlevel] <- upperlevel
      } else {
        x[x > upperlevel] <- upperlevel
      }
    } else {
      if (signed) {
        idx <- which(x < lowerlevel | x > upperlevel)
        x[idx] <- ((x[idx] + 2 ^ (n - 1)) %% (2^n)) - 2 ^ (n - 1)
      } else {
        idx <- which(x > upperlevel)
        x[idx] <- x[idx] %% 2^n
      }
    }

    width <- 2 * v / 2^n
    y <- x * width
    if (!signed) y <- y - v

    y
  }

  if (is.complex(u)) {
    real_part <- decode_it(Re(u))
    imag_part <- decode_it(Im(u))
    y <- complex(real = real_part, imaginary = imag_part)
    # complex() returns a vector not a matrix
    if (is.matrix(u)) {
      y <- matrix(y, nrow(u), ncol(u), byrow = TRUE)
    }
  } else {
    y <- decode_it(u)
  }

  y

}
