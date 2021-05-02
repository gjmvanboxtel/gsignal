# uencode.R
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
# 20191207 Geert van Boxtel          First version for v0.1.0
#------------------------------------------------------------------------------

#' Uniform encoder
#'
#' Quantize and encode floating-point inputs to integer outputs.
#'
#' \code{y <- uencode(u, n)} quantizes the entries in a multidimensional array
#' of floating-point numbers \code{u} and encodes them as integers using
#' \eqn{2^{n}}-level quantization. \code{n} must be an integer between 2 and 32
#' (inclusive). Inputs can be real or complex, double- or single-precision. The
#' output \code{y} and the input \code{u} are arrays of the same size. The
#' elements of the output \code{y} are unsigned integers with magnitudes in the
#' range 0 to  \eqn{2^{n} - 1}. Elements of the input \code{u} outside of the
#' range -1 to 1 are treated as overflows and are saturated.
#' \itemize{
#'   \item For entries in the input u that are less than -1, the value of the
#'   output of uencode is 0.
#'   \item For entries in the input u that are greater than 1, the value of the
#'   output of uencode is \eqn{2^{n}-1}.
#' }
#'
#' \code{y <- uencode(u, n, v)} allows the input \code{u} to have entries with
#' floating-point values in the range \code{-v} to \code{v} before saturating
#' them (the default value for \code{v} is 1). Elements of the input \code{u}
#' outside of the range \code{-v} to \code{v} are treated as overflows and are
#' saturated:
#' \itemize{
#'   \item For input entries less than \code{-v}, the value of the output of
#'   uencode is 0.
#'   \item For input entries greater than \code{v}, the value of the output of
#'   uencode is \eqn{2^{n} - 1}.
#' }
#'
#' \code{y <- uencode(u, n, v, signed)} maps entries in a multidimensional array
#' of floating-point numbers \code{u} whose entries have values in the range
#' \code{-v} to \code{v} to an integer output \code{y}. Input entries outside
#' this range are saturated. The integer type of the output depends on the
#' number of quantization levels \eqn{2^{n}} and the value of \code{signed},
#' which can be one of the following:
#' \itemize{
#'   \item TRUE: Outputs are signed integers with magnitudes in the range
#'   \eqn{-2^{n} / 2} to \eqn{(2^{n} / 2) - 1}.
#'   \item FALSE (default): Outputs are unsigned integers with magnitudes in the
#'   range 0 to \eqn{2^{n} - 1}.
#' }
#'
#' @param u Input, a multidimensional array of numbers, real or complex, single
#'   or double precision.
#' @param n Number of levels used in \eqn{2^{n}}-level quantization. \code{n}
#'   must be between 2 and 32
#' @param v Limit on the range of \code{u} to the range from \code{-v} to
#'   \code{v} before saturating them. Default 1.
#' @param signed Logical indicating signed or unsigned output. See Details.
#'   Default: FALSE.
#'
#' @return Multidimensional array of the same size as \code{u} containing signed
#'   or unsigned integers.
#'
#' @examples
#'
#' u <- seq(-1, 1, 0.01)
#' y <- uencode(u, 3)
#' plot(u, y)
#'
#' @author Georgios Ouzounis, \email{ouzounis_georgios@@hotmail.com}.\cr
#' Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#'
#' @export

uencode <- function(u, n, v = 1, signed = FALSE) {

  if (!isScalar(n) || n < 2 || n > 32 || !isWhole(n))
    stop("n must be an integer in the range 2 to 32")
  if (!isPosscal(v) || !isWhole(v) || v <= 0)
    stop("v must be a positive integer")
  if (!is.logical(signed))
    stop("signed must be a logical")

  # R does not do comparisons with complex numbers (that makes some sense).
  # In Matlab, only the real parts of complex numbers are used for comparison.
  # This is done here as well. It would perhaps also make sense
  # to use the magnitide of the numbers to be compared...
  if (is.complex(u)) u <- Re(u)

  y <- u
  y[] <- 0

  width <- 2 * v / 2 ^ n

  y[u >= v] <- (2^n) - 1
  idx <- (u > -v) & (u < v)
  y[idx] <- floor((u[idx] + v) / width)
  if (signed) y <- y - 2 ^ (n - 1)
  y
}
