# conv.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Octave version: Copyright (C) 1994-2017 John W. Eaton
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
# 2020209  GvB       setup for gsignal v0.1.0
#------------------------------------------------------------------------------

#' Convolution and polynomial multiplication
#'
#' Convolve two vectors \code{a} and \code{b}.
#'
#' The convolution of two vectors, \code{a} and \code{b}, represents the area of
#' overlap under the points as \code{B} slides across \code{a}. Algebraically,
#' convolution is the same operation as multiplying polynomials whose
#' coefficients are the elements of \code{a} and \code{b}.
#'
#' The function \code{conv} uses the \code{\link{filter}} function, NOT
#' \code{fft}, which may be faster for large vectors.
#'
#' @param a,b Input, coerced to vectors, can be different lengths or data types.
#' @param shape Subsection of convolution, partially matched to \code{"full"}
#'   (full convolution - default), \code{"same"} (central part of the
#'   convolution of the same size as \code{a}), or \code{"valid"} (only those
#'   parts of the convolution that are computed without the zero-padded edges)
#'
#' @return Output vector with length equal to \code{length (a) + length (b) -
#'   1}. When the parameter \code{shape} is set to \code{"valid"}, the length of
#'   the output is \code{max(length(a) - length(b) + 1, 0)}, except when
#'   length(b) is zero. In that case, the length of the output vector equals
#'   \code{length(a)}.
#'
#'   When \code{a} and \code{b} are the coefficient vectors of two polynomials,
#'   the convolution represents the coefficient vector of the product
#'   polynomial.
#'
#' @examples
#' u <- rep(1L, 3)
#' v <- c(1, 1, 0, 0, 0, 1, 1)
#' w <- conv(u, v)
#'
#' ## Create vectors u and v containing the coefficients of the polynomials
#' ## x^2 + 1 and 2x + 7.
#' u <- c(1, 0, 1)
#' v <- c(2, 7)
#' ## Use convolution to multiply the polynomials.
#' w <- conv(u, v)
#' ## w contains the polynomial coefficients for 2x^3 + 7x^2 + 2x + 7.
#'
#' ## Central part of convolution
#' u <- c(-1, 2, 3, -2, 0, 1, 2)
#' v <- c(2, 4, -1, 1)
#' w <- conv(u, v, 'same')
#'
#' @author Tony Richardson, \email{arichard@@stark.cc.oh.us}, adapted by John W.
#'   Eaton.\cr Conversion to R by Geert van Boxtel,
#'   \email{G.J.M.vanBoxtel@@gmail.com}.
#'
#' @export

conv <- function(a, b, shape = c("full", "same", "valid")) {

  a <- as.vector(a)
  b <- as.vector(b)
  shape <- match.arg(shape)

  la <- la_orig <- length(a)
  lb <- lb_orig <- length(b)

  ly <- la + lb - 1

  if (ly == 0) {
    y <- NULL
  } else {

    ## Use shortest vector as the coefficient vector to filter.
    if (la > lb) {
      tmp <- a
      a <- b
      b <- tmp
      lb <- la
    }
    x <- b

    ## Pad longer vector to convolution length.
    if (ly > lb) {
      x <- postpad(x, ly)
    }

    y <- filter(a, 1, x)

    if (shape == "same") {
      idx <- ceiling((ly - la) / 2)
      y <- y[(idx + 1):(idx + la)]
    } else if (shape == "valid") {
      len <- la_orig - lb_orig
      if (lb_orig + len < lb_orig) {
        y <- NULL
      } else {
        y <- y[(lb_orig:(lb_orig + len))]
      }
    }
  }
  y
}
