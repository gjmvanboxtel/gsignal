# xcorr2.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Octave version:
# Copyright (C) 2000 Dave Cogdell <cogdelld@asme.org>
# Copyright (C) 2000 Paul Kienzle <pkienzle@users.sf.net>
# Copyright (C) 2012 Carne Draug <carandraug+dev@gmail.com>
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
# 2020313  GvB       setup for gsignal v0.1.0
#------------------------------------------------------------------------------

#' 2-D cross-correlation
#'
#' Compute the 2D cross-correlation of matrices \code{a} and \code{b}.
#'
#' If \code{b} is not specified, computes autocorrelation of \code{a},
#' i.e., same as \code{xcorr2 (a, a)}.
#'
#' @param a Input matrix, coerced to numeric. Must not be missing.
#' @param b Input matrix, coerced to numeric. Default: \code{a}.
#' @param scale Character string. Specifies the type of scaling applied to the
#'   correlation matrix. matched to one of:
#'   \describe{
#'     \item{"none"}{no scaling}
#'     \item{"biased"}{Scales the raw cross-correlation by the maximum number of
#'     elements of \code{a} and \code{b} involved in the generation of any
#'     element of the output matrix.}
#'     \item{"unbiased"}{Scales the raw correlation by dividing each element in
#'     the cross-correlation matrix by the number of products \code{a} and
#'     \code{b} used to generate that element. }
#'     \item{"coeff"}{Scales the normalized cross-correlation on the range of [0
#'     1] so that a value of 1 corresponds to a correlation coefficient of 1. }
#'   }
#'
#' @return 2-D cross-correlation or autocorrelation matrix, returned as a matrix
#'
#' @seealso \code{\link{conv2}}, \code{\link{xcorr}}.
#'
#' @examples
#' m1 <- matrix(c(17, 24,  1,  8, 15,
#'                23,  5,  7, 14, 16,
#'                 4,  6, 13, 20, 22,
#'                10, 12, 19, 21,  3,
#'                11, 18, 25,  2,  9), 5, 5, byrow = TRUE)
#' m2 <- matrix(c(8, 1, 6,
#'                3, 5, 7,
#'                4, 9, 2), 3, 3, byrow = TRUE)
#' R <- xcorr2(m1, m2)
#'
#' @author Dave Cogdell, \email{cogdelld@@asme.org},\cr
#' Paul Kienzle, \email{pkienzle@@users.sf.net},\cr
#' Carne Draug, \email{carandraug+dev@@gmail.com}.\cr
#' Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#'
#' @export

xcorr2 <- function(a, b = a,
                   scale = c("none", "biased", "unbiased", "coeff")) {

  if (!is.matrix(a) || !is.matrix(b)) {
    stop("input matrices must must have 2 dimensions")
  }
  scale <- match.arg(scale)

  ## compute correlation
  ma <- nrow(a)
  na <- ncol(a)
  mb <- nrow(b)
  nb <- ncol(b)
  R <- conv2(a, Conj(b[rev(1:mb), rev(1:nb)]))

  # bias routines by Dave Cogdell (cogdelld@asme.org)
  # optimized by Paul Kienzle (pkienzle@users.sf.net)
  # coeff routine by Carnë Draug (carandraug+dev@gmail.com)
  if (scale == "biased") {
    R <- R / (min(ma, mb) * min(na, nb))
  } else if (scale == "unbiased") {
    lo <- min(na, nb)
    hi <- max(na, nb)
    row  <- c(1:(lo - 1), rep(lo, (hi - lo + 1)), rev(1:(lo - 1)))
    lo  <- min(ma, mb)
    hi <- max(ma, mb)
    col <- c(1:(lo - 1), rep(lo, (hi - lo + 1)), rev(1:(lo - 1)))
    bias <- outer(col, row)
    R <- R / bias
  } else if (scale == "coeff") {
      a <- Re(a)
      b <- Re(b)
      a <- conv2(a^2, matrix(1L, nrow(b), ncol(b)))
      b <- ssq(as.vector(b))
      R <- R / sqrt(a * b)
  }
  R
}
