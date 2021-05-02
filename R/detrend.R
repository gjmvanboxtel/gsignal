# detrend.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Octave signal package:
# Copyright (C) 1995-2017 Kurt Hornik
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
# 20200104  GvB       setup for gsignal v0.1.0
#------------------------------------------------------------------------------

#' Remove Polynomial Trend
#'
#' \code{detrend} removes the polynomial trend of order \code{p} from the data
#' \code{x}.
#'
#' @param x Input vector or matrix. If \code{x} is a matrix, the trend is
#'   removed from the columns.
#' @param p Order of the polynomial. Default: 1. The order of the polynomial can
#'   also be given as a string, in which case \code{p} must be either
#'   \code{"constant"} (corresponds to \code{p = 0}) or \code{"linear"}
#'   (corresponds to \code{p = 1}).
#'
#' @return The detrended data, of same type and dimensions as \code{x}
#'
#' @examples
#' t <- 0:20
#' x <- 3 * sin(t) + t
#' y <- detrend(x)
#' plot(t, x, type = "l", ylim = c(-5, 25), xlab = "", ylab = "")
#' lines(t, y, col = "red")
#' lines(t, x - y, lty = 2)
#' legend('topleft', legend = c('Input Data', 'Detrended Data', 'Trend'),
#'  col = c(1, 2 ,1), lty = c(1, 1, 2))
#'
#' @author Kurt Hornik, \email{Kurt.Hornik@@wu-wien.ac.at}.\cr
#' Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#
#' @export

detrend <- function(x, p = 1) {

  if (!is.numeric(x)) {
    stop("input argument x must be a numeric vector, array or matrix")
  } else {
    x <- as.matrix(x)
  }

  if (is.character(p)) {
    p <- match.arg(p, c("constant", "linear"))
    if (p == "constant") {
      p <- 0
    } else if (p == "linear") {
      p <- 1
    } else {
      stop(paste("input argument p must be 'constant',",
                  "'linear', or a positive integer"))
    }
  } else {
    if (!isPosscal(p) || !isWhole(p)) {
      stop(paste("input argument p must be 'constant',",
                  "'linear', or a positive integer"))
    }
  }

  dims <- dim(x)
  if (dims[1] == 1) {
    x <- t(x)
  }

  r <- nrow(x)
  b <- (as.matrix(1:r) %*% matrix(1L, 1, (p + 1))) ^
     (matrix(1L, r, 1) %*% as.matrix(t(0:p)))
  y <- x - b %*% pracma::mldivide(b, x)

  class(y) <- class(x)
  attributes(y) <- attributes(x)
  y
}
