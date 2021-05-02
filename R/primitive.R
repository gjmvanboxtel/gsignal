# primitive.R
# Copyright (C) 2020 Geert van Bxtel <G.J.M.vanBoxtel@gmail.com>
# Original OCtave code:
# Copyright (C) 2013 Juan Pablo Carbajal
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
# 20201123  GvB       setup for gsignal v0.1.0
#------------------------------------------------------------------------------

#' Primitive
#'
#' Calculate the indefinitive integral of a function.
#'
#' This function is a fancy way of calculating the cumulative sum.
#'
#' @param FUN the function to calculate the primitive of.
#' @param t points at which the function \code{FUN} is evaluated, specified as a
#'   vector of ascending values
#' @param C constant of integration. Default: 0
#'
#' @return Vector of integrated function values.
#'
#' @examples
#' f <- function(t) sin(2 * pi * 3 * t)
#' t <- c(0, sort(runif(100)))
#' F <- primitive (f, t, 0)
#' t_true <- seq(0, 1, length.out = 1e3)
#' F_true <- (1 - cos(2 * pi * 3 * t_true)) / (2 * pi * 3)
#' plot (t, F, xlab = "", ylab = "")
#' lines (t_true, F_true, col = "red")
#' legend("topright", legend = c("Numerical primitive", "True primitive"),
#'   lty = c(0, 1), pch = c(1, NA), col = 1:2)
#'
#' @author Juan Pablo Carbajal.\cr
#' Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}
#'
#' @seealso \code{\link{cumsum}}
#
#' @export

primitive <- function(FUN, t, C = 0) {

  FUN <- match.fun(FUN)
  if (!is.vector(t) || !all(sort(t) == t)) {
    stop("t must be a vector of ascending values")
  }
  if (!isScalar(C)) {
    stop("C must be a scalar")
  }

  F_prev <- t0 <- NULL
  i_chunk <- function(t, f, init) {
    if (is.null(init)) {
      F_prev <<- NULL
      t0     <<- 0
    } else if (is.null(F_prev)) {
      F_prev <<- init
    } else {
      F_prev <<- F_prev + pracma::quadgk(f, t0, t)
      t0     <<- t
    }
    F <- F_prev
    invisible(F)
  }

  i_chunk(0, 0, NULL)
  y <- unlist(lapply(t, i_chunk, f = FUN, init = C))
  y
}
