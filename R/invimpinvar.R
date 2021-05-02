# invimpinvar.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Original Octave function:
# Copyright (C) 2007 R.G.H. Eschauzier <reschauzier@yahoo.com>
# Copyright (C) 2011 Carne Draug <carandraug+dev@gmail.com>
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
# 20200622  GvB       setup for gsignal v0.1.0
#------------------------------------------------------------------------------

#' Inverse impulse invariance method
#'
#' Convert digital filter with coefficients b and a to analog, conserving
#' impulse response.
#'
#' Because \code{invimpinvar} is generic, it can also accept input of class
#' \code{\link{Arma}}.
#'
#' @param b coefficients of numerator polynomial
#' @param a coefficients of denominator polynomial
#' @param fs sampling frequency (Default: 1 Hz)
#' @param tol tolerance. Default: 0.0001
#' @param ... additional arguments (not used)
#'
#' @return A list of class \code{\link{Arma}} containing numerator and
#'   denominator polynomial filter coefficients of the A/D converted filter.
#'
#' @examples
#' f <- 2
#' fs <- 10
#' but <- butter(6, 2 * pi * f, 'low', 's')
#' zbut <- impinvar(but, fs)
#' sbut <- invimpinvar(zbut, fs)
#' all.equal(but, sbut, tolerance = 1e-7)
#'
#' @author R.G.H. Eschauzier, \email{reschauzier@@yahoo.com},\cr
#'  Carne Draug, \email{carandraug+dev@@gmail.com}.\cr
#'  Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}
#'
#' @seealso \code{\link{impinvar}}
#'
#' @references Thomas J. Cavicchi (1996) Impulse invariance and multiple-order
#'   poles. IEEE transactions on signal processing, Vol 40 (9): 2344--2347.
#'
#' @rdname invimpinvar
#' @export

invimpinvar <- function(b, ...) UseMethod("invimpinvar")

#' @rdname invimpinvar
#' @export

invimpinvar.Arma <- function(b, ...)
  invimpinvar(b$b, b$a, ...)

#' @rdname invimpinvar
#' @export

invimpinvar.default <- function(b, a, fs = 1, tol = 0.0001, ...) {

  if (!isPosscal(fs)) {
    stop("fs must be a positive scalar")
  }
  if (!isPosscal(tol)) {
    stop("tol must be a positive scalar")
  }
  ts <- 1 / fs

  b <- c(b, 0)
  rpk_in <- residue(b, a)
  n <- length(rpk_in$r)

  if (length(rpk_in$k) > 1) {
    stop("Order numerator > order denominator")
  }

  r_out  <- rep(0L, n)
  sm_out <- rep(0L, n)

  i <- 1
  while (i <= n) {
    m <- 1
    first_pole <- rpk_in$p[i]
    while (i < n && abs(first_pole - rpk_in$p[i + 1]) < tol) {
      i <- i + 1
      m <- m + 1
    }
    rpk_out      <- inv_z_res(rpk_in$r[(i - m + 1):i], first_pole, ts)
    rpk_in$k <- rpk_in$k - rpk_out$k
    sm_out[(i - m + 1):i] <- rpk_out$p
    r_out[(i - m + 1):i]  <- rpk_out$r

    i <- i + 1
  }
  ba <- inv_residue(r_out, sm_out, 0, tol)
  a    <- zapIm(ba$a)
  b    <- zapIm(ba$b)
  b <- polyreduce(zapsmall(b))
  Arma(b, a)
}

## Inverse function of z_res (see impinvar source)
inv_z_res <- function(r_in, p_in, ts) {

  n    <- length(r_in)
  r_out <- rep(0L, n)

  j <- n
  while (j > 1) {
    r_out[j]   <- r_in[j] / ((ts * p_in)^j)
    r_in[1:j] <- r_in[1:j] - r_out[j] * rev(h1_z_deriv(j - 1, p_in, ts))
    j <- j - 1
  }

  r_out[1] <- r_in[1] / ((ts * p_in))
  k_out    <- r_in[1] / p_in
  sm_out   <- log(p_in) / ts

  list(r = r_out, p = sm_out, k = k_out)
}

# Source code of h1_deriv and h1_z_deriv in impinvar.R
