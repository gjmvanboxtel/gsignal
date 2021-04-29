# impinvar.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Original Octave function:
# Copyright (C) 1994-2017 John W. Eaton
# Copyright (C) 2007 Ben Abbott
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
# 20200616  GvB       setup for gsignal v0.1.0
#------------------------------------------------------------------------------

#' Impulse invariance method for A/D filter conversion
#'
#' Convert analog filter with coefficients b and a to digital, conserving
#' impulse response.
#'
#' Because \code{impinvar} is generic, it can also accept input of class
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
#' freqz(zbut, n = 1024, fs = fs)
#'
#' @author Tony Richardson, \email{arichard@@stark.cc.oh.us},\cr
#'  Ben Abbott, \email{bpabbott@@mac.com},\cr
#'   adapted by John W. Eaton.\cr
#'   Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}
#'
#' @seealso \code{\link{invimpinvar}}
#'
#' @rdname impinvar
#' @export

impinvar <- function(b, ...) UseMethod("impinvar")

#' @rdname impinvar
#' @export

impinvar.Arma <- function(b, ...)
  impinvar(b$b, b$a, ...)

#' @rdname impinvar
#' @export

impinvar.default <- function(b, a, fs = 1, tol = 0.0001, ...) {

  if (!isPosscal(fs)) {
    stop("fs must be a positive scalar")
  }
  if (!isPosscal(tol)) {
    stop("tol must be a positive scalar")
  }
  ts <- 1 / fs

  rpk_in <- residue(b, a)
  n <- length(rpk_in$r)

  if (length(rpk_in$k) > 0) {
    stop("Order numerator >= order denominator")
  }

  r_out <- rep(0L, n)
  p_out <- rep(0L, n)
  k_out <- 0

  i <- 1
  while (i <= n) {
    m <- 1
    first_pole <- rpk_in$p[i]
    while (i < n && abs(first_pole - rpk_in$p[i + 1]) < tol) {
      i <- i + 1
      m <- m + 1
    }
    rpk_out <- z_res(rpk_in$r[(i - m + 1):i], first_pole, ts)
    k_out                <- k_out + rpk_out$k
    p_out[(i - m + 1):i] <- rpk_out$p
    r_out[(i - m + 1):i] <- rpk_out$r

    i <- i + 1
  }

  ba <- inv_residue(r_out, p_out, k_out, tol)
  a <- zapIm(ba$a)
  b <- zapIm(ba$b)

  b <- b[1:(length(b) - 1)]
  Arma(b, a)
}

z_res <- function(r_in, sm, ts) {

  p_out <- exp(ts * sm)
  n     <- length(r_in)
  r_out <- rep(0L, n)

  k_out    <- r_in[1] * ts
  r_out[1] <- r_in[1] * ts * p_out

  if (n > 1) {
    for (i in 2:n) {
      r_out[1:i] <- r_out[1:i] + r_in[i] *
        rev(h1_z_deriv(i - 1, p_out, ts))
    }
  }

  list(r = r_out, p = p_out, k = k_out)
}

# The following functions are
# Copyright (C) 2007 R.G.H. Eschauzier <reschauzier@yahoo.com>
# Conversion to R by Geert van Boxtel

h1_deriv <- function(n) {

  b  <- pracma::fact(n) *
    sapply(0:n, function(k) pracma::nchoosek(n, k))
  b  <- b * (-1)^n
  b
}

h1_z_deriv <- function(n, p, ts) {
  d <- (-1)^n
  for (i in 1:(n - 1)) {
    d <- c(d, 0)
    d <- d + prepad(pracma::polyder(d), i + 1, 0, 2)
  }
  b <- rep(0L, n + 1)
  for (i in 1:n) {
    b  <- b + d[i] * prepad(h1_deriv(n - i + 1), n + 1, 0, 2)
  }
  b <- b * ts ^ (n + 1) / pracma::fact(n)
  b <- b * p^seq(n + 1, 1, -1)
  b
}

inv_residue <- function(r_in, p_in, k_in, tol) {

  n <- length(r_in)
  k <- 0
  if (length(k_in) == 1) {
    k <- k_in[1]
  } else if (length(k_in) > 1) {
    stop("Order numerator > order denominator")
  }
  a_out <- poly(p_in)
  b_out  <- rep(0L, n + 1)
  b_out <- b_out + k * a_out
  i <- 1
  while (i <= n) {
    term   <- c(1, -p_in[i])
    p      <- r_in[i] * pracma::deconv(a_out, term)$q
    p      <- prepad(p, n + 1, 0, 2)
    b_out <- b_out + p
    m          <- 1
    mterm      <- term
    first_pole <- p_in[i]
    while (i < n && abs(first_pole - p_in[i + 1]) < tol) {
      i <- i + 1
      m <- m + 1
      mterm  <- conv(mterm, term)
      p      <- r_in[i] * pracma::deconv(a_out, mterm)$q
      p      <- prepad(p, n + 1, 0, 2)
      b_out  <- b_out + p
    }
    i <- i + 1
  }
  Arma(b_out, a_out)
}
