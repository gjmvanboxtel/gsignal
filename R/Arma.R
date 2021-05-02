# Arma.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Original 'signal' version:
# Copyright (C) 2006 EPRI Solutions, Inc.
# by Tom Short, tshort@eprisolutions.com
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
# 20200127  GvB       setup for gsignal v0.1.0
# 20200402  GvB       Adapted to Octave filter conversion functions
# 20200425  GvB       as.Arma.Zpg(): adapted zero pole gain to z p g
# 20210603  GvB       as.Arma.Sos(): g = x$g instead of g = 1
#------------------------------------------------------------------------------

#' Autoregressive moving average (ARMA) model
#'
#' Create an ARMA model representing a filter or system model, or
#' convert other forms to an ARMA model.
#'
#' The ARMA model is defined by:
#' \deqn{a(L)y(t) = b(L)x(t)}
#' The ARMA model can define an analog or digital model. The AR and MA
#' polynomial coefficients follow the convention in 'Matlab' and 'Octave' where
#' the coefficients are in decreasing order of the polynomial (the opposite of
#' the definitions for \code{\link[stats]{filter}}filter and
#' \code{\link[base]{polyroot}}). For an analog model,
#' \if{latex}{
#'  \deqn{H(s) = (b_1 s^{(m-1)} + b_2 s^{(m-2)} + \ldots + b_m) / (a_1 s^{(n-1)}
#'  + a_2 s^{(n-2)} + \ldots + a_n)}
#' }
#' \if{html}{\preformatted{
#'   H(s) = (b[1]*s^(m-1) + b[2]*s^(m-2) + ... + b[m]) / (a[1]*s^(n-1) +
#'   a[2]*s^(n-2) + ... + a[n])
#' }}
#' For a z-plane digital model,
#' \if{latex}{
#'  \deqn{H(z) = (b_1 + b_2 z^{-1} + \ldots + b_m z^{(-m+1)}) / (a_1 + a_2
#'  z^{-1} + \ldots + a_n z^{(-n+1)})}
#' }
#' \if{html}{\preformatted{
#'   H(z) = (b[1] + b[2]*z^(-1) + … + b[m]*z^(-m+1)) / (a[1] + a[2]*z^(-1) + … +
#'   a[n]*z^(-n+1))
#' }}
#'
#'
#' \code{as.Arma} converts from other forms, including \code{Zpg} and \code{Ma}.
#'
#' @param b moving average (MA) polynomial coefficients.
#' @param a autoregressive (AR) polynomial coefficients.
#' @param x model or filter to be converted to an ARMA representation.
#' @param ...	additional arguments (ignored).
#'
#' @return A list of class \code{'Arma'} with the following list elements:
#' \describe{
#'   \item{b}{moving average (MA) polynomial coefficients}
#'   \item{a}{autoregressive (AR) polynomial coefficients}
#' }
#'
#' @seealso See also \code{\link{Zpg}}, \code{\link{Ma}}, \code{\link{filter}},
#'   and various filter-generation functions like \code{\link{butter}} and
#'   \code{\link{cheby1}} that return Arma models.
#'
#' @examples
#' filt <- Arma(b = c(1, 2, 1)/3, a = c(1, 1))
#' zplane(filt)
#'
#' @author Tom Short, \email{tshort@@eprisolutions.com},\cr
#' adapted by Geert van Boxtel, \email{gjmvanboxtel@@gmail.com}.
#'
#' @rdname Arma
#' @export

Arma <- function(b, a) {
  res <- list(b = b, a = a)
  class(res) <- "Arma"
  res
}

#' @rdname  Arma
#' @export
as.Arma <- function(x, ...) UseMethod("as.Arma")

#' @rdname Arma
#' @usage
#' ## S3 method for class 'Arma'
#' as.Arma(x, ...)
#' @export
as.Arma.Arma <- function(x, ...) x

#' @rdname Arma
#' @usage
#' ## S3 method for class 'Ma'
#' as.Arma(x, ...)
#' @export
as.Arma.Ma <- function(x, ...) {
  Arma(b = unclass(x), a = 1)
}

#' @rdname Arma
#' @usage
#' ## S3 method for class 'Sos'
#' as.Arma(x, ...)
#' @export
as.Arma.Sos <- function(x, ...) {
  ba <- sos2tf(x$sos, x$g)
  Arma(ba$b, ba$a)
}

#' @rdname Arma
#' @usage
#' ## S3 method for class 'Zpg'
#' as.Arma(x, ...)
#' @export
as.Arma.Zpg <- function(x, ...) {
  ba <- zp2tf(x$z, x$p, x$g)
  Arma(ba$b, ba$a)
}
