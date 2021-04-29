# Zpg.R
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
# 20200127  GvB       setup for gsignal v0.1.0
# 20200402  GvB       adapted to Octave filter conversion functions
# 20200406  GvB       change parameter names to z, p, g; added example
#------------------------------------------------------------------------------

#' Zero pole gain model
#'
#' Create an zero pole gain model of an ARMA filter, or convert other forms to a
#' Zpg model.
#'
#' \code{as.Zpg} converts from other forms, including \code{Arma} and \code{Ma}.
#'
#' @param z complex vector of the zeros of the model.
#' @param p complex vector of the poles of the model.
#' @param g overall gain of the model.
#' @param x model to be converted.
#' @param ...	additional arguments (ignored).
#'
#' @return A list of class Zpg with the following list elements:
#' \describe{
#'   \item{z}{complex vector of the zeros of the model}
#'   \item{p}{complex vector of the poles of the model}
#'   \item{g}{gain of the model}
#' }
#'
#' @seealso See also \code{\link{Arma}}
#'
#' @examples
#' ## design notch filter at pi/4 radians = 0.5/4 = 0.125 * fs
#' w = pi/4
#' # 2 poles, 2 zeros
#' # zeroes at r = 1
#' r <- 1
#' z1 <- r * exp(1i * w)
#' z2 <- r * exp(1i * -w)
#' # poles at r = 0.9
#' r = 0.9
#' p1 <- r * exp(1i * w)
#' p2 <- r * exp(1i * -w)
#'
#' zpg <- Zpg(c(z1, z2), c(p1, p2), 1)
#' zplane(zpg)
#' freqz(zpg)
#'
#' ## Sharper edges: increase distance between zeros and poles
#' r = 0.8
#' p1 <- r * exp(1i * w)
#' p2 <- r * exp(1i * -w)
#' zpg <- Zpg(c(z1, z2), c(p1, p2), 1)
#' zplane(zpg)
#' freqz(zpg)
#'
#' @author Tom Short, \email{tshort@@eprisolutions.com},\cr
#'  adapted by Geert van Boxtel, \email{gjmvanboxtel@@gmail.com}.
#' @rdname Zpg
#' @export

Zpg <- function(z, p, g) {
  res <- list(z = z, p = p, g = g)
  class(res) <- "Zpg"
  res
}

#' @rdname Zpg
#' @export
as.Zpg <- function(x, ...) UseMethod("as.Zpg")

#' @rdname Zpg
#' @usage
#' ## S3 method for class 'Arma'
#' as.Zpg(x, ...)
#' @export
as.Zpg.Arma <- function(x, ...) {
  zpk <- tf2zp(x$b, x$a)
  Zpg(zpk$z, zpk$p, zpk$k)
}

#' @rdname Zpg
#' @usage
#' ## S3 method for class 'Ma'
#' as.Zpg(x, ...)
#' @export
as.Zpg.Ma <- function(x, ...) {
  as.Zpg(as.Arma(x))
}

#' @rdname Zpg
#' @usage
#' ## S3 method for class 'Sos'
#' as.Zpg(x, ...)
#'
#' @export
as.Zpg.Sos <- function(x, ...) {

  zpk <- sos2zp(x$sos, x$g)
  Zpg(zpk$z, zpk$p, zpk$k)
}

#' @rdname Zpg
#' @usage
#' ## S3 method for class 'Zpg'
#' as.Zpg(x, ...)
#' @export
as.Zpg.Zpg <- function(x, ...) x
