# sos2tf.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Original Octave version:
# Copyright (C) 2005 Julius O. Smith III <jos@ccrma.stanford.edu>
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
# 20200330  GvB       setup for gsignal v0.1.0
# 20200406  GvB       validated
# 20210306  GvB       initialize a, b with sos[1, ] instead of 1 (bug in Octave
#                     signal?)
# 20210326  GvB       return object of class 'Arma'
# 20210506  GvB       use matrix() instead of as.matrix()
#------------------------------------------------------------------------------

#' Sos to transfer function
#'
#' Convert digital filter second-order section data to transfer function form.
#'
#' @param sos Second-order section representation, specified as an nrow-by-6
#'   matrix, whose rows contain the numerator and denominator coefficients of
#'   the second-order sections:\cr \code{sos <- rbind(cbind(B1, A1), cbind(...),
#'   cbind(Bn, An))}, where \code{B1 <- c(b0, b1, b2)}, and \code{A1 <- c(a0,
#'   a1, a2)} for section 1, etc. The b0 entry must be nonzero for each section.
#' @param g Overall gain factor that effectively scales the output \code{b}
#'   vector (or any one of the input \code{Bi} vectors). Default: 1.
#'
#' @return An object of class "Arma" with the following list elements:
#' \describe{
#'   \item{b}{moving average (MA) polynomial coefficients}
#'   \item{a}{autoregressive (AR) polynomial coefficients}
#' }
#'
#' @seealso \code{\link{as.Arma}}, \code{\link{filter}}
#'
#' @examples
#' sos <- rbind(c(1, 1, 1, 1, 0, -1), c(-2, 3, 1, 1, 10, 1))
#' ba <- sos2tf(sos)
#'
#' @author Julius O. Smith III, \email{jos@@ccrma.stanford.edu}.\cr
#' Conversion to R by Geert van Boxtel, \email{gjmvanboxtel@@gmail.com}.
#'
#' @export

sos2tf <- function(sos, g = 1) {

  sos <- matrix(sos, ncol = 6)
  n <- nrow(sos)
  m <- ncol(sos)
  if (n <= 0) {
    stop("sos must have at least 1 row")
  }
  if (m != 6) {
    stop("sos must be a nrow-by-6 matrix")
  }

  b <- sos[1, 1:3]
  a <- sos[1, 4:6]

  if (n > 1) {
    for (i in 2:n) {
      b <- conv(b, sos[i, 1:3])
      a <- conv(a, sos[i, 4:6])
    }
  }

  nb <- length(b)
  while (nb > 0 && b[nb] == 0) {
    b <- b[1:(nb - 1)]
    nb <-  length(b)
  }

  na <- length(a)
  while (na > 0 && a[na] == 0) {
    a <- a[1:(na - 1)]
    na <- length(a)
  }

  b <-  b * prod(g)

  Arma(b = b, a = a)
}
