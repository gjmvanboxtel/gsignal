# tf2sos.R
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
# 20200402  GvB       setup for gsignal v0.1.0
# 20200406  GvB       validated
#------------------------------------------------------------------------------

#' Transfer function to second-order sections form
#'
#' Convert digital filter transfer function data to second-order section form.
#'
#' @param b moving average (MA) polynomial coefficients
#' @param a autoregressive (AR) polynomial coefficients
#'
#' @return A list with the following list elements:
#' \describe{
#'   \item{sos}{Second-order section representation, specified as an nrow-by-6
#'   matrix, whose rows contain the numerator and denominator coefficients of
#'   the second-order sections:\cr \code{sos <- rbind(cbind(B1, A1), cbind(...),
#'   cbind(Bn, An))}, where \code{B1 <- c(b0, b1, b2)}, and \code{A1 <- c(a0,
#'   a1, a2)} for section 1, etc. The b0 entry must be nonzero for each
#'   section.}
#' \item{g}{Overall gain factor that effectively scales the output \code{b}
#'   vector (or any one of the input \code{Bi} vectors).}
#' }
#'
#' @seealso See also \code{\link{filter}}
#'
#' @examples
#' b <- c(1, 0, 0, 0, 0, 1)
#' a <- c(1, 0, 0, 0, 0, .9)
#' sosg <- tf2sos (b, a)
#'
#' @author Julius O. Smith III, \email{jos@@ccrma.stanford.edu}.\cr
#' Conversion to R by Geert van Boxtel, \email{gjmvanboxtel@@gmail.com}.
#'
#' @export

tf2sos <- function(b, a) {

  zpk <- tf2zp(b, a)
  sos <- zp2sos(zpk$z, zpk$p, zpk$g)
  sos
}
