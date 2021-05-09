# sos2zp.R
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
# 20200405  GvB       replaced roots() by pracma::roots()
# 20200406  GvB       validated
# 20210326  GvB       return object of class 'Zpg'
# 20210506  GvB       use matrix() instead of as.matrix(), sort output
#------------------------------------------------------------------------------

#' Sos to zero-pole-gain
#'
#' Convert digital filter second-order section data to zero-pole-gain form.
#'
#' @param sos Second-order section representation, specified as an nrow-by-6
#'   matrix, whose rows contain the numerator and denominator coefficients of
#'   the second-order sections:\cr \code{sos <- rbind(cbind(B1, A1), cbind(...),
#'   cbind(Bn, An))}, where \code{B1 <- c(b0, b1, b2)}, and \code{A1 <- c(a0,
#'   a1, a2)} for section 1, etc. The b0 entry must be nonzero for each section.
#' @param g Overall gain factor that effectively scales the output \code{b}
#'   vector (or any one of the input \code{B_i} vectors). Default: 1.
#'
#'@return A list of class "Zpg" with the following list elements:
#' \describe{
#'   \item{z}{complex vector of the zeros of the model (roots of \code{B(z)})}
#'   \item{p}{complex vector of the poles of the model (roots of \code{A(z)})}
#'   \item{k}{overall gain (\code{B(Inf)})}
#' }
#'
#' @seealso \code{\link{filter}}
#'
#' @examples
#' sos <- rbind(c(1, 0, 1, 1, 0, -0.81), c(1, 0, 0, 1, 0, 0.49))
#' zpk <- sos2zp(sos)
#'
#' @author Julius O. Smith III \email{jos@@ccrma.stanford.edu}.\cr
#' Conversion to R by, Geert van Boxtel \email{G.J.M.vanBoxtel@@gmail.com}
#'
#' @export

sos2zp <- function(sos, g = 1) {

  sos <- matrix(sos, ncol = 6)
  n <- nrow(sos)
  m <- ncol(sos)
  if (m != 6) {
    stop("sos must be a nrow-by-6 matrix")
  }

  gains <- sos[, 1]             # All b0 coeffs
  g <- prod(gains) * g          # pole-zero gain
  if (g == 0) {
    stop("one or more section gains is zero")
  }
  sos[, 1:3] <- sos[, 1:3] / c(gains, gains, gains)

  z <- p <- rep(0L, 2 * n)
  for (i in seq_len(n)) {
    ndx <- (2 * i - 1):(2 * i)
    zi <- pracma::roots(sos[i, 1:3])
    z[ndx] <- zi
    pi <- pracma::roots(sos[i, 4:6])
    p[ndx] <- pi
  }

  Zpg(z = sort(z), p = sort(p), g = g)
}
