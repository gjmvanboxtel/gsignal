# residuez.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Original Octave function:
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
# 20200804  GvB       setup for gsignal v0.1.0
#------------------------------------------------------------------------------

#' Z-transform partial fraction expansion
#'
#' Finds the residues, poles, and direct term of a Partial Fraction Expansion of
#' the ratio of two polynomials.
#'
#' \code{residuez} converts a discrete time system, expressed as the ratio of
#' two polynomials, to partial fraction expansion, or residue, form.
#'
#' @param b coefficients of numerator polynomial
#' @param a coefficients of denominator polynomial
#'
#' @return A list containing
#' \describe{
#'   \item{r}{vector of filter pole residues of the partial fraction}
#'   \item{p}{vector of partial fraction poles}
#'   \item{k}{vector containing FIR part, if any (empty if \code{length(b) <
#'   length(a)})}
#' }
#'
#' @seealso \code{\link{residue}}, \code{\link{residued}}
#'
#' @examples
#' b0 <- 0.05634
#' b1 <- c(1,  1)
#' b2 <- c(1, -1.0166, 1)
#' a1 <- c(1, -0.683)
#' a2 <- c(1, -1.4461, 0.7957)
#' b <- b0 * conv(b1, b2)
#' a <- conv(a1, a2)
#' res <- residuez(b, a)
#'
#' @author Julius O. Smith III, \email{jos@@ccrma.stanford.edu}.\cr
#'  Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}
#'
#' @export

residuez <- function(b, a) {

  if (!is.vector(b) || !is.vector(a)) {
    stop("b and a must be vectors")
  }

  rpk <- residue(rev(b), rev(a))
  p <- 1 / rpk$p
  m <- mpoles(p)
  r <- rpk$r * ((-p)^m)
  if (!is.null(rpk$k)) {
    k <- Conj(rev(rpk$k))
  } else {
    k <- NULL
  }

  r <- zapIm(r)
  p <- zapIm(p)
  list(r = r, p = p, k = k)
}
