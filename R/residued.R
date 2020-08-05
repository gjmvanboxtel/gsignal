# residued.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Original Octave function:
# Copyright (C) 2005 Julius O. Smith III <jos@ccrma.stanford.edu>
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
# See also: http://www.gnu.org/licenses/gpl-2.0.txt
#
# Version history
# 20200805  GvB       setup for gsignal v0.1.0
#---------------------------------------------------------------------------------------------------------------------

#' delayed z-transform partial fraction expansion
#'
#' Finds the residues, poles, and direct term of a Partial Fraction Expansion of
#' the ratio of two polynomials.
#'
#' In the usual PFE function \code{residuez}, the IIR part (poles \code{p} and
#' residues \code{r}) is driven in parallel with the FIR part (\code{f}). In
#' this variant, the IIR part is driven by the output of the FIR part. This
#' structure can be more accurate in signal modeling applications.
#'
#' @param b coefficients of numerator polynomial
#' @param a coefficients of denominator polynomial
#'
#' @return A \code{\link{list}} containing
#' \describe{
#'   \item{r}{vector of filter pole residues of the partial fraction}
#'   \item{p}{vector of partial fraction poles}
#'  \item{k}{vector containing FIR part, if any (empty if \code{length(b) <
#'  length(a)})}
#' }
#'
#' @seealso \code{\link{residue}}, \code{\link{residuez}}
#'
#' @references \url{https://ccrma.stanford.edu/~jos/filters/residued.html}
#'
#' @examples
#' b <- c(2, 6, 6, 2)
#' a <- c(1, -2, 1)
#' resd <- residued(b, a)
#' resz <- residuez(b, a)
#'
#' @author Original Octave version by Julius O. Smith III
#'   \email{jos@@ccrma.stanford.edu}. Conversion to R by Geert van Boxtel
#'   \email{G.J.M.vanBoxtel@@gmail.com}
#'
#' @export

residued <- function(b, a) {
  
  ## RESIDUED - return residues, poles, and FIR part of B(z)/A(z)
  ##
  ## Let nb = length(b), na = length(a), and N=na-1 = no. of poles.
  ## If nb<na, then f will be empty, and the returned filter is
  ##
  ##             r(1)                      r(N)
  ## H(z) = ----------------  + ... + ----------------- = R(z)
  ##        [ 1-p(1)/z ]^e(1)         [ 1-p(N)/z ]^e(N)
  ##
  ## This is the same result as returned by RESIDUEZ.
  ## Otherwise, the FIR part f will be nonempty,
  ## and the returned filter is
  ##
  ## H(z) = f(1) + f(2)/z + f(3)/z^2 + ... + f(nf)/z^M + R(z)/z^M
  ##
  ## where R(z) is the parallel one-pole filter bank defined above,
  ## and M is the order of F(z) = length(f)-1 = nb-na.
  ##
  ## Note, in particular, that the impulse-response of the parallel
  ## (complex) one-pole filter bank starts AFTER that of the the FIR part.
  ## In the result returned by RESIDUEZ, R(z) is not divided by z^M,
  ## so its impulse response starts at time 0 in parallel with f(n).
  ##
  ## J.O. Smith, 9/19/05
  
  num <- Conj(b)
  den <- Conj(a)
  nb <- length(num)
  na <- length(den)
  k <- NULL
  if (na <= nb) {
    k <- filter(num, den, c(1L, rep(0L, nb - na)))
    num <- num - conv(den, k)
    num <- num[(nb - na + 2):nb]
  }
  rpk <- residuez(num, den)
  if (!is.null(rpk$k)) {
    stop('rpk$f not empty as expected')
  }
  
  r <- zapIm(rpk$r)
  p <- zapIm(rpk$p)
  list(r = r, p = p, k = k)
}
