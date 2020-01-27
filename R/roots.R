# roots.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Original Octave function Copyright (C) 1994-2017 John W. Eaton
# Author: KH <Kurt.Hornik@wu-wien.ac.at>
# Created: 24 December 1993
# Adapted-By: jwe
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
# 20200127  GvB       setup for gsignal v0.1.0
#---------------------------------------------------------------------------------------------------------------------

#' Polynomial roots
#' 
#' Compute the roots of a polynomial represented by \code{p}
#' 
#' The \code{roots} function solves polynomial equations of the form \code{p(1)
#' * x^(n-1) + ... + p(n-1) * x + p(n) = 0}. Polynomial equations contain a
#' single variable with nonnegative exponents.
#' 
#' 
#' @param p Polynomial coefficients with coefficients given in order from
#'   highest to lowest polynomial power. This is the Matlab/Octave convention;
#'   it is opposite of the convention used by polyroot.
#' @param method Either “polyroot” (default) which uses \code{polyroot} for its
#'   computations internally (and should be more accurate) or “eigen” which
#'   uses eigenvalues of the companion matrix for its computation. The latter
#'   returns complex values in case of real valued solutions in less cases.
#' 
#' @return A complex vector with the roots of the polynomial (converted to
#'   numeric if all imaginary components are 0).
#' 
#' @examples
#' roots(1:3)
#' polyroot(3:1) # should be the same
#' poly(roots(1:3))
#' 
#' roots(1:3, method="eigen") # using eigenvalues
#' 
#' @seealso \code{\link{poly}}
#' 
#' @author Original Octave version by Kurt Hornik. Conversion to R by Tom Short.
#'   Adapted by Geert van Boxtel
#
#' @export

roots <- function(p, method = c("polyroot", "eigen")) {
  
  method <- match.arg(method)
  
  p <- as.vector(p)
  if(!(is.numeric(p) || is.complex(p))) {
    stop("p must be a numeric or complex vector")
  }

  n <- length(p)
  maxp <- max(abs(p))
  if (is.null(p) || maxp == 0) {
    r <- NULL
  } else {
    if(method == "polyroot") {
      r <- base::polyroot(rev(p))
    } else if (method == "eigen") {
      f <- which((p / maxp) != 0)
      m <- length(f)
      p <- p[f[1]:f[m]]
      l <- length(p)
      if (l > 1) {
        A <- rbind((-p[2:l] / p[1]), diag(1, l - 1)[1:(l - 2),])   #companion matrix
        r <- eigen(A)$values                                       #roots are eigenvalues of cm
      } else {
        r <- rep(0L, n-f[m])
      }
    } else {
      stop("method should be 'polyroot' or 'eigen'")
    }
    if (all(Im(z <- zapsmall(r)) == 0)) r <- Re(r)
  }
  r[!is.na(r) & !is.infinite(r)]
}