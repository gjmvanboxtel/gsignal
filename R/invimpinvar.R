# invimpinvar.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Original Octave function:
# Copyright (C) 2007 R.G.H. Eschauzier <reschauzier@yahoo.com>
# Copyright (C) 2011 Carnë Draug <carandraug+dev@gmail.com>
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
# 20200622  GvB       setup for gsignal v0.1.0
#---------------------------------------------------------------------------------------------------------------------

#' Inverse impulse invariance method
#' 
#' Converts digital filter with coefficients b and a to analog, conserving impulse response.
#' 
#' This function does the inverse of impinvar so that the following example
#' should restore the original values of a and b.
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
#' @author Original Octave version R.G.H. Eschauzier
#'   \email{reschauzier@@yahoo.com}, Carnë Draug
#'   \email{carandraug+dev@@gmail.com}. Conversion to R by Geert van Boxtel
#'   \email{G.J.M.vanBoxtel@@gmail.com}
#'   
#' @seealso \code{\link{impinvar}}
#' 
#' @references Thomas J. Cavicchi (1996) ``Impulse invariance and multiple-order
#'   poles''. IEEE transactions on signal processing, Vol 40 (9): 2344--2347.
#'
#' @rdname invimpinvar
#' @export

invimpinvar <- function(b, ...) UseMethod("invimpinvar")

#' @rdname invimpinvar
#' @export

invimpinvar.Arma <- function(b, ...)
  invimpinvar (b$b, b$a, ...)

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

  b <- c(b, 0)                  # so we can calculate in z instead of z^-1
  rpk_in <- residue(b, a)       # partial fraction expansion
  n <- length(rpk_in$r)         # Number of poles/residues
  
  if (length(rpk_in$k) > 1) {   # Greater than one means we cannot do impulse invariance
    stop("Order numerator > order denominator")
  }
  
  r_out  <- rep(0L, n)          # Residues of H(s)
  sm_out <- rep(0L, n)          # Poles of H(s)
  
  i <- 1
  while (i <= n) {
    m <- 1
    first_pole <- rpk_in$p[i]                                        # Pole in the z-domain
    while (i < n && abs(first_pole - rpk_in$p[i + 1]) < tol) {       # Multiple poles at p(i)
      i <- i + 1    # Next residue
      m <- m + 1    # Next multiplicity
    }
    rpk_out      <- inv_z_res(rpk_in$r[(i - m + 1):i], first_pole, ts)   # Find s-domain residues
    rpk_in$k <- rpk_in$k - rpk_out$k                                     # Just to check, should end up zero for physical system
    sm_out[(i - m + 1):i] <- rpk_out$p                                   # Copy s-domain pole(s) to output
    r_out[(i - m + 1):i]  <- rpk_out$r                                   # Copy s-domain residue(s) to output
  
    i <- i + 1 # Next z-domain residue/pole
  }
  ba <- inv_residue(r_out, sm_out , 0, tol)
  a    <- zapIm(ba$a)                                        # Get rid of spurious imaginary part
  b    <- zapIm(ba$b)
  b <- polyreduce(zapsmall(b))  
  Arma(b, a)
}

## Inverse function of z_res (see impinvar source)

inv_z_res <- function(r_in, p_in, ts) {

  n    <- length(r_in)                                                     # multiplicity of the pole
  r_out <- rep(0L, n)

  j <- n
  while (j > 1) {                                                          # Go through residues starting from highest order down
    r_out[j]   <- r_in[j] / ((ts * p_in)^j)                                # Back to binomial coefficient for highest order (always 1)
    r_in[1:j] <- r_in[1:j] - r_out[j] * rev(h1_z_deriv(j -1 , p_in, ts))   # Subtract highest order result, leaving r_in(j) zero
    j <- j - 1
  }

  ## Single pole (no multiplicity)
  r_out[1] <- r_in[1] / ((ts * p_in))
  k_out    <- r_in[1] / p_in
  sm_out   = log(p_in) / ts

  list(r = r_out, p = sm_out, k = k_out)
}

# Source code of h1_deriv and h1_z_deriv in impinvar.R