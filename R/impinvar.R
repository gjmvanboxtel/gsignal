# impinvar.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Original Octave function:
# Copyright (C) 1994-2017 John W. Eaton
# Copyright (C) 2007 Ben Abbott
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
# 20200616  GvB       setup for gsignal v0.1.0
#---------------------------------------------------------------------------------------------------------------------

#' Impulse invariance method for A/D filter conversion
#' 
#' Converts analog filter with coefficients b and a to digital, conserving
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
#' @author Original Octave version by Tony Richardson
#'   \email{arichard@@stark.cc.oh.us}, Ben Abbott \email{bpabbott@@mac.com},
#'   adapted by John W. Eaton. Conversion to R by Geert van Boxtel
#'   \email{G.J.M.vanBoxtel@@gmail.com}
#'   
#' @seealso \code{invimpinvar}
#'
#' @rdname impinvar
#' @export

impinvar <- function(b, ...) UseMethod("impinvar")

#' @rdname impinvar
#' @export

impinvar.Arma <- function(b, ...)
  impinvar (b$b, b$a, ...)

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

  rpk_in <- residue(b, a)                # partial fraction expansion
  n <- length(rpk_in$r)                  # Number of poles/residues
  
  if (length(rpk_in$k) > 0) {            # Greater than zero means we cannot do impulse invariance
    stop("Order numerator >= order denominator")
  }
  
  r_out <- rep(0L, n)                        # Residues of H(z)
  p_out <- rep(0L, n)                        # Poles of H(z)
  k_out <- 0                                 # Constant term of H(z)
  
  i <- 1
  while (i <= n) {
    m <- 1
    first_pole <- rpk_in$p[i]                                   # Pole in the s-domain
    while (i < n && abs(first_pole - rpk_in$p[i + 1]) < tol) {  # Multiple poles at p(i)
      i <- i + 1          # Next residue
      m <- m + 1          # Next multiplicity
    }
    rpk_out <- z_res(rpk_in$r[(i - m + 1):i], first_pole, ts)   # Find z-domain residues
    k_out                <- k_out + rpk_out$k                   # Add direct term to output
    p_out[(i - m + 1):i] <- rpk_out$p                           # Copy z-domain pole(s) to output
    r_out[(i - m + 1):i] <- rpk_out$r                           # Copy z-domain residue(s) to output
  
    i <- i + 1       # Next s-domain residue/pole
  }
  
  ba      <- inv_residue(r_out, p_out, k_out, tol)
  a    <- zapIm(ba$a)                                        # Get rid of spurious imaginary part
  b    <- zapIm(ba$b)
  
  ## Shift results right to account for calculating in z instead of z^-1
  b <- b[1:(length(b) - 1)]
  Arma(b, a)
}


## Convert residue vector for single and multiple poles in s-domain (located at sm) to
## residue vector in z-domain. The variable k is the direct term of the result.

z_res <- function (r_in, sm, ts) {

  p_out <- exp(ts * sm)        # z-domain pole
  n     <- length(r_in)        # Multiplicity of the pole
  r_out <- rep(0L, n)          # Residue vector

  ## First pole (no multiplicity)
  k_out    <- r_in[1] * ts         # PFE of z/(z-p) = p/(z-p)+1; direct part
  r_out[1] = r_in[1] * ts * p_out  # pole part of PFE

  if (n > 1) {
    for (i in 2:n) {           # Go through s-domain residues for multiple pole
      r_out[1:i] <- r_out[1:i] + r_in[i] * rev(h1_z_deriv(i - 1, p_out, ts)) # Add z-domain residues
    }
  }

  list(r = r_out, p = p_out, k = k_out)
}


# The following functions are # Copyright (C) 2007 R.G.H. Eschauzier <reschauzier@yahoo.com>
# Port to R by Geert van Boxtel

## Find (z^n)*(d/dz)^n*H1(z), where H1(z)=ts*z/(z-p), ts=sampling period,
## p=exp(sm*ts) and sm is the s-domain pole with multiplicity n+1.
## The result is (ts^(n+1))*(b(1)*p/(z-p)^1 + b(2)*p^2/(z-p)^2 + b(n+1)*p^(n+1)/(z-p)^(n+1)),
## where b(i) is the binomial coefficient bincoeff(n,i) times n!. Works for n>0.
h1_deriv <- function (n) {

  b  <- pracma::fact(n) * sapply(0:n, function(k) pracma::nchoosek(n, k))      # Binomial coefficients: [1], [1 1], [1 2 1], [1 3 3 1], etc.
  b  <- b * (-1)^n
  b
}

## Find {-zd/dz}^n*H1(z). I.e., first differentiate, then multiply by -z, then differentiate, etc.
## The result is (ts^(n+1))*(b(1)*p/(z-p)^1 + b(2)*p^2/(z-p)^2 + b(n+1)*p^(n+1)/(z-p)^(n+1)).
## Works for n>0.

h1_z_deriv <- function (n, p, ts) {

  ## Build the vector d that holds coefficients for all the derivatives of H1(z)
  ## The results reads d(n)*z^(1)*(d/dz)^(1)*H1(z) + d(n-1)*z^(2)*(d/dz)^(2)*H1(z) +...+ d(1)*z^(n)*(d/dz)^(n)*H1(z)
  d <- (-1)^n                                           # Vector with the derivatives of H1(z)
  for (i in 1:(n-1)) {
    d <- c(d, 0)                                        # Shift result right (multiply by -z)
    d <- d + prepad(pracma::polyder(d), i + 1, 0, 2)    # Add the derivative
  }

  ## Build output vector
  b <- rep(0L, n + 1)
  for (i in 1:n) {
    b  <- b + d[i] * prepad(h1_deriv(n-i+1), n+1, 0, 2)
  }

  b <- b * ts^(n + 1) / pracma::fact(n)

  ## Multiply coefficients with p^i, where i is the index of the coeff.
  b <- b * p^seq(n + 1, 1, -1)

  b
}

inv_residue <- function(r_in, p_in, k_in, tol) {

  n <- length(r_in) # Number of poles/residues

  k <- 0 # Capture constant term
  if (length(k_in) == 1) {         # A single direct term (order N = order D)
    k <- k_in[1]                   # Capture constant term
  } else if (length(k_in) > 1) {   # Greater than one means non-physical system
    stop("Order numerator > order denominator")
  }

  a_out <- poly(p_in)

  b_out  <- rep(0L, n + 1)
  b_out <- b_out + k * a_out       # Constant term: add k times denominator to numerator
  i <- 1
  while (i <= n) {
    term   <- c(1, -p_in[i])                          # Term to be factored out
    p      <- r_in[i] * pracma::deconv(a_out, term)$q # Residue times resulting polynomial
    p      <- prepad(p, n + 1, 0, 2)                  # Pad for proper length
    b_out <- b_out + p

    m          <- 1
    mterm      <- term
    first_pole <- p_in[i]
    while (i < n && abs(first_pole - p_in[i+1]) < tol){  # Multiple poles at p(i)
      i <- i + 1 # Next residue
      m <- m + 1
      mterm  <- conv(mterm, term)                        # Next multiplicity to be factored out
      p      <- r_in[i] * pracma::deconv(a_out, mterm)$q # Resulting polynomial
      p      <- prepad(p, n + 1, 0, 2)                   # Pad for proper length
      b_out  <- b_out + p
    }
    i <- i + 1
  }
  Arma(b_out, a_out)
}
