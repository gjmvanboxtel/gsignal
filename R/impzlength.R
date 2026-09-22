# impzlength.R
# Copyright (C) 2026 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Original Octave function:
# Copyright (C) 2026 Tang Chonghao <chadholton@qq.com>
#s
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
# 20260922  GvB       setup for gsignal v0.4.0
#------------------------------------------------------------------------------

#' Impulse response length
#' 
#'  For a finite impulse response (FIR) filter specified by the numerator
#'  coefficients \code{b}, the length is simply the number of coefficients in
#'  \code{b}.
#'
#' For an infinite impulse response (IIR) filter specified by the numerator
#' \code{b} and denominator \code{a} polynomials in z^-1, the function computes
#' an effective impulse response sequence length.
#'
#' @details
#' The algorithm proceeds as follows:
#' \itemize{
#'    \item{If the filter is FIR, the length is simply length (b)}.
#'    \item{The poles of the transfer function are computed as the roots
#'     of the denominator polynomial \code{a}.}
#'    \item{The multiplicity of the dominant pole (the pole with the largest
#'    magnitude) is determined by counting poles at the same complex coordinate
#'    within tolerance.}
#'    \item{For a stable IIR filter (dominant pole magnitude < 1 - 10^{-5}),
#'     the effective length is estimated as
#'     
#'     \code{floor (M * log10 (tol) / log10 (maxpole)) + delay}
#'     
#'     where M is the multiplicity of the dominant pole and d is the initial
#'      delay (number of leading zeros in b).}
#'    \item{For an unstable IIR filter (dominant pole magnitude > 1 + 10^{-4}),
#'     a heuristic formula is used:
#'     
#'     \code{floor (6 / log10 (maxpole))}}
#'    \item{For filters with poles near the unit circle (oscillatory behavior),
#'     the length is the maximum of: five periods of the slowest oscillation,
#'      and the decay length of damped poles, plus the initial delay. }
#' }
#'
#' @param filt For the default case, the moving-average coefficients of an ARMA
#'   filter (normally called \code{b}), specified as a numeric or complex
#'   vector. Generically, \code{filt} specifies an arbitrary model or filter
#'   operation.
#' @param a the autoregressive (recursive) coefficients of an ARMA filter,
#'   specified as a numeric or complex vector. If \code{a[1]} is not equal to 1,
#'   then filter normalizes the filter coefficients by \code{a[1]}. Therefore,
#'   \code{a[1]} must be nonzero.
#' @param tol specifies the tolerance used to estimate the effective length of
#'   an IIR filter's impulse response. The default tolerance is 5e-5.
#'   Increasing @var{tol} estimates a shorter effective length, while decreasing
#'   @var{tol} produces a longer effective length.
#' @param ... additional arguments (ignored).
#'
#' @return Length of the impulse response, specified as a positive integer. For
#'   stable IIR filters with absolutely summable impulse responses, impzlength
#'   returns an effective length for the impulse response beyond which the
#'   coefficients are essentially zero. You can control this cutoff point by
#'   specifying the optional tol input argument.
#'
#' @examples
#' b <- 1
#' a <- c(1, -0.9)
#' len <- impzlength(b, a)
#'
#' @author Tang Chonghao, \email{chadholton@@qq.com}.\cr Conversion to R by
#'   Geert van Boxtel, \email{gjmvanboxtel@@gmail.com}
#'
#' @rdname impzlength
#' @export

impzlength <- function(filt, ...) UseMethod("impzlength")

#' @rdname impzlength
#' @export

impzlength.Arma <- function(filt, ...)
  impzlength(filt$b, filt$a, ...)

#' @rdname impzlength
#' @export

impzlength.Ma <- function(filt, ...)
  impzlength(filt, a = 1, ...)

#' @rdname impzlength
#' @export

impzlength.Sos <- function(filt, ...)
  impzlength(as.Arma(filt), ...)

#' @rdname impzlength
#' @export

impzlength.Zpg <- function(filt, ...)
  impzlength(as.Arma(filt), ...)


#' @rdname impzlength
#' @export

impzlength.default <- function(filt, a = 1, tol = 5e-5, ...)  {
  
  if (!(is.numeric(filt) || is.complex(filt)) ||
      !(is.numeric(a) || is.complex(a))) {
    stop("a and b must be numeric or complex")
  }
  
  if(!isPosscal(tol)) {
    stop("tol must be a positive scalar")
  }
  
  b <- filt
  
  ## FIR filter: length equals number of numerator coefficients
  if (isScalar(a) || is.null(a)) {
    len <- length(b)
    return(len)
  }
  
  ##  find delay index of first non-zero coefficient in b
  idx <- which(b != 0)[1]
  if (is.na(idx) || length(idx) <= 0) {
    b_delay <- 0
  } else {
    b_delay <- idx - 1
  }
  
  # Find poles of the transfer function
  r <- pracma::roots(a)
  abs_r <- abs(r)
  maxpole <- max(abs_r)
  
  # Find multiplicity of the dominant pole (true repeated roots).
  # Count poles at the same complex coordinate, NOT just same magnitude,
  # to avoid counting complex conjugate pairs as repeated.
  idx_maxmag <- as.numeric(abs(abs_r - maxpole) < tol)
  candidates <- r[idx_maxmag]
  if (length(candidates) > 0) {
    ref <- candidates[1]
    mult <- sum(as.numeric(abs(r - ref) < tol))
  } else {
    mult <- 1
  }
  
  if (maxpole > 1 + 1e-4) {
    ## Unstable IIR filter
    len <- floor (6 / log10(maxpole))
  } else if (maxpole < 1 - 1e-5) {
    ## Stable IIR filter
    len <- floor(mult * log10(tol) / log10(maxpole)) + b_delay
  } else {
    ## Filter has poles near the unit circle (oscillatory).
    ## Compute five periods of the slowest oscillation.
    n_periodic <- 10
    for (i in seq_along(r)) {
      if (abs(r[i] - 1) < 1e-5) {
        r[i] <- -r[i]  ## avoid numerical issues with arg(r) near 0
      }
    }
    rperiodic <- r[which(abs_r >= 1 - 1e-5 & abs_r <= 1 + 1e-4 
                         & abs (Arg(r)) > 0)]
    if (length(rperiodic) > 0) {
      n_periodic <- floor(5 * 2 * pi / min(abs(Arg(rperiodic))))
    }
  
    ## Compute damped component length
    n_damped <- 0
    rdamped = r[which(abs_r < 1 - 1e-5)]
    if (length(rdamped) > 0) {
      max_rdamped <- max(abs(rdamped))
  
      ## Find true repeated roots among damped poles (same complex coordinate)
      idx_damped_max <- which((abs(abs(rdamped) - max_rdamped) < tol))
      candidates_damped <- rdamped[idx_damped_max]
  
      if (length(candidates_damped) > 0) {
        ref_damped <- candidates_damped[1]
        mult_damped <- sum(which(abs(rdamped - ref_damped) < tol))
      } else {
        mult_damped <- 1
      }
  
      n_damped <- floor(mult_damped * log10(tol) / log10 (max_rdamped))
    }
  
    len <- max (n_periodic, n_damped) + b_delay
  }
  
  len

}
