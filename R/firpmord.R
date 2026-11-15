# firpmord.R
# Copyright (C) 2026 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Octave signal package:
# Copyright (C) 2015 The Octave Project Developers
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; see the file COPYING. If not, see
# <https://www.gnu.org/licenses/>.
#
# 20261114 Geert van Boxtel          First version for v0.4.0
#------------------------------------------------------------------------------

#' Parks-McClellan optimal FIR filter order estimation
#'
#' Estimate the filter-order needed for \code{firpm} to design a type-I or
#' type-II linear-phase FIR filter according to the given specifications.
#'
#' The estimation algorithm used is per Ichige et al.[1]. Accuracy tends to
#' decrease as the number of bands increases. Even with two bands (i.e.
#' high-pass or low-pass), the algorithm may under- or over-estimate.
#' 
#' In order to precisely determine the minimum order needed for a particular
#' design, \code{firpmord} could be used to seed an algorithm iterating
#' invocations of \code{firpm}
#' 
#' @param f A vector of real-numbers, increasing in the range (0, fs/2), giving
#'   the frequencies of the left and right edges of each band for which a
#'   specific amplitude response is desired (omitting 0 and fs/2, which are
#'   implied): \code{c(r1, l2, r2,  …)}. Transition-bands are defined implicitly
#'   as the regions between the given bands.
#' @param a A vector of real-numbers giving the ideal amplitude response. An
#'   amplitude value is given for each band specified by f: \code{c(a1, a2, …)}.
#'   1 represents unity-gain, 0 represents infinite attenuation, and −1
#'   represents a phase change of pi radians.
#' @param d A vector of positive real-numbers giving the maximum allowable
#'   linear deviation from the amplitudes given in \code{a}, thus constraining
#'   the actual amplitude response (where defined by \code{f}) to be within
#'   \code{a +/− d}. Note that, though related, \code{d} does not equate to
#'   \code{firpm}’s \code{w} argument.
#' @param fs sampling frequency. Default: 2
#'
#' @return A list of class \code{firpmord} with the following elements:
#' \describe{
#'   \item{n}{Filter order}
#'   \item{f}{Numeric vector of normalized frequency-band edges.}
#'   \item{a}{amplitude response}
#'   \item{w}{weights}
#' }
#'
#' @examples
#' 
#' db2mag <- function(db) 10^(db / 20)
#' # Low-pass:
#' firpmord(f = c(0.4, 0.5), a = c(0, 1), d = c(db2mag(-72), db2mag(-72)))
#' firpmord(c(0.4, 0.5), c(0, 1), db2mag(-72))
#' 
#' # Multi-band-pass:
#' ds <- db2mag(-80); dp <- 0.01
#' firpmord(seq(1, 8) / 10, c(0, 1, 0, 1, 0), c(ds, dp, ds, dp, ds))
#' firpmord(seq(1, 8) / 10, c(0, 1, 0, 1, 0), c(ds, dp))
#'
#' @author The Octave Project Developers.\cr
#'  Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#'
#' @references [1] Ichige, K, Iwaki, M., and Ishii, R. (2000). Accurate
#'   estimation of minimum filter length for optimum FIR digital filters. IEEE
#'   Transactions on Circuits and Systems, 47(10), 1008-1017,
#'   \doi{10.1109/82.877143}
#' 
#' @seealso \code{\link{firpm}}, \code{\link{kaiserord}}
#'
#' @export

firpmord <- function(f, a, d, fs = 2) {
  
  
  if (!is.numeric(fs) || length(fs) != 1 || is.na(fs) || fs <= 0 ||
      is.complex(fs)) {
    stop("fs must be a positive real scalar")
  }
  
  f <- Re(as.vector(f)) / fs
  a <- Re(as.vector(a))
  d <- Re(as.vector(d))
  
  if (any(f <= 0) || any(f >= 0.5) ||
      any(diff(c(0, f)) <= 0)) {
    stop("f elements must increase between 0 and fs/2")
  } else if (length(f) == 0 || length(f) != 2 * length(a) - 2) {
    stop("f length must be twice the length of a, minus 2")
  } else if (length(d) < 1 || length(d) > length(a) || any(d <= 0)) {
    stop("d must be a vector of positive numbers with length(d) <= length(a)")
  }
  
  ## Repeat (as necessary) elements of d to obtain 1 element/band:
  l <- 1L
  if (length(d) < length(a)) {
      d <- rep_len(d, length(a))
  }
  
  ## Scale d w.r.t. a:
  d <- d / (abs(a) + (a == 0))
  
  ## Find the longest length needed to implement any of the
  ## low-pass or high-pass transitions within the given edges.
  ## HP transitions are mirrored in [0,.5] then treated as LP:
  
  fl <- f[seq(1, length(f) - 1, by = 2)]   # Freqs at left of transition.
  fr <- f[seq(2, length(f), by = 2)]       # Freqs at right of transition.
  L <- 2                                   # Anything less is not a filter.
  for (k in seq_along(fl)) {               # (fl, fr have equal length.)
    if (a[k] > a[k + 1]) {                 # Low-pass transition?
      L <- max(L, estimate_lp(fl[k], fr[k], d[k], d[k + 1]))
    } else if (a[k] < a[k + 1]) {          # High-pass transition?
      L <- max(L, estimate_lp(0.5 - fr[k], 0.5 - fl[k], d[k + 1], d[k]))
    }
  }
  
  ## Format length & specs for firpm:
  f <- c(0, f * 2, 1)
  w <- max(d) / d
  n <- ceiling(L) - 1                   # Filter order.
  if (a[length(a)] != 0 && n %% 2 == 1) {  # Fix-up high-pass/band-stop.
    n <- n + 1
  }
  a <- as.vector(kronecker(a, c(1, 1))) # kron(A, [1;1])
  
  result <- list(n = n, f = f, a = a, w = w)
  class(result) <- "firpmord"
  result
}

## Estimate filter-length needed for LP transition (Ichige, 2000):
estimate_lp <- function(fp, fs, dp, ds) {
  
  dF <- fs - fp
  stopifnot(dF > 0)
                                                                 # Ichige eqn.#
  v <- function(d) 2.325 * (-log10(d))^-0.445 * dF^-1.39         # (10)
  g <- function(fp, d) 2 / pi * atan(v(d) * (1 / fp - 1 / (0.5 - dF)))  # (9)
  h <- function(fp, c) 2 / pi * atan(c / dF * (1 / fp - 1 / (0.5 - dF))) # (15)
  
  Nc <- ceiling(1.101 * (-log10(2 * dp))^1.1 / dF + 1)           # (8) ...
  N3 <- ceiling(Nc * (g(fp, dp) + g(0.5 - fs, dp) + 1) / 3)      # (11)
  Nm <- 0.52 * log10(dp / ds) / dF * (-log10(dp))^0.17           # (15)
  DN <- ceiling(Nm * (h(fp, 1.1) - (h(0.5 - fs, 0.29) - 1) / 2)) # (15)
  N3 + DN                                                        # (15)
}
