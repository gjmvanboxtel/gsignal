# firls.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Octave signal package:
# Copyright (C) 2006 Quentin Spencer <qspencer@ieee.org>
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
# 20200706 Geert van Boxtel           First version for v0.1.0
#------------------------------------------------------------------------------

#' Least-squares linear-phase FIR filter design
#'
#' Produce a linear phase filter such that the integral of the weighted mean
#' squared error in the specified bands is minimized.
#'
#' @param n filter order (1 less than the length of the filter). Must be even.
#'   If odd, it is incremented by one.
#' @param f vector of frequency points in the range from 0 to 1, where 1
#'   corresponds to the Nyquist frequency. Each band is specified by two
#'   frequencies, so the vector must have an even length. .
#' @param a vector of the same length as \code{f} containing the desired
#'   amplitude at each of the points specified in \code{f}.
#' @param w weighting function that contains one value for each band that
#'   weights the mean squared error in that band. \code{w} must be half the
#'   length of \code{f}.
#'
#' @return The FIR filter coefficients, a vector of length \code{n + 1}, of
#'   class \code{Ma}.
#'
#' @examples
#' freqz(firls(255, c(0, 0.25, 0.3, 1), c(1, 1, 0, 0)))
#'
#' @seealso \code{\link{Ma}}, \code{\link{filter}}, \code{\link{fftfilt}},
#'   \code{\link{fir1}}
#'
#' @author Quentin Spencer, \email{qspencer@@ieee.org}.\cr
#'   Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#'
#' @export

firls <- function(n, f, a, w = rep(1L, length(a) / 2)) {

  # filter length must be a scalar > 0
  if (!isPosscal(n) || !isWhole(n) || n <= 0) {
    stop("n must be an integer > 0")
  }

  # f, a, w must be real-valued vectors
  if (!is.vector(f) || !is.numeric(f) || !is.vector(a) || !is.numeric(a) ||
     !is.vector(w) || !is.numeric(w)) {
    stop("f, a, and w must be real-valued vectors")
  }

  # test for lengths of f, a, and w
  if (length(f) != length(a)) {
    stop("f and a must have equal lengths")
  } else if (2 * length(w) != length(a)) {
    stop("w must contain one weight per band")
  }

  n <- n + n %% 2
  M <- n / 2

  ww <- kronecker(w, c(-1, 1))
  omega <- f * pi
  i1 <- seq(1, length(omega), 2)
  i2 <- seq(2, length(omega), 2)

  ## Generate the matrix Q (see [1])
  cos_ints <- rbind(omega, sin((1:n) %o% omega))
  q <- c(1, 1 / (1:n)) * (cos_ints %*% ww)
  Q <- stats::toeplitz(q[1:(M + 1)]) +
    pracma::hankel(q[1:(M + 1)], q[(M + 1):length(q)])

  ## Derive the vector b (see [1])
  cos_ints2 <- rbind(omega[i1]^2 - omega[i2]^2,
                     cos((1:M) %o% omega[i2]) - cos((1:M) %o% omega[i1])) /
    (c(2, 1:M) %o% (omega[i2] - omega[i1]))
  d <- as.vector(rbind(-w * a[i1], w * a[i2]))
  b <- c(1, 1 / (1:M)) * ((kronecker(cos_ints2, cbind(1, 1)) +
                             cos_ints[1:(M + 1), ]) %*% d)

  ## Having computed the components Q and b of the  matrix equation,
  ## solve for the filter coefficients.
  aa <- pracma::mldivide(Q, b, pinv = FALSE)
  laa <- length(aa)
  coef <- c(aa[seq(laa, 2, -1)], 2 * aa[1], aa[2:laa])
  Ma(coef)
}
