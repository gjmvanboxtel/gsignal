# pei_tseng_notch.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Octave signal package:
# Copyright (C) 2011 Alexander Klein <alexander.klein@math.uni-giessen.de>
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
# 20200622 Geert van Boxtel           First version for v0.1.0
#------------------------------------------------------------------------------

#' Pei-Tseng notch filter
#'
#' Compute the transfer function coefficients of an IIR narrow-band notch
#' filter.
#'
#' The filter construction is based on an all-pass which performs a reversal of
#' phase at the filter frequencies. Thus, the mean of the phase-distorted and
#' the original signal has the respective frequencies removed.
#'
#' @param w vector of critical frequencies of the filter. Must be between 0
#'   and 1 where 1 is the Nyquist frequency.
#' @param bw vector of bandwidths. Bw should be of the same length as \code{w}.
#'
#' @return List of class \code{\link{Arma}} with list elements:
#' \describe{
#'   \item{b}{moving average (MA) polynomial coefficients}
#'   \item{a}{autoregressive (AR) polynomial coefficients}
#' }
#'
#' @examples
#' ## 50 Hz notch filter
#' fs <- 256
#' nyq <- fs / 2
#' notch <- pei_tseng_notch(50 / nyq, 2 / nyq)
#' freqz(notch, fs = fs)
#'
#' @references Pei, Soo-Chang, and Tseng, Chien-Cheng "IIR Multiple Notch Filter
#'   Design Based on Allpass Filter"; 1996 IEEE Tencon, doi:
#'   \doi{10.1109/TENCON.1996.608814}
#'
#' @seealso \code{\link{Arma}}, \code{\link{filter}}
#'
#' @author Alexander Klein, \email{alexander.klein@@math.uni-giessen.de}.\cr
#'  Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#'
#' @export

pei_tseng_notch <- function(w, bw) {

  # check input arguments
  if (!is.vector(w) || !is.vector(bw)) {
    stop("All arguments must be vectors")
  }
  if (length(w) != length(bw)) {
    stop("All arguments must be of equal length")
  }
  if (!all(w > 0 && bw < 1)) {
    stop("All frequencies must be in the range (0, 1)")
  }
  if (!all(bw > 0 && bw < 1)) {
    stop("All bandwidths must be in the range (0, 1)")
  }

  ## Normalize appropriately
  w  <- w * pi
  bw <- bw * pi
  M2 <- 2 * length(w)

  ## Splice center and offset frequencies (Equation 11)
  omega <- as.vector(rbind(w - bw / 2, w))

  ## Splice center and offset phases (Equations 12)
  factors <- seq(1, M2, 2)
  phi     <- as.vector(rbind(-pi * factors + pi / 2, -pi * factors))

  ## Create linear equation
  t_beta <- tan((phi + M2 * omega) / 2)

  Q <- matrix(0L, nrow = M2, ncol = M2)

  for (k in seq_len(M2)) {
    Q [, k] <- sin(k * omega) - t_beta * cos(k * omega)
  }

  ## Compute coefficients of system function (Equations 19, 20) ...
  h_a   <- as.vector(pracma::mldivide(Q, t_beta))
  denom <- c(1, h_a)
  num   <- c(rev(h_a), 1)

  ## ... and transform them to coefficients for difference equations
  a <- denom
  b <- (num + denom) / 2
  Arma(b, a)
}
