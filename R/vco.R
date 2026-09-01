# vco.R
# Copyright (C) 2026 Geert van Boxtel <gjmvanboxtel@gmail.com>
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
# 20260901 Geert van Boxtel          First version for v0.4.0
#------------------------------------------------------------------------------

#' Voltage-controlled oscillator
#'
#' Creates a signal that oscillates at a frequency determined by input \code{x}
#' with a sampling frequency fs. 
#'
#' @param x Input data, specified as a real vector or real matrix. \code{x}
#'   ranges from –1 to 1, where \code{x = –1} corresponds to 0 frequency output,
#'   \code{x = 0} corresponds to \code{fc}, and \code{x = 1} corresponds to
#'   \code{2 * fc}. If \code{x} is a matrix, \code{vco} produces a matrix whose
#'   columns oscillate according to the columns of \code{x}.
#' 
#' @param fc Carrier frequency used to modulate the input signal, specified as a
#'   real positive scalar, or a vector of 2 positive scalars. In the latter
#'   case, \code{fc} give the frequency modulation range limits. Values must be
#'   in the range \code{0} to \code{fs / 2}. Default: \code{fs / 4}.
#'  
#' @param fs Sampling rate, specified as a positive scalar. Default: 1
#'
#' @return Output signal, returned as a real vector or real matrix. It is the
#'   same size as \code{x} and has amplitude equal to 1.
#'
#' @seealso \code{\link{pulstran}}
#'
#' @examples
#'
#' x <- c(0, 0.5, 0.8, -0.2)
#' y1 <- vco (x, 10, 100)
#' y2 <- vco (x, c(5, 15), 100)
#' 
#' fs <- 10000
#' t <- seq(0, 2, 1 / fs)
#' x <- sawtooth(2 * pi * t, 0.75)
#' y <- vco(x, c(0.1, 0.4) * fs, fs)
#' specgram(y, fs = fs)
#' 
#' @author The Octave Project Developers, \email{maintainers@@octave.org}, 
#'   Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#
#' @export

vco <- function(x, fc = fs / 4, fs = 1) {

  if (!is.numeric(x)) {
    stop("x must be numeric")
  }
  
  if (is.vector(x)) {
    x <- matrix(x, ncol = 1)
    vec <- TRUE
  } else if (is.matrix(x)) {
    vec <- FALSE
  } else {
    stop("x must be a numeric vector or matrix")
  }
  
  if (max(abs(x)) > 1) {
    stop("x must be in the range between -1 and 1")
  }
  
  if (!isPosscal(fs)) {
    stop("fs must be a positive scalar")
  }
  
  if (length(fc) == 2) {
    fmin <- fc[1]
    fmax <- fc[2]
    if (fmin > fs / 2 || fmax > fs / 2) {
      stop ("the carrier frequency fc must be less than fs / 2")
    }
    if (fmin >= fmax) {
      stop ("fc[2] must be greater than fc[1]")
    }
    if (!(isPosscal(fmin) && isPosscal(fmax))) {
      stop ("fc must contain positive scalars");
    }
    # Comment in Octave function:
    ## I don't know whether this is the best way, but I think it gives
    ## the proper value for fc and the proper deviation (+/-) from
    ## that value.
    fc <- (fmin + fmax) / 2
    kf <- (fmax - fmin) / 2 / fs * 2 * pi
  } else {
    if (length(fc) == 1 && isPosscal(fc)) {
      kf <- (fc / fs) * 2 * pi
    } else {
      stop ("fc must be a scalar or a vector of two elements")
    }
  }
  
  # Comment in Octave function:
  ## Note, no need to generate a matrix for T, we can rely on
  ## broadcasting instead.
  len <- nrow (x)
  t <- seq(0, (len - 1) / fs, 1 / fs)
  
  y <- apply(x, 2, function(x) cos (2 * pi * fc * t + kf * cumsum (x)))

  if (vec) {
    y <- as.vector(y)
  }
  dimnames(y) <- dimnames(x)
  y
}
