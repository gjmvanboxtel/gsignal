# welchwin.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Octave signal package:
# Copyright (C) 2007 Muthiah Annamalai <muthiah.annamalai@uta.edu>
# Copyright (C) 2008-2009 Mike Gross <mike@appl-tech.com>
# Copyright (C) 2008-2009 Peter V. Lanspeary <pvl@mecheng.adelaide.edu.au>
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
# 20191210 Geert van Boxtel          First version for v0.1.0
#------------------------------------------------------------------------------

#' Welch window
#'
#' Return the filter coefficients of a Welch window of length \code{n}.
#'
#' The Welch window is a polynomial window consisting of a single parabolic
#' section:
#' \deqn{w(k) = 1 - (k / N - 1)^2, n=0,1, ... n-1}.
#' The optional argument specifies a "symmetric" window (the default) or a
#' "periodic" window. A symmetric window has zero at each end and maximum in the
#' middle, and the length must be an integer greater than 2. The variable
#' \code{N} in the formula above is \code{(n-1)/2}. A periodic window wraps
#' around the cyclic interval \code{0,1, ... m-1}, and is intended for use
#' with the DFT. The length must be an integer greater than 1. The variable
#' \code{N} in the formula above is \code{n/2}.
#'
#' @param n Window length, specified as a positive integer.
#' @param method Character string. Window sampling method, specified as:
#' \describe{
#'   \item{"symmetric"}{(Default). Use this option when using windows for filter
#'   design.}
#'   \item{"periodic"}{This option is useful for spectral analysis because it
#'   enables a windowed signal to have the perfect periodic extension implicit
#'   in the discrete Fourier transform. When 'periodic' is specified, the
#'   function computes a window of length \code{n + 1} and returns the first
#'   \code{n} points.}
#' }
#'
#' @return Welch window, returned as a vector.
#'
#' @examples
#'
#' w <- welchwin(64)
#' plot (w, type = "l", xlab = "Samples", ylab =" Amplitude")
#'
#' ws = welchwin(64,'symmetric')
#' wp = welchwin(63,'periodic')
#' plot (ws, type = "l", xlab = "Samples", ylab =" Amplitude")
#' lines(wp, col="red")
#'
#' @author Muthiah Annamalai, \email{muthiah.annamalai@@uta.edu},\cr
#' Mike Gross, \email{mike@@appl-tech.com},\cr
#' Peter V. Lanspeary, \email{pvl@@mecheng.adelaide.edu.au}.\cr
#' Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#
#' @export

welchwin <- function(n, method = c("symmetric", "periodic")) {

  if (!isPosscal(n) || ! isWhole(n) || n <= 0)
    stop("n must be an integer strictly positive")
  method <- match.arg(method)

  if (method == "periodic") {
    N <- n / 2
    nmin <- 2
  } else if (method == "symmetric") {
    N <- (n - 1) / 2
    nmin <- 3
  } else {
    stop("method must be either 'periodic' or 'symmetric'")
  }

  ## Periodic window is not properly defined for m < 2.
  ## Symmetric window is not properly defined for m < 3.
  if (n < nmin) {
    stop(paste("n must be an integer greater than", nmin))
  }

  k <- 0:(n - 1)
  w <- 1 - ((k - N) / N)^2
  w
}
