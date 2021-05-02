# tukeywin.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Octave signal package:
# Copyright (C) 2007 Laurent Mazet <mazet@crm.mot.com>
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
# 20200101 Geert van Boxtel          First version for v0.1.0
#------------------------------------------------------------------------------

#' Tukey (tapered cosine) window
#'
#' Return the filter coefficients of a Tukey window (also known as the
#' cosine-tapered window) of length \code{n}.
#'
#' The Tukey window, also known as the tapered cosine window, can be regarded as
#' a cosine lobe that is convolved with a rectangular window.  \code{r} defines
#' the ratio between the constant section and and the cosine section. It has to
#' be between 0 and 1. The function returns a Hann window for \code{r} equal to
#' 1 and a rectangular window for \code{r} equal to 0.
#'
#' @param n Window length, specified as a positive integer.
#' @param r Cosine fraction, specified as a real scalar. The Tukey window is a
#'   rectangular window with the first and last \code{r / 2} percent of the
#'   samples equal to parts of a cosine. For example, setting \code{r = 0.5}
#'   (default) produces a Tukey window where 1/2 of the entire window length
#'   consists of segments of a phase-shifted cosine with period 2r = 1. If you
#'   specify r <= 0, an n-point rectangular window is returned. If you specify r
#'   >= 1, an n-point von Hann window is returned.
#'
#' @return Tukey window, returned as a vector.
#'
#' @examples
#'
#' n <- 128
#' t0 <- tukeywin(n, 0)        # Equivalent to a rectangular window
#' t25 <- tukeywin(n, 0.25)
#' t5 <- tukeywin(n)           # default r = 0.5
#' t75 <- tukeywin(n, 0.75)
#' t1 <- tukeywin(n, 1)         # Equivalent to a Hann window
#' plot(t0, type = "l", xlab = "Samples", ylab =" Amplitude", ylim=c(0,1.2))
#' lines(t25, col = 2)
#' lines(t5, col = 3)
#' lines(t75, col = 4)
#' lines(t1, col = 5)
#'
#' @author Laurent Mazet, \email{mazet@@crm.mot.com}.\cr
#' Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#
#' @export

tukeywin <- function(n, r = 1 / 2) {

  if (!isPosscal(n) || ! isWhole(n) || n <= 0)
    stop("n must be an integer strictly positive")
  if (r > 1) {
    r <- 1
  } else if (r < 0) {
    r <- 0
  }

  if (r == 0) {
    w <- rep(1L, n)
  } else if (r == 1) {
    w <- hann(n)
  } else {
    if (n == 1) {
      w <- 1
    } else {
      # cosine-tapered window
      k <- seq(0, 1, length.out = n)[1:(n / 2)]
      w <- (1 + cos(pi * (2 * k / r - 1))) / 2
      idx <- (floor(r * (n - 1) / 2) + 2):length(w)
      if (!any(idx > length(w))) {
        w[idx] <- 1
      }
      w <- c(w, rep(1, n %% 2),  rev(w))
    }
  }
  w
}
