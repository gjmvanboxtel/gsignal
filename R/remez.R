# remez.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Octave signal package:
# Copyright (C) 1995, 1998 Jake Janovetz <janovetz@uiuc.edu>
# Copyright (C) 1999 Paul Kienzle <pkienzle@users.sf.net>
# Copyright (C) 2000 Kai Habel <kahacjde@linux.zrz.tu-berlin.de>
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
# 20200803 Geert van Boxtel           First version for v0.1.0
#------------------------------------------------------------------------------

#' Parks-McClellan optimal FIR filter design
#'
#' Parks-McClellan optimal FIR filter design using the Remez exchange algorithm.
#'
#' @param n filter order (1 less than the length of the filter).
#' @param f normalized frequency points, strictly increasing vector in the range
#'   [0, 1], where 1 is the Nyquist frequency. The number of elements in the
#'   vector is always a multiple of 2.
#' @param a vector of desired amplitudes at the points specified in \code{f}.
#'   \code{f} and \code{a} must be the same length. The length must be an even
#'   number.
#' @param w vector of weights used to adjust the fit in each frequency band. The
#'   length of \code{w} is half the length of \code{f} and \code{a}, so there is
#'   exactly one weight per band. Default: 1.
#' @param ftype filter type, matched to one of \code{"bandpass"} (default),
#'   \code{"differentiatior"}, or \code{"hilbert"}.
#' @param density determines how accurately the filter will be constructed. The
#'   minimum value is 16 (default), but higher numbers are slower to compute.
#'
#' @return The FIR filter coefficients, a vector of length \code{n + 1}, of
#'   class \code{Ma}
#'
#' @references \url{https://en.wikipedia.org/wiki/Fir_filter}
#'
#' @examples
#' ## low pass filter
#' f1 <- remez(15, c(0, 0.3, 0.4, 1), c(1, 1, 0, 0))
#' freqz(f1)
#'
#' ## band pass
#' f <- c(0, 0.3, 0.4, 0.6, 0.7, 1)
#' a <- c(0, 0, 1, 1, 0, 0)
#' b <- remez(17, f, a)
#' hw <- freqz(b, 512)
#' plot(f, a, type = "l", xlab = "Radian Frequency (w / pi)",
#'      ylab = "Magnitude")
#' lines(hw$w/pi, abs(hw$h), col = "red")
#' legend("topright", legend = c("Ideal", "Remez"), lty = 1,
#'        col = c("black", "red"))
#'
#' @seealso \code{\link{Ma}}, \code{\link{filter}}, \code{\link{fftfilt}},
#'   \code{\link{fir1}}
#'
#' @author Jake Janovetz, \email{janovetz@@uiuc.edu},\cr
#'   Paul Kienzle, \email{pkienzle@@users.sf.net},\cr
#'   Kai Habel, \email{kahacjde@@linux.zrz.tu-berlin.de}.\cr
#'   Conversion to R Tom Short\cr
#'   adapted by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#'
#' @references Rabiner, L.R., McClellan, J.H., and Parks, T.W. (1975). FIR
#'   Digital Filter Design Techniques Using Weighted Chebyshev Approximations,
#'   IEEE Proceedings, vol. 63, pp. 595 - 610.\cr
#'   \url{https://en.wikipedia.org/wiki/Parks-McClellan_filter_design_algorithm}
#'
#' @export

remez <- function(n, f, a, w = rep(1.0, length(f) / 2),
                  ftype = c("bandpass", "differentiator", "hilbert"),
                  density = 16) {

  ftype <- as.integer(factor(match.arg(ftype),
                             c("bandpass", "differentiator", "hilbert")))

  if (!isPosscal(n) || !isWhole(n) || n < 4) {
    stop("Filter length n must be an integer greater than 3")
  }
  if (!is.vector(f) || length(f) %% 2 == 1) {
    stop("f must be a vector of even length")
  }
  if (any(diff(f) < 0)) {
    stop("f must be a vector of increasing numbers")
  }
  if (any(f < 0) || any(f > 1)) {
    stop("f must be in the range [0,1]")
  }
  if (length(a) != length(f)) {
    stop("length(a) must equal length(f)")
  }
  if (2 * length(w) != length(f)) {
    stop("length(w) must be half of length(f)")
  }
  if (density < 16) {
    stop("density is too low, must be greater than or equal to 16")
  }

  z <- .Call("_gsignal_remez",
          h = as.double(rep(0, n + 1)),
          as.integer(n + 1),
          as.integer(length(f) / 2),
          as.double(f / 2),
          as.double(a),
          as.double(w),
          as.integer(ftype),
          as.integer(density),
          PACKAGE = "gsignal")

  Ma(z)
}
