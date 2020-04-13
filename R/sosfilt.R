# sosfilt.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
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
# 20200413  GvB       setup for gsignal v0.1.0
#---------------------------------------------------------------------------------------------------------------------

#' Second-order filtering
#' 
#' One-dimensional second-order (biquadratic) IIR digital filtering
#' 
#' @param sos Second-order section representation, specified as an nrow-by-6
#'   matrix, whose rows contain the numerator and denominator coefficients of
#'   the second-order sections:\cr \code{sos <- rbind(cbind(B1, A1), cbind(...),
#'   cbind(Bn, An))}, where \code{B1 <- c(b0, b1, b2)}, and \code{A1 <- c(a0,
#'   a1, a2)} for section 1, etc. The b0 entry must be nonzero for each section.
#' @param x The input data to be filtered, coerced to a vector.
#' 
#' @return The filtered signal, normally of the same length as the input signal
#'   \code{x}, returned as a vector
#' 
#' @examples
#' fs <- 1000                                           # sampling frequency
#' t <- seq(0, 1, 1/fs)                                 # 1 second sample
#' s <- sin(2* pi * t * 6)                              # 6 Hz sinus
#' x <- s + rnorm(length(t))                            # add noise
#' plot(t, x, type = "l", col="light gray")
#' lines(t, s, col="black")
#' #bf <- signal::butter(3, 0.02)                        # low-pass 0.02 * 1000 = 20 Hz
#' bf <- Arma(c(2.914649e-05, 8.743948e-05, 8.743948e-05, 2.914649e-05),
#'            c(1.0000000, -2.8743569,  2.7564832, -0.8818931))
#' sosg <- as.Sos(bf)                                   # convert to second order sections
#' sos <- sosg$sos
#' sos[1, 1:3] <- sos[1, 1:3] * sosg$g                  # apply gain factor
#' y <- sosfilt(matrix(sos, ncol=6), x)                 # apply filter
#' lines(t, y, col="red")
#' yy <- filtfilt(sosg, x)                              # filtfilt and filter handle gain factor
#' lines(t, yy, col="blue")
#' 
#' @seealso \code{\link{filter}}, \code{\link{filtfilt}}
#' 
#' @author Geert van Boxtel \email{G.J.M.vanBoxtel@@gmail.com}.
#' 
#' @export

sosfilt <- function(sos, x) {
  
  if (is.vector(sos)) {
    if (length(sos) == 6) {
      sos <- matrix(sos, ncol = 6)
    } else {
      stop('sos must a matrix with 6 columns')
    }
  } else if (is.matrix(sos)) {
    if(ncol(sos) != 6) {
      stop('sos must a matrix with 6 columns')
    }
  } else {
    stop('sos must a matrix with 6 columns')
  }

  x <- as.vector(x)

  if (any(sos[, 4] == 0)) {
    y <- rep(NA, length(x))
  } else {
    y <- .Call("_gsignal_sosfilt", PACKAGE = "gsignal", sos, x)
  }
  
  y
}