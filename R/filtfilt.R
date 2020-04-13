# filtfilt.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Octave version:
# Copyright (C) 1999 Paul Kienzle <pkienzle@users.sf.net>
# Copyright (C) 2007 Francesco Potortì <pot@gnu.org>
# Copyright (C) 2008 Luca Citi <lciti@essex.ac.uk>
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
# 20200217  GvB       setup for gsignal v0.1
# 20200413  GvB       added S3 method method for Sos
#---------------------------------------------------------------------------------------------------------------------

#' Zero-phase digital filtering
#' 
#' Forward and reverse filter the signal.
#' 
#' Forward and reverse filtering the signal corrects for phase distortion
#' introduced by a one-pass filter, though it does square the magnitude response
#' in the process. That’s the theory at least. In practice the phase correction
#' is not perfect, and magnitude response is distorted, particularly in the stop
#' band.
#' 
#' In Matlab filtfilt reduces filter startup transients by carefully choosing
#' initial conditions, and by prepending onto the input sequence a short,
#' reflected piece of the input sequence. The current (1.4.1) Octave version uses a
#' slightly different method to choose initial conditions. Neither of these
#' methods have been implemented in the current version (mainly because I did
#' not entirely understand them). Here, a reflected sequence of the input signal
#' is added to the beginning and end of the signal, and tapered to zero, as per
#' the recommendations on the Matlab website. This is different from the current
#' (0.7-6) R signal package, which pads the input signal with zeroes.
#' 
#' @param filt For the default case, the moving-average coefficients of an ARMA
#'   filter (normally called ‘b’). Generically, \code{filt} specifies an arbitrary
#'   filter operation.
#' @param a the autoregressive (recursive) coefficients of an ARMA filter.
#' @param x the input signal to be filtered. If \code{x} is a matrix, all
#' coulums
#' @param ... additional arguments (ignored).
#' 
#' @return The filtered signal, normally of the same length of the input signal
#'   \code{x}, returned as a vector or matrix
#' 
#' @examples
#' #bf <- butter(3, 0.1)                                # 10 Hz low-pass filter
#' bf <- Arma(c(0.002898195, 0.008694584, 0.008694584, 0.002898195),
#'  c(1.0000000, -2.3740947,  1.9293557, -0.5320754))   # change later
#' t <- seq(0, 1, len = 100)                            # 1 second sample
#' x <- sin(2* pi * t * 2.3) + 0.25 * rnorm(length(t))  # 2.3 Hz sinusoid+noise
#' z <- filter(bf, x)                                   # apply filter
#' plot(t, x, type = "l")
#' lines(t, z, col = "red")
#' zz <- filtfilt(bf, x)
#' lines(t, zz, col="blue")
#' legend("bottomleft", legend = c("original", "filter", "filtfilt"), lty = 1,
#'  col = c("black", "red", "blue"))
#' 
#' @seealso \code{\link{filter}}
#' 
#' @author Paul Kienzle \email{pkienzle@@users.sf.net}, Francesco Potortì
#'   \email{pot@@gnu.org}, Luca Citi \email{lciti@@essex.ac.uk}, port to R by
#'   Geert van Boxtel \email{G.J.M.vanBoxtel@@gmail.com}.
#' 
#' @rdname filtfilt
#' @export

filtfilt <- function(filt, ...) UseMethod("filtfilt")

#' @rdname filtfilt
#' @method filtfilt default
#' @export

filtfilt.default <- function(filt, a, x, ...) {
  
  ff_single <- function(x, filt, a) {
    if (missing(a)) {
      fl <- length(filt)
    } else {
      fl <- max(length(filt), length(a))
    }
    if (fl >= length(x)) {
      stop("filter length must be shorter than series length")
    }
    hann <- hann(2 * fl + 1)
    x <- c((hann[1:fl] * rev(x[1:fl])), x, (hann[(fl + 2):length(hann)] * rev(x[(length(x) - fl + 1):length(x)])))
    y <- filter(filt, a, x)
    y <- rev(filter(filt, a, rev(y)))
    y[(fl + 1):(length(y) - fl)]
  }
  
  d <- dim(x)
  if (is.null(d) || (length(d) == 2 && d[2] == 1)) {
    y <- ff_single(x, filt, a)
  } else if (length(d) == 2 && d[2] > 1) {
    y <- apply(x, 2, ff_single, filt = filt, a = a)
  } else {
    stop('Incorrect array dimensions')
  }
  y
}

#' @rdname filtfilt
#' @method filtfilt Arma
#' @export
filtfilt.Arma <- function(filt, x, ...) # IIR
  filtfilt(filt$b, filt$a, x, ...)

#' @rdname filtfilt
#' @method filtfilt Ma
#' @export
filtfilt.Ma <- function(filt, x, ...) # FIR
  filtfilt(unclass(filt), 1, x, ...)

#' @rdname filtfilt
#' @method filtfilt Sos
#' @export
filtfilt.Sos <- function(filt, x, ...) { # Second-order sections

  ff_single <- function(x, sos) {
    fl <- 3
    hann <- hann(2 * fl + 1)
    x <- c((hann[1:fl] * rev(x[1:fl])), x, (hann[(fl + 2):length(hann)] * rev(x[(length(x) - fl + 1):length(x)])))
    y <- sosfilt(sos, x)
    y <- rev(sosfilt(sos, rev(y)))
    y[(fl + 1):(length(y) - fl)]
  }
  
  if (filt$g != 1) {
    filt$sos[1, 1:3] <- filt$sos[1, 1:3] * filt$g
  }
  
  d <- dim(x)
  if (is.null(d) || (length(d) == 2 && d[2] == 1)) {
    y <- ff_single(x, filt$sos)
  } else if (length(d) == 2 && d[2] > 1) {
    y <- apply(x, 2, ff_single, sos = filt$sos)
  } else {
    stop('Incorrect array dimensions')
  }
  y
}

#' @rdname filtfilt
#' @method filtfilt Zpg
#' @export
filtfilt.Zpg <- function(filt, x, ...) # zero-pole-gain form
  filtfilt(as.Arma(filt), x)
