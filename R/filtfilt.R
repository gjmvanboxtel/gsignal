# filtfilt.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Octave version:
# Copyright (C) 1999 Paul Kienzle <pkienzle@users.sf.net>
# Copyright (C) 2007 Francesco Potortì <pot@gnu.org>
# Copyright (C) 2008 Luca Citi <lciti@essex.ac.uk>
#
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
# 20200217  GvB       setup for gsignal v0.1
# 20200413  GvB       added S3 method method for Sos
# 20210402  GvB       use padding and Gustafsson method for initial conditions
# 20210712  GvB       copy attributes of input x to output y
# 20220330  GvB       corrected bug in nfact (default and Sos)
#------------------------------------------------------------------------------

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
#' Before filtering the input signal is extended with a reflected part of both
#' ends of the signal. The length of this extension is 3 times the filter order.
#' The Gustafsson [1] method is then used to specify the initial conditions used
#' to further handle the edges of the signal.
#'
#' @param filt For the default case, the moving-average coefficients of an ARMA
#'   filter (normally called \code{b}). Generically, \code{filt} specifies an
#'   arbitrary filter operation.
#' @param a the autoregressive (recursive) coefficients of an ARMA filter,
#'   specified as a vector. If \code{a[1]} is not equal to 1, then filter
#'   normalizes the filter coefficients by \code{a[1]}. Therefore, \code{a[1]}
#'   must be nonzero.
#' @param x the input signal to be filtered. If \code{x} is a matrix, all
#'   colums are filtered.
#' @param ... additional arguments (ignored).
#'
#' @return The filtered signal, normally of the same length of the input signal
#'   \code{x}, returned as a vector or matrix.
#'
#' @examples
#' bf <- butter(3, 0.1)                                 # 10 Hz low-pass filter
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
#' @seealso \code{\link{filter}}, \code{\link{filter_zi}}, \code{\link{Arma}},
#'   \code{\link{Sos}}, \code{\link{Zpg}}
#'
#' @author Paul Kienzle, \email{pkienzle@@users.sf.net},\cr Francesco Potortì,
#'   \email{pot@@gnu.org},\cr Luca Citi, \email{lciti@@essex.ac.uk}.\cr
#'   Conversion to R and adapted by Geert van Boxtel
#'   \email{G.J.M.vanBoxtel@@gmail.com}.
#'
#' @references [1] Gustafsson, F. (1996). Determining the initial states in
#'   forward-backward filtering. IEEE Transactions on Signal Processing, 44(4),
#'   988 - 992.
#'
#' @rdname filtfilt
#' @export

filtfilt <- function(filt, ...) UseMethod("filtfilt")

#' @rdname filtfilt
#' @method filtfilt default
#' @export

filtfilt.default <- function(filt, a, x, ...) {

  if (!is.vector(filt) || ! is.vector(a) ||
      !is.numeric(filt) || !is.numeric(a)) {
    stop("b and a must be numeric vectors")
  }
  la <- length(a)
  lb <- length(filt)
  lab <- max(la, lb)
  nfact <- max(1, 3 * (lab - 1))  #length of edge transients

  # Compute initial conditions as per [1]
  if (lab > 1) {
    zi <- filter_zi(filt, a)
  } else {
    zi <- NULL
  }

  #save attributes of x
  atx <- attributes(x)
  
  if (is.vector(x)) {
    x <- as.matrix(x, ncol = 1)
    vec <- TRUE
  } else {
    vec <- FALSE
  }
  nrx <- nrow(x)
  ncx <- ncol(x)
  # nfact <- min(nfact - 1, nrx)
  # corrected bug 20220328
  nfact <- min(nfact - 1, nrx - 1)

  y <- matrix(0, nrx, ncx)

  for (icol in seq_len(ncx)) {
    if (nfact > 0) {
      temp <- c(x[seq(nfact + 1, 2, -1), icol], x[, icol],
                x[seq(nrx - 1, nrx - nfact, -1), icol])
      temp <- filter(filt, a, temp, zi * temp[1])$y
      temp <- rev(temp)
      temp <- rev(filter(filt, a, temp, zi * temp[1])$y)
    } else {
      temp <- x[, icol]
      temp <- filter(filt, a, temp)
      temp <- rev(temp)
      temp <- rev(filter(filt, a, temp))
    }

    y[, icol] <- temp[(nfact + 1):(length(temp) - nfact)]
  }

  if (vec) {
    y <- as.vector(y)
  }
  # set attributes of y nd return
  attributes(y) <- atx
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

  if (!is.matrix(filt$sos) || !is.numeric(filt$sos)) {
    stop("sos must be a numeric matrix")
  }
  nfact <- max(1, 3 * length(as.Zpg(filt)$p)) # filter order

  # Compute initial conditions as per [1]
  if (nfact > 1) {
    zi <- filter_zi(filt)
  } else {
    zi <- NULL
  }

  #save attributes of x
  atx <- attributes(x)
  
  if (is.vector(x)) {
    x <- as.matrix(x, ncol = 1)
    vec <- TRUE
  } else {
    vec <- FALSE
  }
  nrx <- nrow(x)
  ncx <- ncol(x)
  # nfact <- min(nfact - 1, nrx)
  # corrected bug 20220328
  nfact <- min(nfact - 1, nrx - 1)
  
  y <- matrix(0, nrx, ncx)

  for (icol in seq_len(ncx)) {
    if (nfact > 0) {
      temp <- c(2 * x[1, icol] - x[seq(nfact + 1, 2, -1), icol], x[, icol],
                2 * x[nrx, icol] - x[seq(nrx - 1, nrx - nfact, -1), icol])
      temp <- filter(filt, temp, zi * temp[1])$y
      temp <- rev(temp)
      temp <- rev(filter(filt, temp, zi * temp[1])$y)
    } else {
      temp <- x[, icol]
      temp <- filter(filt, temp)
      temp <- rev(temp)
      temp <- rev(filter(filt, temp))
    }
    y[, icol] <- temp[(nfact + 1):(length(temp) - nfact)]
  }

  if (vec) {
    y <- as.vector(y)
  }
  # set attributes of y and return
  attributes(y) <- atx
  y
}

#' @rdname filtfilt
#' @method filtfilt Zpg
#' @export
filtfilt.Zpg <- function(filt, x, ...) # zero-pole-gain form
  filtfilt(as.Arma(filt), x)
