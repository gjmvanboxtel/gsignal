# filter.R
# Copyright (C) 2020-2021 Geert van Boxtel <gjmvanboxtel@gmail.com>
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
# 20200208  GvB       setup for gsignal v0.1.0
# 20200413  GvB       added S3 method for Sos
# 20210319  GvB       new setup using filter.cpp to handle initial conditions
#                     the function now also has the options to return the
#                     final conditions
# 20210515  GvB       check return value of rfilter
# 20210712  GvB       copy attributes of input x to output y
# 20210724  GvB       allow filtering complex signals and coefficients
#------------------------------------------------------------------------------

#' Filter a signal
#'
#' Apply a 1-D digital filter compatible with 'Matlab' and 'Octave'.
#'
#' The filter is a direct form II transposed implementation of the standard
#' linear time-invariant difference equation:
#' \if{latex}{
#'   \deqn{\sum_{k=0}^{N} a(k+1) y(n-k) + \sum_{k=0}^{M} b(k+1) x(n-k) = 0; 1
#'   \le n \le length(x)}
#' }
#' \if{html}{\preformatted{
#'   N                  M
#'  SUM a(k+1)y(n-k) + SUM b(k+1)x(n-k) = 0;   1 <= n <= length(x)
#'  k=0                k=0
#' }}
#' where \code{N = length(a) - 1} and \code{M = length(b) - 1}.
#'
#' The initial and final conditions for filter delays can be used to filter data
#' in sections, especially if memory limitations are a consideration. See the
#' examples.
#'
#' @param filt For the default case, the moving-average coefficients of an ARMA
#'   filter (normally called \code{b}), specified as a numeric or complex
#'   vector. Generically, \code{filt} specifies an arbitrary filter operation.
#' @param a the autoregressive (recursive) coefficients of an ARMA filter,
#'   specified as a numeric or complex vector. If \code{a[1]} is not equal to 1,
#'   then filter normalizes the filter coefficients by \code{a[1]}. Therefore,
#'   \code{a[1]} must be nonzero.
#' @param x the input signal to be filtered, specified as a numeric or complex
#'   vector or matrix. If \code{x} is a matrix, each column is filtered.
#' @param zi If \code{zi} is provided, it is taken as the initial state of the
#'   system and the final state is returned as zf. The state vector is a vector
#'   or a matrix (depending on \code{x}) whose length or number of rows is equal
#'   to the length of the longest coefficient vector \code{b} or \code{a} minus
#'   one. If \code{zi} is not supplied (NULL), the initial state vector is set
#'   to all zeros. Alternatively, \code{zi} may be the character string
#'   \code{"zf"}, which specifies to return the final state vector even though
#'   the initial state vector is set to all zeros. Default: NULL.
#' @param ... additional arguments (ignored).
#'
#' @return The filtered signal, of the same dimensions as the input signal. In
#'   case the \code{zi} input argument was specified, a list with two elements
#'   is returned containing the variables \code{y}, which represents the output
#'   signal, and \code{zf}, which contains the final state vector or matrix.
#'
#' @examples
#' bf <- butter(3, 0.1)                                 # 10 Hz low-pass filter
#' t <- seq(0, 1, len = 100)                            # 1 second sample
#' x <- sin(2* pi * t * 2.3) + 0.25 * rnorm(length(t))  # 2.3 Hz sinusoid+noise
#' z <- filter(bf, x)                                   # apply filter
#' plot(t, x, type = "l")
#' lines(t, z, col = "red")
#'
#' ## specify initial conditions
#' ## from Python scipy.signal.lfilter() documentation
#' t <- seq(-1, 1, length.out =  201)
#' x <- (sin(2 * pi * 0.75 * t * (1 - t) + 2.1)
#'       + 0.1 * sin(2 * pi * 1.25 * t + 1)
#'       + 0.18 * cos(2 * pi * 3.85 * t))
#' h <- butter(3, 0.05)
#' lab <- max(length(h$b), length(h$a)) - 1
#' zi <- filtic(h$b, h$a, rep(1, lab), rep(1, lab))
#' z1 <- filter(h, x)
#' z2 <- filter(h, x, zi * x[1])
#' plot(t, x, type = "l")
#' lines(t, z1, col = "red")
#' lines(t, z2$y, col = "green")
#' legend("bottomright", legend = c("Original signal",
#'         "Filtered without initial conditions",
#'         "Filtered with initial conditions"),
#'        lty = 1, col = c("black", "red", "green"))
#'
#' @seealso \code{\link{filter_zi}}, \code{\link{sosfilt}} (preferred because it
#'   avoids numerical problems).
#'
#' @author Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#'
#' @rdname filter
#' @export

filter <- function(filt, ...) UseMethod("filter")

#' @rdname filter
#' @method filter default
#' @export

filter.default <- function(filt, a, x, zi = NULL, ...) {
  
  if (!is.vector(filt) || ! is.vector(a)) {
    stop("b and a must be numeric vectors")
  }
  if (is.numeric(filt) && is.numeric(a)) {
    real_coefs <- TRUE
  } else if (is.complex(filt) || is.complex(a)) {
    real_coefs <- FALSE
  } else {
    stop("b and a must be numeric or complex")
  }
  la <- length(a)
  lb <- length(filt)
  lab <- max(la, lb)
  
  #save attributes of x
  atx <- attributes(x)
  
  if (is.null(x)) {
    return(NULL)
  }
  if (is.numeric(x)) {
    real_x <- TRUE
  } else if (is.complex(x)) {
    real_x <- FALSE
  } else {
    stop("x must be a numeric or complex vector or matrix")
  }
  if (is.vector(x)) {
    x <- as.matrix(x, ncol = 1)
    vec <- TRUE
  } else {
    vec <- FALSE
  }
  nrx <- NROW(x)
  ncx <- NCOL(x)
  if (is.null(nrx) || nrx <= 0) {
    return(x)
  }
  
  rzf <- (is.character(zi) && zi == "zf")
  if (!is.null(zi) && !is.numeric(zi) && ! is.complex(zi) && !rzf) {
    stop("invalid value for zi")
  }
  if (is.null(zi) || rzf) {
    if (is.numeric(x)) {
      zi <- matrix(0, lab - 1, ncx)
    } else if(is.complex(x)) {
      zi <- matrix(0 + 0i, lab - 1, ncx)
    }
    if (is.null(zi)) {
      rzf <- FALSE
    }
  } else if (!is.null(zi)) {
    rzf <- TRUE
  }
  if (is.vector(zi)) {
    zi <- as.matrix(zi, ncol = 1)
  }
  nrzi <- NROW(zi)
  nczi <- NCOL(zi)
  if (nrzi != lab - 1) {
    stop("zi must be of length max(length(a), length(b)) - 1")
  }
  if (nczi != ncx) {
    stop("number of columns of zi and x must agree")
  }
  
  while (length(a) > 1 && a[1] == 0) {
    a <- a[2:length(a)]
  }
  if (length(a) < 1 || a[1] == 0) {
    stop("There must be at least one nonzero element in the vector a")
  }
  if (a[1] != 1) {
    # Normalize the coefficients so a[1] == 1.
    filt <- filt / a[1]
    a <- a / a[1]
  }
  if (la <= 1 && nrzi <= 0) {
    retval <- filt[1] * x
    if (vec) retval <- as.vector(retval)
    return(retval)
  }
  
  if (real_coefs) {
    if (real_x) {
      # real filter coefs, real x
      y <- matrix(0, nrx, ncx)
      zf <- matrix(0, lab - 1, nczi)
      for (icol in seq_len(ncx)) {
        l <- .Call("_gsignal_rfilter", PACKAGE = "gsignal",
                   filt, a, x[, icol], zi[, icol])
        if (length(l) <= 0) {
          stop("Error filtering data")
        }
        y[, icol] <- l[["y"]]
        zf[, icol] <- l[["zf"]]
      }
    } else {
      # real filter coefs, complex x
      y <- matrix(0 + 0i, nrx, ncx)
      zf <- matrix(0 + 0i, lab - 1, nczi)
      for (icol in seq_len(ncx)) {
        re <- .Call("_gsignal_rfilter", PACKAGE = "gsignal",
                    filt, a, Re(x[, icol]), Re(zi[, icol]))
        if (length(re) <= 0) {
          stop("Error filtering data")
        }
        im <- .Call("_gsignal_rfilter", PACKAGE = "gsignal",
                    filt, a, Im(x[, icol]), Im(zi[, icol]))
        if (length(im) <= 0) {
          stop("Error filtering data")
        }
        y[, icol] <- re[["y"]] + 1i * im[["y"]]
        zf[, icol] <- re[["zf"]] + 1i * im[["zf"]]
      }
    }
  } else {
    if (real_x) {
      # complex filter coefs, real x
      y <- matrix(0 + 0i, nrx, ncx)
      zf <- matrix(0 + 0i, lab - 1, nczi)
      for (icol in seq_len(ncx)) {
        re <- .Call("_gsignal_rfilter", PACKAGE = "gsignal",
                    Re(filt), Re(a), x[, icol], Re(zi[, icol]))
        if (length(re) <= 0) {
          stop("Error filtering data")
        }
        im <- .Call("_gsignal_rfilter", PACKAGE = "gsignal",
                    Im(filt), Im(a), x[, icol], Im(zi[, icol]))
        if (length(im) <= 0) {
          stop("Error filtering data")
        }
        y[, icol] <- re[["y"]] + 1i * im[["y"]]
        zf[, icol] <- re[["zf"]] + 1i * im[["zf"]]
      }
    } else {
      # complex filter coefs, complex x
      # avoid one filtering operation
      y <- matrix(0 + 0i, nrx, ncx)
      zf <- matrix(0 + 0i, lab - 1, nczi)
      for (icol in seq_len(ncx)) {
        l1 <- .Call("_gsignal_rfilter", PACKAGE = "gsignal",
                    Re(filt), Re(a), Re(x[, icol]), Re(zi[, icol]))
        if (length(l1) <= 0) {
          stop("Error filtering data")
        }
        l2 <- .Call("_gsignal_rfilter", PACKAGE = "gsignal",
                    Im(filt), Im(a), Im(x[, icol]), Im(zi[, icol]))
        if (length(l2) <= 0) {
          stop("Error filtering data")
        }
        l3 <- .Call("_gsignal_rfilter", PACKAGE = "gsignal",
                    Re(filt) + Im(filt), Re(a) + Im(a),
                    Re(x[, icol]) + Im(x[, icol]), Re(zi[, icol]) + Im(zi[, icol]))
        if (length(l3) <= 0) {
          stop("Error filtering data")
        }
        y[, icol] <- l1[["y"]] - l2[["y"]] +
          + 1i * (l3[["y"]] - l1[["y"]] - l2[["y"]])
        zf[, icol] <- l1[["zf"]] - l2[["zf"]] +
          + 1i * (l3[["zf"]] - l1[["zf"]] - l2[["zf"]])
      }
    }
  }
  
  if (vec) {
    y <- as.vector(y)
    zf <- as.vector(zf)
  }
  if (!real_x || !real_coefs) {
    y <- zapIm(y)
    zf <- zapIm(zf)
  }
  # set attributes of y nd return
  attributes(y) <- atx
  
  if (rzf) {
    retval <- list(y = y, zf = zf)
  } else {
    retval <- y
  }
  
  retval
  
}

#' @rdname filter
#' @method filter Arma
#' @export
filter.Arma <- function(filt, x, ...) # IIR
  filter(filt$b, filt$a, x, ...)

#' @rdname filter
#' @method filter Ma
#' @export
filter.Ma <- function(filt, x, ...) # FIR
  filter(unclass(filt), 1, x, ...)

#' @rdname filter
#' @method filter Sos
#' @export
filter.Sos <- function(filt, x, ...) { # Second-order sections
  if (filt$g != 1) {
    filt$sos[1, 1:3] <- filt$sos[1, 1:3] * filt$g
  }
  sosfilt(filt$sos, x, ...)
}

#' @rdname filter
#' @method filter Zpg
#' @export
filter.Zpg <- function(filt, x, ...) # zero-pole-gain form
  filter(as.Arma(filt), x, ...)

