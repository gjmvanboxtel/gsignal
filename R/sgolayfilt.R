# sgolayfilt.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Original Matlab/Octave version:
# Copyright (C) 2001 Paul Kienzle <pkienzle@users.sf.net>
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
# 20200322    GvB       setup for gsignal v0.1.0
# 20210405    GvB       if x is a matrix, filter its columns
# 20220328    GvB       copy dimnames of x to output object
#------------------------------------------------------------------------------

#' Savitzky-Golay filtering
#'
#' Filter a signal with a Savitzky-Golay FIR filter.
#'
#' Savitzky-Golay smoothing filters are typically used to "smooth out" a noisy
#' signal whose frequency span (without noise) is large. They are also called
#' digital smoothing polynomial filters or least-squares smoothing filters.
#' Savitzky-Golay filters perform better in some applications than standard
#' averaging FIR filters, which tend to filter high-frequency content along with
#' the noise. Savitzky-Golay filters are more effective at preserving high
#' frequency signal components but less successful at rejecting noise.
#'
#' Savitzky-Golay filters are optimal in the sense that they minimize the
#' least-squares error in fitting a polynomial to frames of noisy data.
#'
#' @param x the input signal to be filtered, specified as a vector or as a
#'   matrix. If \code{x} is a matrix, each column is filtered.
#' @param p Polynomial filter order; must be smaller than \code{n}.
#' @param n Filter length; must a an odd positive integer.
#' @param m Return the m-th derivative of the filter coefficients. Default: 0
#' @param ts Scaling factor. Default: 1
#' @param filt Filter characteristics, usually the result of a call to
#'   \code{sgolay}
#' @param ... Additional arguments (ignored)
#'
#' @return The filtered signal, of the same dimensions as the input signal.
#'
#' @examples
#' # Compare a 5 sample averager, an order-5 butterworth lowpass
#' # filter (cutoff 1/3) and sgolayfilt(x, 3, 5), the best cubic
#' # estimated from 5 points.
#' bf <- butter(5, 1/3)
#' x <- c(rep(0, 15), rep(10, 10), rep(0, 15))
#' sg <- sgolayfilt(x)
#' plot(sg, type="l", xlab = "", ylab = "")
#' lines(filtfilt(rep(1, 5) / 5, 1, x), col = "red") # averaging filter
#' lines(filtfilt(bf, x), col = "blue")              # butterworth
#' points(x, pch = "x")                              # original data
#' legend("topleft", c("sgolay (3,5)", "5 sample average", "order 5
#' Butterworth", "original data"), lty=c(1, 1, 1, NA),
#' pch = c(NA, NA, NA, "x"), col = c(1, "red", "blue", 1))
#'
#' @seealso \code{\link{sgolay}}
#'
#' @author Paul Kienzle, \email{pkienzle@@users.sf.net}.\cr
#' Conversion to R Tom Short,\cr
#' adapted by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#'
#' @rdname sgolayfilt
#' @export

filter.sgolayFilter <- function(filt, x, ...) {
  sgolayfilt(x, filt, ...)
}

#' @rdname sgolayfilt
#' @export

sgolayfilt <- function(x, p = 3, n = p + 3 - p %% 2, m = 0, ts = 1) {

  if (is.null(x)) {
    return(NULL)
  }
  if (!is.numeric(x)) {
    stop("x must be a numeric vector or matrix")
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
  y <- matrix(0, nrx, ncx)

  ## The first k rows of F are used to filter the first k points
  ## of the data set based on the first n points of the data set.
  ## The last k rows of F are used to filter the last k points
  ## of the data set based on the last n points of the dataset.
  ## The remaining data is filtered using the central row of F.
  ## As the filter coefficients are used in the reverse order of what
  ## seems the logical notation, reverse F[k+1,] so that antisymmetric
  ## sequences are used with the right sign.
  if ("sgolayFilter" %in% class(p) || (!is.null(dim(p)) && dim(p) > 1)) {
    Fm <- p
    n <- nrow(Fm)
  } else {
    Fm <- sgolay(p, n, m, ts)
  }
  k <- floor(n / 2)
  z <- filter(Fm[(k + 1), n:1], 1, x)
  for (icol in seq_len(ncx)) {
    y[, icol] <- c(Fm[1:k, ] %*% x[1:n, icol],
                   z[n:nrx, icol],
                   Fm[(k + 2):n, ] %*% x[(nrx - n + 1):nrx, icol]
                  )
  }
  if (vec) {
    y <- as.vector(y)
  }
  dimnames(y) <- dimnames(x)
  y
}
