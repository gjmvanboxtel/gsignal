# xcorr.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Octave version:
# Copyright (C) 1999-2001 Paul Kienzle <pkienzle@users.sf.net>
# Copyright (C) 2004 <asbjorn.sabo@broadpark.no>
# Copyright (C) 2008, 2010 Peter Lanspeary <peter.lanspeary@.adelaide.edu.au>
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
# 2020313  GvB       setup for gsignal v0.1.0
#------------------------------------------------------------------------------

#' Cross-correlation
#'
#' Estimate the cross-correlation between two sequences or the autocorrelation
#' of a single sequence
#'
#' Estimate the cross correlation R_xy(k) of vector arguments \code{x} and
#' \code{y} or, if \code{y} is omitted, estimate autocorrelation R_xx(k) of
#' vector \code{x}, for a range of lags \code{k} specified by the argument
#' \code{maxlag}. If \code{x} is a matrix, each column of \code{x} is correlated
#' with itself and every other column.
#'
#' The cross-correlation estimate between vectors \code{x} and \code{y} (of
#' length \code{N}) for lag \code{k} is given by
#' \if{latex}{
#'   \deqn{R_{xy}(k) = \sum_{i=1}^{N} x_{i+k} Conj(y_i)}
#' }
#' \if{html}{\preformatted{
#'             N
#'      Rxy = SUM x(i+k) . Conj(y(i))
#'            i=1
#' }}
#' where data not provided (for example \code{x[-1], y[N+1]}) is zero. Note the
#' definition of cross-correlation given above. To compute a cross-correlation
#' consistent with the field of statistics, see xcov.
#'
#' The cross-correlation estimate is calculated by a "spectral" method
#' in which the FFT of the first vector is multiplied element-by-element
#' with the FFT of second vector.  The computational effort depends on
#' the length N of the vectors and is independent of the number of lags
#' requested.  If you only need a few lags, the "direct sum" method may
#' be faster.
#'
#' @param x Input, numeric or complex vector or matrix. Must not be missing.
#' @param y Input, numeric or complex vector data.  If \code{x} is a matrix (not
#'   a vector), \code{y} must be omitted. \code{y} may be omitted if \code{x} is
#'   a vector; in this case \code{xcorr} estimates the autocorrelation of
#'   \code{x}.
#' @param maxlag Integer scalar. Maximum correlation lag. If omitted, the
#'   default value is \code{N-1}, where \code{N} is the greater of the lengths
#'   of \code{x} and \code{y} or, if \code{x} is a matrix, the number of rows in
#'   \code{x}.
#' @param scale Character string. Specifies the type of scaling applied to the
#'   correlation vector (or matrix). matched to one of:
#'   \describe{
#'     \item{"none"}{return the unscaled correlation, R}
#'     \item{"biased"}{return the biased average, R / N}
#'     \item{"unbiased"}{return the unbiased average, R(k) / (N - |k|)}
#'     \item{"coeff"}{return the correlation coefficient, R / (rms(x) .
#'     rms(y))}, where \code{k} is the lag, and \code{N} is the length of
#'     \code{x}
#'   }
#'  If omitted, the default value is \code{"none"}. If \code{y} is supplied but
#'  does not have the same length as \code{x}, scale must be \code{"none"}.
#'
#' @return A list containing the following variables:
#' \describe{
#'   \item{R}{array of correlation estimates}
#'   \item{lags}{vector of correlation lags \code{[-maxlag:maxlag]}}
#' }
#' The array of correlation estimates has one of the following forms:
#' \enumerate{
#'   \item Cross-correlation estimate if X and Y are vectors.
#'   \item Autocorrelation estimate if is a vector and Y is omitted.
#'   \item If \code{x} is a matrix, \code{R} is a matrix containing the
#'   cross-correlation estimate of each column with every other column. Lag
#'   varies with the first index so that \code{R} has \code{2 * maxlag + 1} rows
#'   and \eqn{P^2} columns where \code{P} is the number of columns in \code{x}.
#' }
#' @seealso \code{\link{xcov}}.
#'
#' @examples
#' ## Create a vector x and a vector y that is equal to x shifted by 5
#' ## elements to the right. Compute and plot the estimated cross-correlation
#' ## of x and y. The largest spike occurs at the lag value when the elements
#' ## of x and y match exactly (-5).
#' n <- 0:15
#' x <- 0.84^n
#' y <- pracma::circshift(x, 5)
#' rl <- xcorr(x, y)
#' plot(rl$lag, rl$R, type="h")
#'
#' ## Compute and plot the estimated autocorrelation of a vector x.
#' ## The largest spike occurs at zero lag, when x matches itself exactly.
#' n <- 0:15
#' x <- 0.84^n
#' rl <- xcorr(x)
#' plot(rl$lag, rl$R, type="h")
#'
#' ## Compute and plot the normalized cross-correlation of vectors
#' ## x and y with unity peak, and specify a maximum lag of 10.
#' n <- 0:15
#' x <- 0.84^n
#' y <- pracma::circshift(x, 5)
#' rl <- xcorr(x, y, 10, 'coeff')
#' plot(rl$lag, rl$R, type="h")
#'
#' @author Paul Kienzle, \email{pkienzle@@users.sf.net},\cr
#'   Asbjorn Sabo, \email{asbjorn.sabo@@broadpark.no},\cr
#'   Peter Lanspeary. \email{peter.lanspeary@@adelaide.edu.au}.\cr
#'    Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#'
#' @export

xcorr <- function(x, y = NULL,
                  maxlag = if (is.matrix(x)) nrow(x) - 1
                  else max(length(x), length(y)) - 1,
                  scale = c("none", "biased", "unbiased", "coeff")) {

  if (is.array(x)) {
    ld <- length(dim(x))
    if (ld == 1) {
      x <- as.vector(x)
    } else if (ld == 2) {
      x <- as.matrix(x)
    } else {
      stop("multidimensional arrays are not supported (x)")
    }
  }
  if (!is.null(y) && is.array(y)) {
    ld <- length(dim(y))
    if (ld == 1) {
      y <- as.vector(y)
    } else if (ld == 2) {
      y <- as.matrix(y)
    } else {
      stop("multidimensional arrays are not supported (y)")
    }
  }

  if (is.vector(x)) {
    N <- max(length(x), length(y))
  } else {
    N <- nrow(x)
  }
  scale <- match.arg(scale)

  ## check argument values
  if (!(is.vector(x) || is.matrix(x)) || !(is.numeric(x) || is.complex(x))) {
    stop("x must be a numeric or complex vector or matrix")
  }
  if (!is.null(y) && is.matrix(y) && !(is.numeric(y) || is.complex(y))) {
    stop("y must be a vector")
  }
  if (!is.null(y) && !is.vector(x)) {
    stop("x must be a vector if y is specified")
  }
  if (!isPosscal(maxlag) || !isWhole(maxlag)) {
    stop("maxlag must be a non-negative integer")
  }

  # Correlations for lags in excess of +/-(N-1)
  #  (a) are not calculated by the FFT algorithm,
  #  (b) are all zero;
  # so provide them by padding the results (with zeros) before returning.
  if (maxlag > (N - 1)) {
    pad_result <- maxlag - (N - 1)
    maxlag <- N - 1
  } else {
    pad_result <- 0
  }

  if (is.vector(x) && is.vector(y) &&
      length(x) != length(y) && scale != "none") {
    stop("scale must be 'none' if length(x) != length(y)")
  }

  P <- ncol(x)
  M <- nextpow2(N + maxlag)
  if (!is.vector(x)) {
    # correlate each column "i" with all other "j" columns
    R <- matrix(0L, 2 * maxlag + 1, P^2)
    # do FFTs of padded column vectors
    pre <- stats::mvfft(postpad(prepad(x, N + maxlag), M))
    post <- Conj(stats::mvfft(postpad(x, M)))
    # do autocorrelations (each column with itself)
    cor <- imvfft(post * pre)
    R[, seq(1, P^2, P + 1)] <- cor[1:(2 * maxlag + 1), ]
    # do the cross correlations
    for (i in 1:(P - 1)) {
      j <- (i + 1):P
      if (length(j) > 1) {
        cor <- imvfft(pre[, i * rep(1L, length(j))] * post[, j])
        R[, ((i - 1) * P + j)] <- cor[1:(2 * maxlag + 1), ]
        R[, ((j - 1) * P + i)] <-
          Conj(pracma::flipud(cor[1:(2 * maxlag + 1), ]))
      } else {
        cor <- ifft(pre[, i * rep(1L, length(j))] * post[, j])
        R[, ((i - 1) * P + j)] <- cor[1:(2 * maxlag + 1)]
        R[, ((j - 1) * P + i)] <- Conj(rev(cor[1:(2 * maxlag + 1)]))
      }

    }
  } else if (is.null(y)) {
    # compute autocorrelation of a single vector
    post <- stats::fft(postpad(x, M))
    cor <- ifft(post * Conj(post))
    R <- c(Conj(cor[seq((maxlag + 1), 2, -1)]), cor[1:(maxlag + 1)])
  } else {
    # compute cross-correlation of x and y
    pre <- stats::fft(postpad(prepad(x, length(x) + maxlag), M))
    post <- stats::fft(postpad(y, M))
    cor <- ifft(pre * Conj(post))
    R <- cor[1:(2 * maxlag + 1)]
  }

  # if inputs are real, outputs should be real, so ignore the
  # insignificant complex portion left over from the FFT
  if (is.numeric(x) && (is.null(y) || is.numeric(y))) {
    dr <- dim(R)
    R <- as.numeric(R)
    dim(R) <- dr
  }

  # correct for bias
  if (scale == "biased") {
    R <- R / N
  } else if (scale == "unbiased") {
    R <- R / (c((N - maxlag):(N - 1), N,
                rev((N - maxlag):(N - 1)))) * rep(1L, NCOL(R))
  } else if (scale == "coeff") {
    ## R = R ./ R(maxlag+1) works only for autocorrelation
    ## For cross correlation coeff, divide by rms(X)*rms(Y).
    if (!is.vector(x)) {
      ## for matrix (more than 1 column) X
      rms <- sqrt(ssq(x))
      R <- R / (rep(1L, nrow(R)) * rms)
    } else if (is.null(y)) {
      ##  for autocorrelation, R(zero-lag) is the mean square.
      R <- R / R[maxlag + 1]
    } else {
      ##  for vectors X and Y
      R <- R / sqrt(ssq(x) * ssq(y))
    }
  }

  ## Pad result if necessary
  ##  (most likely is not required, use "if" to avoid unnecessary code)
  ## At this point, lag varies with the first index in R;
  ##  so pad **before** the transpose.
  if (pad_result) {
    R_pad <- matrix(0L, pad_result, ncol(R))
    R <- cbind(R_pad, R, R_pad)
  }

  maxlag <- maxlag + pad_result
  lags <- -maxlag:maxlag

  list(R = R, lags = lags)
}
