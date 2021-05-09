# dct.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Original Octave code:
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
# 20201015  GvB       setup for gsignal v0.1.0
# 20210506  GvB       use matrix() instead of as.matrix
#------------------------------------------------------------------------------

#' Inverse Discrete Cosine Transform
#'
#' Compute the inverse unitary discrete cosine transform of a signal.
#'
#' The discrete cosine transform (DCT) is closely related to the discrete
#' Fourier transform. You can often reconstruct a sequence very accurately from
#' only a few DCT coefficients. This property is useful for applications
#' requiring data reduction.
#'
#' @param x input discrete cosine transform, specified as a numeric vector or
#'   matrix. In case of a vector it represents a single signal; in case of a
#'   matrix each column is a signal.
#' @param n transform length, specified as a positive integer scalar. Default:
#'   \code{NROW(x)}.
#'
#' @return Inverse discrete cosine transform, returned as a vector or matrix.
#'
#' @examples
#' x <- seq_len(100) + 50 * cos(seq_len(100) * 2 * pi / 40)
#' X <- dct(x)
#'
#' # Find which cosine coefficients are significant (approx.)
#' # zero the rest
#' nsig <- which(abs(X) < 1)
#' N <- length(X) - length(nsig) + 1
#' X[nsig] <- 0
#'
#' # Reconstruct the signal and compare it to the original signal.
#' xx <- idct(X)
#' plot(x, type = "l")
#' lines(xx, col = "red")
#' legend("bottomright", legend = c("Original", paste("Reconstructed, N =", N)),
#'        lty = 1, col = 1:2)
#'
#' @author Paul Kienzle, \email{pkienzle@@users.sf.net}.\cr
#' Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#'
#' @seealso \code{\link{dct}}
#'
#' @export

idct <- function(x, n = NROW(x)) {

  # check parameters
  if (!(is.vector(x) || is.matrix(x)) || !(is.numeric(x) || is.complex(x))) {
    stop("x must be a numeric or complex vector or matrix")
  } else {
    realx <- is.numeric(x)
  }

  if (is.vector(x)) {
    vec <- TRUE
    x <- matrix(x, ncol = 1)
  } else {
    vec <- FALSE
  }
  nr <- nrow(x)
  ns <- ncol(x)

  if (!isPosscal(n) || !isWhole(n)) {
    stop("n must be a positive integer")
  }


  if (n != nr) {
    x <- postpad(x, n)
  }

  if (realx && n %% 2 == 0) {
    w <- c(sqrt(n / 4), sqrt(n / 2) * exp((1i * pi / 2 / n) *
                                            seq_len(n - 1))) %o% rep(1, ns)
    y <- imvfft(w * x)
    y[c(seq(1, n, 2), seq(n, 1, -2)), ] <- 2 * Re(y)
  } else if (n == 1) {
    y <- x
  } else {
    ## reverse the steps of dct using inverse operations
    ## 1. undo post-fft scaling
    w <- c(sqrt(4 * n), sqrt(2 * n) * exp((1i * pi / 2 / n) *
                                            seq_len(n - 1))) %o% rep(1, ns)
    y <- x * w
    ## 2. reconstruct fft result and invert it
    w <- exp(-1i * pi * seq(n - 1, 1, -1) / n) %o% rep(1, ns)
    y <- imvfft(rbind(y, rep(0, ns),
                      matrix(y[seq(n, 2, -1), ], ncol = ns) * w))
    ## 3. keep only the original data; toss the reversed copy
    y <- y[1:n, ]
  }

  if (realx) {
    y <- Re(y)
  }

  if (vec) {
    y <- as.vector(y)
  }
  y
}
