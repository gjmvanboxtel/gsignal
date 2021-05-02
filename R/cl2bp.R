# cl2bp.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Octave signal package:
# Copyright (C) 1995, Author: Ivan Selesnick, Rice University
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
# 20200808 Geert van Boxtel           First version for v0.1.0
#------------------------------------------------------------------------------

#' Constrained L2 bandpass FIR filter design
#'
#' Constrained least square band-pass FIR filter design without specified
#' transition bands.
#'
#' This is a fast implementation of the algorithm cited below. Compared to
#' \code{remez}, it offers implicit specification of transition bands, a higher
#' likelihood of convergence, and an error criterion combining features of both
#' L2 and Chebyshev approaches
#'
#' @param m  degree of cosine polynomial, resulting in a filter of length
#'   \code{2 * m + 1}. Must be an even number. Default: 30.
#' @param w1,w2 bandpass filter cutoffs in the range \code{0 <= w1 < w2 <= pi},
#'   where pi is the Nyquist frequency.
#' @param up vector of 3 upper bounds for c(stopband1, passband, stopband2).
#' @param lo vector of 3 lower bounds for c(stopband1, passband, stopband2).
#' @param L search grid size; larger values may improve accuracy, but greatly
#'   increase calculation time. Default: 2048, maximum: 1e6.
#'
#' @return The FIR filter coefficients, a vector of length \code{2 * m + 1}, of
#'   class \code{Ma}.
#'
#' @references Selesnick, I.W., Lang, M., and Burrus, C.S. (1998) A modified
#'   algorithm for constrained least square design of multiband FIR filters
#'   without specified transition bands. IEEE Trans. on Signal Processing,
#'   46(2), 497-501. \cr
#'   \url{https://www.ece.rice.edu/dsp/software/cl2.shtml}
#'
#' @examples
#' w1 <- 0.3 * pi
#' w2 <- 0.6 * pi
#' up <- c(0.02, 1.02, 0.02)
#' lo <- c(-0.02, 0.98, -0.02)
#' h  <- cl2bp(30, w1, w2, up, lo, 2^11)
#' freqz(h)
#'
#' @seealso \code{\link{Ma}}, \code{\link{filter}}, \code{\link{remez}}
#'
#' @author Ivan Selesnick, Rice University, 1995,
#'   downloaded from \url{https://www.ece.rice.edu/dsp/software/cl2.shtml}.\cr
#'   Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#'
#' @export

cl2bp <- function(m = 30, w1, w2, up, lo, L = 2048) {

  # parameter checking
  if (!isPosscal(m) || !isWhole(m) || m %% 2 != 0) {
    stop("polynomial degree m must be a positive even integer")
  }
  if (!is.numeric(w1) || w1 < 0 || w1 >= pi ||
      !is.numeric(w2) || w2 <= 0 || w2 > pi ||
      w1 >= w2) {
    stop("bandpass filter cutoffs must be in the range 0 <= w1 < w2 <= pi")
  }
  if (!is.vector(up) || length(up) != 3 ||
      !is.vector(lo) || length(lo) != 3) {
    stop(paste("the up and lo vectors must contain 3 values",
               "[stopband1, passband, stopband2]"))
  }
  if (!isPosscal(L) || L > 1000000) {
    stop("L must be a positive scalar <= 1000000")
  }

  # ----- calculate Fourier coefficients and upper ---
  # ----- and lower bound functions ------------------

  q1 <- round(L * w1 / pi)
  q2 <- round(L * (w2 - w1) / pi)
  q3 <- L + 1 - q1 - q2
  u <- c(up[1] * rep(1, q1), up[2] * rep(1, q2), up[3] * rep(1, q3))
  l <- c(lo[1] * rep(1, q1), lo[2] * rep(1, q2), lo[3] * rep(1, q3))
  w <- (0:L) * pi / L
  Z <- rep(0L, 2 * L - 1 - 2 * m)
  r <- sqrt(2)
  c <- c((w2 - w1) * r, 2 * (sin(w2 * (1:m)) - sin(w1 * (1:m))) / (1:m)) / pi
  a <- c       	# best L2 cosine coefficients
  mu <- NULL    # Lagrange multipliers
  SN <- 1e-9   	# Small Number
  kmax <- NULL; uvo <- 0
  kmin <- NULL; lvo <- 0

  counter <- 0
  while (TRUE) {
    counter <- counter + 1
    if ((uvo > SN / 2) | (lvo > SN / 2)) {
      # ----- include old extremal ----------------
      if (uvo > lvo) {
        kmax <- c(kmax, okmax[k1]); okmax <- okmax[-k1]
      } else {
        kmin <- c(kmin, okmin[k2]); okmin <- okmin[-k2]
      }
    } else {
      # ----- calculate A -------------------------
      A <- stats::fft(c(a[1] * r, a[2:(m + 1)], Z, a[seq(m + 1, 2, -1)]))
      A <- Re(A[1:(L + 1)]) / 2
      # ----- find extremals ----------------------
      okmax <- kmax;         okmin <- kmin
      kmax  <- local_max(A); kmin  <- local_max(-A)
      kmax  <- kmax[A[kmax] > u[kmax] - SN / 10]
      kmin  <- kmin[A[kmin] < l[kmin] + SN / 10]
      # ----- check stopping criterion ------------
      Eup <- A[kmax] - u[kmax]
      Elo <- l[kmin] - A[kmin]
      E <- max(c(Eup, Elo, 0))
      if (E < SN) break
    }

    # ----- calculate new multipliers -----------
    n1 <- length(kmax); n2 <- length(kmin)
    O  <- rbind(matrix(1L, n1, m + 1), matrix(-1L, n2, m + 1))
    G  <- O * cos(w[c(kmax, kmin)] %o% (0:m))
    G[, 1] <- G[, 1] / r
    d  <- c(u[kmax],  -l[kmin])
    mu <- pracma::mldivide(G %*% t(G), G %*% c - d, pinv = FALSE)
    # ----- remove negative multiplier ----------
    min_mu <- min(mu); K <- which.min(mu)
    while (min_mu < 0) {
      G <- G[-K, ]
      d <- d[-K]
      mu <- pracma::mldivide(G %*% t(G), G %*% c - d, pinv = FALSE)
      if (K > n1) {
        kmin <- kmin[- (K - n1)]; n2 <- n2 - 1
      } else {
        kmax <- kmax[-K]; n1 <- n1 - 1
      }
      min_mu <- min(mu); K <- which.min(mu)
    }
    # ----- determine new coefficients ----------
    a <- c - t(G) %*% mu

    if (length(okmax) > 0) {
      Aokmax <- a[1] / r + cos(w[okmax] %o% (1:m)) %*% a[2:(m + 1)]
      uvo <- max(c(Aokmax - u[okmax], 0))
      k1 <- which.max(c(Aokmax - u[okmax], 0))
    } else {
      uvo <- 0
    }
    if (length(okmin) > 0) {
      Aokmin <- a[1] / r + cos(w[okmin] %o% (1:m)) %*% a[2:(m + 1)]
      lvo <- max(c(l[okmin] - Aokmin, 0))
      k2 <- which.max(c(l[okmin] - Aokmin, 0))
    } else {
      lvo <- 0
    }

  }

  h <- c(a[seq(m + 1, 2, -1)], a[1] * r,  a[2:(m + 1)]) / 2

  # return filter coefficients as Ma structure
  Ma(h)

}

local_max <- function(x) {
  # finds location of local maxima
  N <- length(x)
  b1 <- x[1:(N - 1)] <= x[2:N]
  b2 <- x[1:(N - 1)] >  x[2:N]
  k <- which(b1[1:(N - 2)] & b2[2:(N - 1)]) + 1
  if (x[1] > x[2]) k <- c(k, 1)
  if (x[N] > x[N - 1]) k <- c(k, N)
  k <- sort(k)
  k
}
