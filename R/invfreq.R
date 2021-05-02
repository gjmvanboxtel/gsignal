# invfreq.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Octave signal package:
# Copyright (C) 1986, 2000, 2003 Julius O. Smith III <jos@ccrma.stanford.edu>
# Copyright (C) 2007 Rolf Schirmacher <Rolf.Schirmacher@MuellerBBM.de>
# Copyright (C) 2003 Andrew Fitting
# Copyright (C) 2010 Pascal Dupuis <Pascal.Dupuis@uclouvain.be>
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
# 20201118 Geert van Boxtel          First version for v0.1.0
#------------------------------------------------------------------------------

#' Inverse Frequency Response
#'
#' Identify filter parameters from frequency response data.
#'
#' Given a desired (one-sided, complex) spectrum \code{h(w)} at equally spaced
#' angular frequencies \eqn{w = (2 \pi k) / N}, k = 0, ... N-1, this function
#' finds the filter \code{B(z)/A(z)} or \code{B(s)/A(s)} with \code{nb} zeroes
#' and \code{na} poles. Optionally, the fit-errors can be weighted with respect
#' to frequency according to the weights \code{wt}.
#'
#' @param h Frequency response, specified as a vector
#' @param w Angular frequencies at which \code{h} is computed, specified as a
#'   vector
#' @param nb,na Desired order of the numerator and denominator polynomials,
#'   specified as positive integers.
#' @param wt Weighting factors, specified as a vector of the same length as
#'   \code{w}. Default: \code{rep(1, length(w))}
#' @param plane \code{"z"} (default) for discrete-time spectra; \code{"s"} for
#'   continuous-time spectra
#' @param method minimization method used to solve the normal equations, one of:
#'   \describe{
#'     \item{"ols"}{ordinary least squares (default)}
#'     \item{"tls"}{total least squares}
#'     \item{"qr"}{QR decomposition}
#'   }
#' @param norm logical indicating whether frequencies must be normalized to
#'   avoid matrices with rank deficiency. Default: TRUE
#'
#' @return A list of class \code{'Arma'} with the following list elements:
#' \describe{
#'   \item{b}{moving average (MA) polynomial coefficients}
#'   \item{a}{autoregressive (AR) polynomial coefficients}
#' }
#'
#' @examples
#' order <- 6  # order of test filter
#' fc <- 1/2   # sampling rate / 4
#' n <- 128    # frequency grid size
#' ba <- butter(order, fc)
#' hw <- freqz(ba, n)
#' BA = invfreq(hw$h, hw$w, order, order)
#' HW = freqz(BA, n)
#' plot(hw$w, abs(hw$h), type = "l", xlab = "Frequency (rad/sample)",
#'   ylab = "Magnitude")
#' lines(HW$w, abs(HW$h), col = "red")
#' legend("topright", legend = c("Original", "Measured"), lty  = 1, col = 1:2)
#' err <- norm(hw$h - HW$h, type = "2")
#' title(paste('L2 norm of frequency response error =', err))
#'
#' @author Julius O. Smith III, Rolf Schirmacher, Andrew Fitting, Pascal
#'   Dupuis.\cr Conversion to R by Geert van Boxtel,
#'   \email{G.J.M.vanBoxtel@@gmail.com}.
#'
#' @references
#'   \url{https://ccrma.stanford.edu/~jos/filters/FFT_Based_Equation_Error_Method.html}
#'
#' @rdname invfreq
#' @export

invfreq <- function(h, w, nb, na,
                    wt = rep(1, length(w)),
                    plane = c("z", "s"),
                    method = c("ols", "tls", "qr"),
                    norm = TRUE) {

  # Parameter checking
  if (!is.vector(h)) {
    stop("h must be a vector")
  }
  if (!is.vector(w)) {
    stop("h must be a vector")
  }
  nw <- length(w)
  if (length(h) != nw) {
    stop("h and f must be of equal length")
  }
  if (!isPosscal(nb) || !isWhole(nb) || nb <= 0 ||
     !isPosscal(na) || !isWhole(na) || na <= 0) {
    stop("na and nb must be positive integers > 0")
  }
  if (length(nb) > 1) {
    zb <- nb[2]
    nb <- nb[1]
  } else {
    zb <- 0
  }
  n <- max(na, nb)
  ma <- na + 1; mb <- nb + 1
  if (!is.vector(wt) || length(wt) != nw) {
    stop("wt must be a vector of the same length as w")
  }
  plane <- match.arg(plane)
  method <- match.arg(method)
  norm <- is.logical(norm)
  # End of parameter checking

  Ruu <- matrix(0, mb, mb)
  Ryy <- matrix(0, na, na)
  Ryu <- matrix(0, na, mb)
  Pu <- rep(0, mb)
  Py <- rep(0, na)

  s <- 1i * w
  if (plane == "z") {
    if (max(w) > pi || min(w) < 0) {
      # frequency is outside the range 0 to pi
      w <- seq(0, pi, length.out = length(h))
      s <- 1i * w
    }
    s <- exp(-s)
  } else if (plane == "s") {
    wmax <- max(w)
    if (wmax > 1e6 && n > 5 && !norm) {
      warning("Be careful, there are risks of generating singular matrices")
      warning("Use norm = TRUE to avoid it")
    }
    if (norm) {
      s <- 1i * w / wmax
    }
  }

  for (k in seq_len(nw)) {
    Zk <- s[k]^seq(0, n)
    Hk <- h[k]
    aHks <- Hk * Conj(Hk)
    Rk <- (wt[k] * Zk) %*% t(Zk)
    rRk <- Re(Rk)
    Ruu <- Ruu + rRk[1:mb, 1:mb]
    Ryy <- Ryy + aHks * rRk[2:ma, 2:ma]
    Ryu <- Ryu + Re(Hk * Rk[2:ma, 1:mb])
    Pu <- Pu + wt[k] * Re(Conj(Hk) * Zk[1:mb])
    Py <- Py + (wt[k] * aHks) * Re(Zk[2:ma])
  }
  Rr <- matrix(1, length(s), mb + na)
  Zk <- s
  for (k in seq_len(min(na, na))) {
    Rr[, (1 + k)] <- Zk
    Rr[, (mb + k)] <- -Zk * h
    Zk <- Zk * s
  }
  from <- 1 + min(na, nb)
  to <- max(na, nb) - 1
  if (to >= from) {
    for (k in seq(from, to)) {
      if (k <= nb) {
        Rr[, (1 + k)] <- Zk
      }
      if (k <= na) {
        Rr[, (mb + k)] <- -Zk * h
      }
      Zk <- Zk * s
    }
  }
  k <- k + 1
  if (k <= nb) {
    Rr[, (1 + k)] <- Zk
  }
  if (k <= na) {
    Rr[, (mb + k)] <- -Zk * h
  }

  ## complex to real equation system -- this ensures real solution
  Rr <- Rr[, (1 + zb):ncol(Rr)]
  Rr <- rbind(Re(Rr), Im(Rr))
  Pr <- c(Re(h), Im(h))

  if (method == "ols") {
    R <- qr.R(qr(cbind(Rr, Pr)))
    Theta <- pracma::mldivide(R[1:(nrow(R) - 1), 1:(ncol(R) - 1)],
                              R[1:(nrow(R) - 1), ncol(R)])
  } else if (method == "tls") {
    SVD <- svd(cbind(Rr, Pr))
    V <- SVD$v
    Theta <-  -V[1:(nrow(V) - 1), ncol(V)] / V[nrow(V), ncol(V)]
  } else if (method == "qr") {
    R <- qr.R(qr(cbind(Rr, Pr)))
    eb <- mb - zb
    sa <- eb + 1
    SVD <- svd(R[sa:nrow(R), sa:ncol(R)])
    V <- SVD$v
    Theta <- -V[1:(nrow(V) - 1), ncol(V)] / V[nrow(V), ncol(V)]
    Theta <- c(
      pracma::mldivide(R[1:eb, 1:eb],
                       (R[1:eb, ncol(R)] - R[1:eb, sa:(ncol(R) - 1)] %*%
                          Theta)), Theta)
  } else {
    stop("unknown method")  #can never happen
  }

  B <- c(rep(0, zb), Theta[1:(mb - zb)])
  A <- c(1, Theta[mb - zb + (1:na)])

  if (plane == "s") {
    B <- rev(B)
    A <- rev(A)
    if (norm) {
      Zk <- wmax ^ seq(n, 0, -1)
      for (k in seq(nb, 1 + zb, -1)) {
        B[k] <- B[k] / Zk[k]
      }
      for (k in seq(na, 1, -1)) {
        A[k] <- A[k] / Zk[k]
      }
    }
  }
  Arma(B, A)
}
#' @rdname invfreq
#' @export

invfreqs <- function(h, w, nb, na, wt = rep(1, length(w)),
                     method = c("ols", "tls", "qr"), norm = TRUE) {
  invfreq(h, w, nb, na, wt, plane = "s", method, norm)
}

#' @rdname invfreq
#' @export

invfreqz <- function(h, w, nb, na, wt = rep(1, length(w)),
                     method = c("ols", "tls", "qr"), norm = TRUE) {
  invfreq(h, w, nb, na, wt, plane = "z", method, norm)
}
