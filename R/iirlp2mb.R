# iirlp2mb.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Octave signal package:
# Copyright (C) 2011 Alan J. Greenberger <alanjg@ptd.net>
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
# 20200604  Geert van Boxtel              First version for v0.1.0
#------------------------------------------------------------------------------

#' IIR lowpass filter to IIR multiband
#'
#' Transform an IIR lowpass filter prototype to an IIR multiband filter.
#'
#' The utility of a prototype filter comes from the property that all other
#' filters can be derived from it by applying a scaling factor to the components
#' of the prototype. The filter design need thus only be carried out once in
#' full, with other filters being obtained by simply applying a scaling factor.
#' Especially useful is the ability to transform from one bandform to another.
#' In this case, the transform is more than a simple scale factor. Bandform here
#' is meant to indicate the category of passband that the filter possesses. The
#' usual bandforms are lowpass, highpass, bandpass and bandstop, but others are
#' possible. In particular, it is possible for a filter to have multiple
#' passbands. In fact, in some treatments, the bandstop filter is considered to
#' be a type of multiple passband filter having two passbands. Most commonly,
#' the prototype filter is expressed as a lowpass filter, but other techniques
#' are possible[1].
#'
#' Filters with multiple passbands may be obtained by applying the general
#' transformation described in [2].
#'
#' Because \code{iirlp2mb} is generic, it can be extended to accept other
#' inputs.
#'
#' @param b numerator polynomial of prototype low pass filter
#' @param a denominator polynomial of prototype low pass filter
#' @param Wo (normalized angular frequency)/pi to be transformed
#' @param Wt vector of (norm. angular frequency)/pi transform targets
#' @param type one of "pass" or "stop". Specifies to filter to produce: bandpass
#'   (default) or bandstop.
#' @param ... additional arguments (not used)
#'
#' @return List of class \code{\link{Arma}} numerator and denominator
#'   polynomials of the resulting filter.
#'
#' @examples
#' ## Design a prototype real IIR lowpass elliptic filter with a gain of about
#' ## –3 dB at 0.5pi rad/sample.
#' el <- ellip(3, 0.1, 30, 0.409)
#' ## Create a real multiband filter with two passbands.
#' mb1 <- iirlp2mb(el, 0.5, c(.2, .4, .6, .8), 'pass')
#' ## Create a real multiband filter with two stopbands.
#' mb2 <- iirlp2mb(el, 0.5, c(.2, .4, .6, .8), 'stop')
#' ## Compare the magnitude responses of the filters.
#' hfl <- freqz(el)
#' hf1 <- freqz(mb1)
#' hf2 <- freqz(mb2)
#' plot(hfl$w, 20 * log10(abs(hfl$h)), type = "l",
#'     xlab = "Normalized frequency (* pi rad/sample)",
#'     ylab = "Magnitude (dB)")
#' lines(hf1$w, 20 * log10(abs(hf1$h)), col="red")
#' lines(hf2$w, 20 * log10(abs(hf2$h)), col="blue")
#' legend('bottomleft',
#'        legend = c('Prototype', 'Two passbands', 'Two Stopbands'),
#'        col=c("black", "red", "blue"), lty = 1)
#'
#' @references [1] \url{https://en.wikipedia.org/wiki/Prototype_filter}\cr
#' [2] \url{https://en.wikipedia.org/wiki/Prototype_filter#Lowpass_to_multi-band}
#'
#' @author Alan J. Greenberger, \email{alanjg@@ptd.net}.\cr
#'   Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#'
#' @rdname iirlp2mb
#' @export

iirlp2mb <- function(b, ...) UseMethod("iirlp2mb")

#' @rdname iirlp2mb
#' @export

iirlp2mb.Arma <- function(b, Wo, Wt, type, ...) {
  iirlp2mb(b$b, b$a, Wo, Wt, type, ...)
}

#' @rdname iirlp2mb
#' @export

iirlp2mb.Zpg <- function(b, Wo, Wt, type, ...) {
  ba <- as.Arma(b)
  iirlp2mb(ba$b, ba$a, Wo, Wt, type, ...)
}

#' @rdname iirlp2mb
#' @export

iirlp2mb.Sos <- function(b, Wo, Wt, type, ...) {
  ba <- as.Arma(b)
  iirlp2mb(ba$b, ba$a, Wo, Wt, type, ...)
}

#' @rdname iirlp2mb
#' @export

iirlp2mb.default <- function(b, a, Wo, Wt, type = c("pass", "stop"), ...) {

  # input validation
  type <- match.arg(type)
  if (type == "pass") {
    pass_stop <- -1
  } else if (type == "stop") {
    pass_stop <- 1
  }
  if (!isPosscal(Wo) || Wo > 1) {
    stop(paste("Frequency value Wo of prototype filter",
               "must be a scalar between 0 and 1"))
  }
  if (any(Wt < 0) || any(Wt > 1)) {
    stop("Frequency values Wt of target filter must be between 0 and 1")
  }
  Wt <- unique(sort(Wt))

  ## The first stage allpass denominator computation
  K <- apd(pi * Wo)

  ## The second stage allpass computation
  phi <- pi * Wt
  P <- apd(phi)
  PP <- rev(P)

  AllpassDen <- P - (K[2] * PP)
  AllpassDen <- AllpassDen / AllpassDen[1]      # normalize
  AllpassNum <- pass_stop * rev(AllpassDen)
  ba <- transform(b, a, AllpassNum, AllpassDen, pass_stop)
  ba
}
###############################################################################
# Helper functions for iirlp2mb, not exported from the namespace

# all pass denominator
apd <- function(phi) {
  Pkm1 <- 1                 # P0 initial condition from [FFM] eq. 22
  for (k in seq_along(phi)) {
    P <- pk(Pkm1, k, phi[k])
    Pkm1 <- P
  }
  P
}

# kth iteration of P(z)
pk <- function(Pkm1, k, phik) {

  Pk <- rep(0L, k + 1)
  sin_k <- sin(phik / 2)
  cos_k <- cos(phik / 2)
  for (i in 1:k) {
    Pk[i] <- Pk[i] +  sin_k * Pkm1[i] - ((-1)^k * cos_k * Pkm1[k + 1 - i])
    Pk[i + 1] <- Pk[i + 1] +
      sin_k * Pkm1[i] +
      ((-1)^k * cos_k * Pkm1[k + 1 - i])
  }
  Pk <- Pk / Pk[1]
  Pk
}

# Regenerate ith power of P from stored PPower
ppower <- function(Ppower, i, powcols) {

  if (i == 0) {
    p  <- 1
  } else {
    p  <- NULL
    for (j in 1:powcols) {
      if (is.na(Ppower[i, j])) break
      p <-  cbind(p, Ppower[i, j])
    }
  }
  p
}

# add polynomials of possibly different length
polysum <- function(p1, p2) {

  n1 <- length(p1)
  n2 <- length(p2)
  if (n1 > n2) {
    ## pad p2
    p2 <- c(p2, rep(0L, n1 - n2))
  } else if (n2 > n1) {
    ## pad p1
    p1 <- c(p1, rep(0L, n2 - n1))
  }
  poly <- p1 + p2
  poly
}

transform <- function(B, A, PP, P, pass_stop) {

  na <- length(A)
  nb <- length(B)
  n  <- max(na, nb)
  np <- length(P)
  powcols <- np + (np - 1) * (n - 2)
  Ppower <- matrix(NA, nrow = n - 1, ncol = powcols)
  Ptemp <- P
  for (i in 1:(n - 1)) {
    for (j in seq_along(Ptemp)) {
      Ppower[i, j]  <- Ptemp[j]
    }
    Ptemp <- conv(Ptemp, P)
  }

  ## Compute numerator and denominator of transformed filter
  Num <- Den <- NULL
  for (i in 1:n) {
    if ((n - i) == 0) {
      p_pownmi <- 1
    } else {
      p_pownmi <- ppower(Ppower, n - i, powcols)
    }
    if (i == 1) {
      pp_powim1 <- 1
    } else {
      pp_powim1 <- rev(ppower(Ppower, i - 1, powcols))
    }
    if (i <= nb) {
      Bterm <- (pass_stop ^ (i - 1)) * B[i] * conv(pp_powim1, p_pownmi)
      Num <- polysum(Num, Bterm)
    }
    if (i <= na) {
      Aterm <- (pass_stop ^ (i - 1)) * A[i] * conv(pp_powim1, p_pownmi)
      Den <- polysum(Den, Aterm)
    }
  }
  ## Scale both numerator and denominator to have Den(1) = 1
  Den <- Den / Den[1]
  Num <- Num / Den[1]

  Arma(Num, Den)
}
