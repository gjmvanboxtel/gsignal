# zp2sos.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Original Octave version:
# Copyright (C) 2005 Julius O. Smith III <jos@ccrma.stanford.edu>
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
# 20200331  GvB       setup for gsignal v0.1.0
# 20200401  GvB       handle NULL input
# 20200403  GvB       use 'relaxed' tolerance 1e-7 for cplxreal,
#                     flipud sos for compatibility with Matlab/Octave
# 20200406  GvB       validated
# 20210326  GvB       renamed k to g; added 'order' argument;
#                     return class 'Sos' argument
#------------------------------------------------------------------------------

#' Zero-pole-gain to second-order section format
#'
#' Convert digital filter zero-pole-gain data to second-order section form.
#'
#' @param z complex vector of the zeros of the model (roots of \code{B(z)})
#' @param p complex vector of the poles of the model (roots of \code{A(z)})
#' @param g overall gain (\code{B(Inf)}). Default: 1
#' @param order row order, specified as:
#' \describe{
#'   \item{"up"}{order the sections so the first row contains the poles farthest
#'   from the unit circle.}
#'   \item{"down" (Default)}{order the sections so the first row of \code{sos}
#'   contains the poles closest to the unit circle.}
#' }
#' The ordering influences round-off noise and the probability of overflow.
#'
#' @return A list with the following list elements:
#' \describe{
#'   \item{sos}{Second-order section representation, specified as an nrow-by-6
#'   matrix, whose rows contain the numerator and denominator coefficients of
#'   the second-order sections:\cr \code{sos <- rbind(cbind(B1, A1), cbind(...),
#'   cbind(Bn, An))}, where \code{B1 <- c(b0, b1, b2)}, and \code{A1 <- c(a0,
#'   a1, a2)} for section 1, etc. The b0 entry must be nonzero for each
#'   section.}
#'   \item{g}{Overall gain factor that effectively scales the output \code{b}
#'   vector (or any one of the input \code{Bi} vectors).}
#' }
#'
#' @seealso \code{\link{as.Sos}}, \code{\link{filter}}, \code{\link{sosfilt}}
#'
#' @examples
#' zpk <- tf2zp (c(1, 0, 0, 0, 0, 1), c(1, 0, 0, 0, 0, .9))
#' sosg <- zp2sos (zpk$z, zpk$p, zpk$g)
#'
#' @author Julius O. Smith III, \email{jos@@ccrma.stanford.edu}.\cr
#' Conversion to R by Geert van Boxtel, \email{gjmvanboxtel@@gmail.com}
#'
#' @export

zp2sos <- function(z, p, g = 1, order = c("down", "up")) {

  order <- match.arg(order)

  if (!is.null(z)) {
    zcr <- cplxreal(z, tol = 1e-7)
    if (is.null(zcr$zc)) {
      zc <- nzc <- 0
    } else {
      zc <- zcr$zc
      nzc <- length(zc)
    }
    if (is.null(zcr$zr)) {
      zr <- nzr <- 0
    } else {
      zr <- zcr$zr
      nzr <- length(zr)
    }
  } else {
    zc <- zr <- 0
    nzc <- nzr <- 0
  }

  if (!is.null(p)) {
    pcr <- cplxreal(p, tol = 1e-7)
    if (is.null(pcr$zc)) {
      pc <- npc <- 0
    } else {
      pc <- pcr$zc
      npc <- length(pc)
    }
    if (is.null(pcr$zr)) {
      pr <- npr <- 0
    } else {
      pr <- pcr$zr
      npr <- length(pr)
    }
  } else {
    pc <- pr <- 0
    npc <- npr <- 0
  }

  # Pair up real zeros
  if (nzr > 0) {
    if (nzr %% 2 == 1) {
      zr <- c(zr, 0)
      nzr <- nzr + 1
    }
    nzrsec <- nzr / 2
    zrms <- -zr[seq(1, nzr - 1, 2)] - zr[seq(2, nzr, 2)]
    zrp <- zr[seq(1, nzr - 1, 2)] * zr[seq(2, nzr, 2)]
  } else {
    nzrsec <- 0
  }

  # Pair up real poles:
  if (npr > 0) {
    if (npr %% 2 == 1) {
      pr <- c(pr, 0)
      npr <- npr + 1
    }
    nprsec <- npr / 2
    prms <- -pr[seq(1, npr - 1, 2)] - pr[seq(2, npr, 2)]
    prp <- pr[seq(1, npr - 1, 2)] * pr[seq(2, npr, 2)]
  } else {
    nprsec <- 0
  }

  nsecs <- max(nzc + nzrsec, npc + nprsec)
  if (nsecs <= 0) nsecs <- 1

  # Convert complex zeros and poles to real 2nd-order section form:
  zcm2r <- -2 * Re(zc)
  zca2 <- abs(zc)^2
  pcm2r <- -2 * Re(pc)
  pca2 <- abs(pc)^2

  sos <- matrix(0L, nsecs, 6)
  sos[, 1] <- rep(1L, nsecs)    # all 2nd-order polynomials are monic
  sos[, 4] <- rep(1L, nsecs)

  nzrl <- nzc + nzrsec    # index of last real zero section
  nprl <- npc + nprsec    # index of last real pole section

  for (i in seq_len(nsecs)) {

    if (i <= nzc) {            # lay down a complex zero pair:
      sos[i, 2:3] <- c(zcm2r[i], zca2[i])
    } else if (i <= nzrl) {    # lay down a pair of real zeros:
      sos[i, 2:3] <- c(zrms[i - nzc], zrp[i - nzc])
    }

    if (i <= npc) {              # lay down a complex pole pair:
      sos[i, 5:6] <- c(pcm2r[i], pca2[i])
    } else if (i <= nprl) {    # lay down a pair of real poles:
      sos[i, 5:6] <- c(prms[i - npc], prp[i - npc])
    }
  }

  if (order == "down") {
    rv <- Sos(sos = sos, g = g)
  } else {
    rv <- Sos(sos = sos[rev(seq_len(nsecs)), ], g = g)
  }
  rv
}
