# phasez.R
# Copyright (C) 2026 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Original Octave function:
# Copyright (C) 2023 Leonardo Araujo <leolca@gmail.com>
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
# 20261009  GvB       setup for gsignal v0.4.0
#------------------------------------------------------------------------------

#' Phase response of digital filter
#'
#' Compute the z-plane frequency response of an ARMA model or rational IIR
#' filter.
#'
#' @param filt for the default case, the moving-average coefficients of an ARMA
#'   model or filter. Generically, \code{filt} specifies an arbitrary model or
#'   filter operation.
#' @param a the autoregressive (recursive) coefficients of an ARMA filter.
#' @param n	number of points at which to evaluate the frequency response. If
#'   \code{n} is a vector with a length greater than 1, then evaluate the
#'   frequency response at these points. For fastest computation, \code{n}
#'   should factor into a small number of small primes. Default: 512.
#' @param whole	FALSE (the default) to evaluate around the upper half of the
#'   unit circle or TRUE to evaluate around the entire unit circle.
#' @param fs sampling frequency in Hz. If not specified (default = 2 * pi), the
#'   frequencies are in radians.
#' @param x	object to be printed or plotted.
#' @param w vector of frequencies
#' @param phi phase response, specified as a vector.
#' @param ...	 for methods of \code{phasez}, arguments are passed to the default
#'   method. For \code{phasez_plot}, additional arguments are passed through to
#'   plot.
#' @param object object of class \code{"phasez"} for \code{summary}
#'
#' @return For \code{phasez}, a list of class \code{'phasez'} with items:
#' \describe{
#'   \item{phi}{complex array of phase responses at frequencies \code{w}.}
#'   \item{w}{array of frequencies.}
#'   \item{u}{units of (angular) frequency; either rad/s or Hz.}
#' }
#'
#' @examples
#' b <- c(1, 0, -1)
#' a <- c(1, 0, 0, 0, 0.25)
#' phasez(b, a)
#'
#' ph <- phasez(b, a)
#' summary(ph)
#' 
#' sos <- cheby2(40, 50, 0.4, output = "Sos")
#' ph <- phasez(sos)
#' summary(ph)
#' ph
#'
#' @author Leonardo Araujo \email{leolca@@gmail.com}.\cr
#'  Port to R by Geert van Boxtel, \email{gjmvanboxtel@@gmail.com}
#'  
#' @references Oppenheim, A., and Schafer, R. Discrete-Time Signal
#'   Processing. 3rd edition, Pearson, 2009.
#'
#' @rdname phasez
#' @export

phasez <- function(filt, ...) UseMethod("phasez")

#' @rdname phasez
#' @export
phasez.default <- function(filt, a = 1, n = 512,
                          whole = ifelse((is.numeric(filt) && is.numeric(a)),
                                         FALSE, TRUE),
                          fs = 2 * pi, ...)  {

  if (!(is.vector(filt) && is.vector(a))) {
    stop("filt and a must be vectors")
  }
  if (anyNA(filt) || anyNA(a)) {
    stop("filt and a must not contain missing values")
  } else {
    b <- filt
  }
  if(!is.vector(n) || anyNA(n) || any(n < 0)) {
    stop("n must be positive")
  }
  if (!is.logical(whole)) {
    whole <- FALSE
  }
  if (!(isScalar(fs) && fs > 0)) {
    stop("fs must be a scalar > 0")
  }
  if (fs == 2 * pi) {
    u <- "rad/s"
  } else {
    u <- "Hz"
  }


  hw <- gsignal::freqz(filt, a, n, whole, fs, ...)
  phi <- Arg(hw$h)
  phi[which(is.na(phi))] <- 0
  phi <- my_unwrap(phi)
  res <- list(phi = phi, w = hw$w, u = hw$u)
  class(res) <- "phasez"
  res
}

#' @rdname phasez
#' @export

phasez.Arma <- function(filt, ...)
  phasez.default(filt$b, filt$a, ...)

#' @rdname phasez
#' @export

phasez.Ma <- function(filt, ...) # FIR
  phasez.default(unclass(filt), ...)

#' @rdname phasez
#' @export

phasez.Sos <- function(filt, ...) {
  hw <- freqz.Sos(filt)
  phi <- Arg(hw$h)
  phi[which(is.na(phi))] <- 0
  phi <- my_unwrap(phi)
  res <- list(phi = phi, w = hw$w, u = hw$u)
  class(res) <- "phasez"
  res
}

#' @rdname phasez
#' @export

phasez.Zpg <- function(filt, n = 512, whole = FALSE, fs = 2 * pi, ...)
  phasez.Sos(as.Sos(filt), ...)

#' @rdname phasez
#' @export

print.phasez <- plot.phasez <- function(x, ...)
  phasez_plot(x$w, x$phi, ...)

#' @rdname phasez
#' @export

summary.phasez <- function(object, ...) {

  nm <- deparse(substitute(object))
  phi <- object$phi
  w <- object$w

  rw <- range(w)

  phi[which(is.na(phi))] <- 0
  phase <- my_unwrap(phi)
  maxph <- max(phase)
  wmaxph <- w[which.max(phase)]
  minph <- min(phase)
  wminph <- w[which.min(phase)]
  
  rp <- range(phi)

  structure(list(nm = nm, rw = rw, maxph = maxph, wmaxph = wmaxph,
                 minph = minph, wminph = wminph, rp = rp, u = object$u),
            class = c("summary.phasez", "list"))
}

#' @rdname phasez
#' @export

print.summary.phasez <- function(x, ...) {

  cat(paste0("\nSummary of phasez object '", x$nm, "':\n"))
  rw <- round(x$rw, 3)
  cat(paste("\nFrequencies ranging from", rw[1], "to", rw[2], x$u))
  maxph <- round(x$maxph, 3)
  wmaxph <- round(x$wmaxph, 3)
  cat(paste0("\nMaximum phase ", maxph, " rad (",
             round(maxph * 360 / (2 * pi), 3), " degrees) at frequency ",
             wmaxph, " ", x$u))
  minph <- round(x$minph, 3)
  wminph <- round(x$wminph, 3)
  cat(paste0("\nMinimum phase ", minph, " rad (",
             round(minph * 360 / (2 * pi), 3), " degrees) at frequency ",
             wminph, " ", x$u))
  rp <- round(x$rp, 3)
  rpd <- round(rp * 360 / (2 * pi), 3)
  cat(paste0("\nPhase ranging from ", rp[1], " to ", rp[2], " rad (",
             rpd[1], " to ", rpd[2], " degrees)"))
  cat("\n\n")
}

#' @rdname phasez
#' @export

phasez_plot <- function(w, phi, ...) {

  graphics::plot(w, phi * 360 / (2 * pi), type = "l",
                 xlab = "Frequency", ylab = "", ...)
  graphics::title("Phase (degrees)")
}

my_unwrap <- function(x) {
  lx <- length(x)
  stillunwrap <- TRUE
  while (stillunwrap) {
    dx <- diff(x)
    idx <- which(abs(dx) > pi - 0.05)[1]
    if (!is.na(idx)) {
      if (dx[idx] > 0) {
        x[(idx + 1):lx] <- x[(idx + 1):lx] - pi
      } else {
        x[(idx + 1):lx] <- x[(idx + 1):lx] + pi
      }
    } else {
      stillunwrap <- FALSE
    }
  }
  x
}
