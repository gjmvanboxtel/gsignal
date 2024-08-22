# freqs.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Original Octave function:
# Copyright (C) 2003 Julius O. Smith III <jos@ccrma.stanford.edu>
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
# 20200413  GvB       setup for gsignal v0.1.0
# 20200427  GvB       added S3 methods
# 20240822  GvB       changed setup to match freqz using S3 methods
#------------------------------------------------------------------------------

#' Frequency response of analog filters
#'
#' Compute the s-plane frequency response of an IIR filter.
#'
#' The s-plane frequency response of the IIR filter \code{B(s) / A(s)} is
#' computed as \code{H = polyval(B, 1i * W) / polyval(A, 1i * W)}. If called
#' with no output argument, a plot of magnitude and phase are displayed.
#'
#' @param filt for the default case, moving average (MA) polynomial
#'   coefficients, specified as a numeric vector or matrix. In case of a matrix,
#'   then each row corresponds to an output of the system. The number of columns
#'   of \code{b} must be less than or equal to the length of \code{a}.
#' @param a autoregressive (AR) polynomial coefficients, specified as a vector.
#' @param w angular frequencies, specified as a positive real vector expressed
#'   in rad/second.
#' @param x	object to be printed or plotted.
#' @param object object of class \code{"freqs"} for \code{summary}
#' @param ... for methods of \code{freqs}, arguments are passed to the default
#'   method. For \code{freqs_plot}, additional arguments are passed through to
#'   plot.
#'
#' @return For \code{freqs}, a list of class \code{'freqs'} with items:
#' \describe{
#'   \item{h}{complex array of frequency responses at frequencies \code{f}.}
#'   \item{w}{array of frequencies.}
#' }
#'
#' @examples
#' b <- c(1, 2); a <- c(1, 1)
#' w <- seq(0, 4, length.out = 128)
#' freqs (b, a, w)
#'
#' @author Julius O. Smith III, \email{jos@@ccrma.stanford.edu}.\cr
#' Conversion to R by Geert van Boxtel \email{gjmvanboxtel@@gmail.com}
#'
#' @rdname freqs
#' @export

freqs <- function(filt, ...) UseMethod("freqs")

#' @rdname freqs
#' @export

freqs.default <- function(filt, a, w, ...) {

  h <- pracma::polyval(filt, 1i * w) / pracma::polyval(a, 1i * w)

  res <- list(h = h, w = w)
  class(res) <- "freqs"
  res
}

#' @rdname freqs
#' @export

freqs.Arma <- function(filt, w, ...) # IIR
  freqs.default(filt$b, filt$a, w, ...)

#' @rdname freqs
#' @export

freqs.Ma <- function(filt, w, ...) # FIR
  freqs.default(filt, 1, w, ...)

#' @rdname freqs
#' @export

freqs.Sos <- function(filt, w, ...) # second-order sections
  freqs.Arma(as.Arma(filt), w, ...)

#' @rdname freqs
#' @export

freqs.Zpg <- function(filt, w, ...) # zero-pole-gain
  freqs.Arma(as.Arma(filt), w, ...)

#' @rdname freqs
#' @export

print.freqs <- plot.freqs <- function(x, ...)
  freqs_plot(x, ...)

#' @rdname freqs
#' @export

summary.freqs <- function(object, ...) {
  
  nm <- deparse(substitute(object))
  h <- object$h
  w <- object$w
  
  rw <- range(w)
  
  mag <- 20 * log10(abs(h))
  mmag <- max(mag)
  wmag <- w[which.max(mag)]
  cutoff <-
    w[diff(ifelse((!is.finite(mag) | is.na(mag) | mag < -3), 0, 1)) != 0]
  
  phase <- unwrap(Arg(h))
  rp <- range(phase)
  
  structure(list(nm = nm, rw = rw, mmag = mmag, wmag = wmag,
                 cutoff = cutoff, rp = rp),
            class = c("summary.freqs", "list"))
}

#' @rdname freqs
#' @export

print.summary.freqs <- function(x, ...) {
  
  cat(paste0("\nSummary of freqs object '", x$nm, "':\n"))
  rw <- round(x$rw, 3)
  cat(paste("\nFrequencies ranging from", rw[1], "to", rw[2]))
  mmag <- round(x$mmag, 3)
  wmag <- round(x$wmag, 3)
  cat(paste0("\nMaximum magnitude ", mmag, " dB at frequency ", wmag))
  cutoff <- round(x$cutoff, 3)
  lc <- length(cutoff)
  fr <- ifelse(lc > 1, "frequencies", "frequency")
  cat(paste0("\n-3 dB cutoff at ", fr, " ", cutoff[1]))
  if (lc > 1) {
    for (i in 2:lc) {
      cat(paste(",", cutoff[i]))
    }
  }
  rp <- round(x$rp, 3)
  rpd <- round(rp * 360 / (2 * pi), 3)
  cat(paste0("\nPhase ranging from ", rp[1], " to ", rp[2], " rad (",
             rpd[1], " to ", rpd[2], " degrees)"))
  cat("\n")
}

#' @rdname freqs
#' @export

freqs_plot <- function(x, ...) {
  
  mag <- 20 * log10(abs(x$h))
  phase <- unwrap(Arg(x$h))
  
  op <- graphics::par(mfrow = c(2, 1))
  on.exit(graphics::par(op))
  
  graphics::plot(x$w, mag, type = "l", xlab = "", ylab = "dB", ...)
  graphics::legend("topright", "Magnitude (dB)", lty = 1)
  graphics::title("Frequency response plot by freqs")
  
  graphics::plot(x$w, phase / (2 * pi), type = "l",
                 xlab = "Frequency (rad/s)", ylab = "Phase", ...)
  graphics::legend("topright", "Phase (radians / 2 pi)", lty = 1)
  graphics::title("")
  
}
