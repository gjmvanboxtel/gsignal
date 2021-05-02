# freqz.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Original Octave function:
# Copyright (C) 1994-2017 John W. Eaton
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
# 20200422  GvB       setup for gsignal v0.1.0
# 20200423  GvB       corrected minor bug in print.summary.freqz
#                     print phase in degrees
# 20200425  GvB       Added S3 method for class 'Zpg'
# 20200515  GvB       resolve infinite ylim values in freqz.plot
# 20200616  GvB       pass default parameters to methods
# 20200629  GvB       bug in parameter passing for class 'Ma'
# 20201103  GvB       changed the S3 method handling
#------------------------------------------------------------------------------

#' Frequency response of digital filter
#'
#' Compute the z-plane frequency response of an ARMA model or rational IIR
#' filter.
#'
#' The frequency response of a digital filter can be interpreted as the transfer
#' function evaluated at \eqn{z = e^{j\omega}}.
#'
#' The 'Matlab' and 'Octave' versions of \code{freqz} produce magnitude and
#' phase plots. The \code{freqz} version in the 'signal' package produces
#' separate plots of magnitude in the pass band (max - 3 dB to max) and stop
#' (total) bands, as well as a phase plot. The current version produces slightly
#' different plots. The magnitude plots are separate for stop and pass bands,
#' but the pass band plot has an absolute lower limit of -3 dB instead of max -
#' 3 dB. In addition a \code{summary} method was added that prints out the most
#' important information about the frequency response of the filter.
#'
#' @note When results of \code{freqz} are printed, \code{freqz_plot} will be
#'   called to display frequency plots of magnitude and phase. As with lattice
#'   plots, automatic printing does not work inside loops and function calls, so
#'   explicit calls to print or plot are needed there.
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
#' @param h complex frequency response \eqn{H(e^{j\omega})}, specified as a
#'   vector.
#' @param ...	 for methods of \code{freqz}, arguments are passed to the default
#'   method. For \code{freqz_plot}, additional arguments are passed through to
#'   plot.
#' @param object object of class \code{"freqz"} for \code{summary}
#'
#' @return For \code{freqz}, a list of class \code{'freqz'} with items:
#' \describe{
#'   \item{h}{complex array of frequency responses at frequencies \code{f}.}
#'   \item{w}{array of frequencies.}
#'   \item{u}{units of (angular) frequency; either rad/s or Hz.}
#' }
#'
#' @examples
#' b <- c(1, 0, -1)
#' a <- c(1, 0, 0, 0, 0.25)
#' freqz(b, a)
#'
#' hw <- freqz(b, a)
#' summary(hw)
#'
#' @author John W. Eaton, Paul Kienzle, \email{pkienzle@@users.sf.net}.\cr
#'  Port to R by Tom Short,\cr
#'  adapted by Geert van Boxtel, \email{gjmvanboxtel@@gmail.com}
#'
#' @rdname freqz
#' @export

freqz <- function(filt, ...) UseMethod("freqz")

#' @rdname freqz
#' @export
freqz.default <- function(filt, a = 1, n = 512,
                          whole = ifelse((is.numeric(filt) && is.numeric(a)),
                                         FALSE, TRUE),
                          fs = 2 * pi, ...)  {

  if (!(is.vector(filt) && is.vector(a))) {
    stop("'filt' and 'a' must be vectors")
  }
  b <- filt
  if (!is.logical(whole)) {
    whole <- FALSE
  }
  if (fs == 2 * pi) {
    u <- "rad/s"
  } else {
    u <- "Hz"
  }

  if (length(n) > 1) { ## Explicit frequency vector given
    w <- 2 * pi * n / fs
    hb <- pracma::polyval(rev(b), exp(-1i * w))
    ha <- pracma::polyval(rev(a), exp(-1i * w))
    w <- n
  } else if (whole) {
    w <- fs * (0:(n - 1)) / n
    hb <- stats::fft(postpad(b, n))
    ha <- stats::fft(postpad(a, n))
  } else {
    w <- fs / 2 * (0:(n - 1)) / n
    hb <- stats::fft(postpad(b, 2 * n))[1:n]
    ha <- stats::fft(postpad(a, 2 * n))[1:n]
  }

  h <- hb / ha

  res <- list(h = h, w = w, u = u)
  class(res) <- "freqz"
  res
}

#' @rdname freqz
#' @export

freqz.Arma <- function(filt, n = 512,
                       whole = ifelse(
                         (is.numeric(filt$b) && is.numeric(filt$a)),
                         FALSE, TRUE),
                       fs = 2 * pi, ...)
  freqz.default(filt$b, filt$a, n, whole, fs, ...)

#' @rdname freqz
#' @export

freqz.Ma <- function(filt, n = 512,
                     whole = ifelse(is.numeric(filt), FALSE, TRUE),
                     fs = 2 * pi, ...) # FIR
  freqz.default(unclass(filt), 1, n, whole, fs, ...)

#' @rdname freqz
#' @export

freqz.Sos <- function(filt, n = 512, whole = FALSE, fs = 2 * pi, ...)
  freqz.Arma(as.Arma(filt), n, whole, fs, ...)

#' @rdname freqz
#' @export

freqz.Zpg <- function(filt, n = 512, whole = FALSE, fs = 2 * pi, ...)
  freqz.Arma(as.Arma(filt), n, whole, fs, ...)

#' @rdname freqz
#' @export

print.freqz <- plot.freqz <- function(x, ...)
  freqz_plot(x$w, x$h, ...)

#' @rdname freqz
#' @export

summary.freqz <- function(object, ...) {

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
                 cutoff = cutoff, rp = rp, u = object$u),
            class = c("summary.freqz", "list"))
}

#' @rdname freqz
#' @export

print.summary.freqz <- function(x, ...) {

  cat(paste0("\nSummary of freqz object '", x$nm, "':\n"))
  rw <- round(x$rw, 3)
  cat(paste("\nFrequencies ranging from", rw[1], "to", rw[2], x$u))
  mmag <- round(x$mmag, 3)
  wmag <- round(x$wmag, 3)
  cat(paste0("\nMaximum magnitude ", mmag, " dB at frequency ", wmag, " ", x$u))
  cutoff <- round(x$cutoff, 3)
  lc <- length(cutoff)
  fr <- ifelse(lc > 1, "frequencies", "frequency")
  cat(paste0("\n-3 dB cutoff at ", fr, " ", cutoff[1]))
  if (lc > 1) {
    for (i in 2:lc) {
      cat(paste(",", cutoff[i]))
    }
  }
  cat(paste0(" ", x$u))
  rp <- round(x$rp, 3)
  rpd <- round(rp * 360 / (2 * pi), 3)
  cat(paste0("\nPhase ranging from ", rp[1], " to ", rp[2], " rad (",
             rpd[1], " to ", rpd[2], " degrees)"))
  cat("\n")
}

#' @rdname freqz
#' @export

freqz_plot <- function(w, h, ...) {

  mag <- 20 * log10(abs(h))
  maxmag <- max(mag, na.rm = TRUE)
  if (is.na(maxmag) || maxmag == Inf) maxmag <- 1
  argh <- Arg(h)
  argh[which(is.na(argh))] <- 0
  phase <- unwrap(argh)
  op <- graphics::par(mfrow = c(3, 1), mar = c(4, 4, 1.5, 1))
  on.exit(graphics::par(op))
  graphics::plot(w, mag, type = "l", xlab = "", ylab = "",
                 ylim = c(-4, maxmag), ...)
  graphics::title("Pass band (dB)")
  graphics::abline(h = -3, col = "red", lty = 2)
  graphics::plot(w, mag, type = "l", xlab = "", ylab = "", ...)
  graphics::title("Stop band (dB)")
  graphics::abline(h = -3, col = "red", lty = 2)
  graphics::plot(w, phase * 360 / (2 * pi), type = "l",
                 xlab = "Frequency", ylab = "", ...)
  graphics::title("Phase (degrees)")
}
