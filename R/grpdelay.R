# grpdelay.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Original Octave function:
# Copyright (C) 2000 Paul Kienzle <pkienzle@users.sf.net>
# Copyright (C) 2004 Julius O. Smith III <jos@ccrma.stanford.edu>
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
#------------------------------------------------------------------------------

#' Group delay
#'
#' Compute the average delay of a filter (group delay).
#'
#' If the denominator of the computation becomes too small, the group delay is
#' set to zero. (The group delay approaches infinity when there are poles or
#' zeros very close to the unit circle in the z plane.)
#'
#' @param filt for the default case, the moving-average coefficients of an ARMA
#'   model or filter. Generically, filt specifies an arbitrary model or filter
#'   operation.
#' @param a the autoregressive (recursive) coefficients of an ARMA filter.
#' @param n	number of points at which to evaluate the frequency response. If
#'   \code{n} is a vector with a length greater than 1, then evaluate the
#'   frequency response at these points. For fastest computation, \code{n}
#'   should factor into a small number of small primes. Default: 512.
#' @param whole	FALSE (the default) to evaluate around the upper half of the
#'   unit circle or TRUE to evaluate around the entire unit circle.
#' @param fs sampling frequency in Hz. If not specified, the frequencies are in
#'   radians.
#' @param x	object to be plotted.
#' @param xlab,ylab,type as in plot, but with more sensible defaults.
#' @param ...	 for methods of grpdelay, arguments are passed to the default
#'   method. For plot.grpdelay, additional arguments are passed through to plot.
#'
#' @return A list of class \code{grpdelay} with items:
#' \describe{
#'   \item{gd}{the group delay, in units of samples. It can be converted to
#'   seconds by multiplying by the sampling period (or dividing by the sampling
#'   rate fs).}
#'   \item{w}{frequencies at which the group delay was calculated.}
#'   \item{ns}{number of points at which the group delay was calculated.}
#'   \item{Hzflag}{TRUE for frequencies in Hz, FALSE for frequencies in
#'   radians.}
#'
#' }
#'
#' @examples
#' # Two Zeros and Two Poles
#' b <- poly(c(1 / 0.9 * exp(1i * pi * 0.2), 0.9 * exp(1i * pi * 0.6)))
#' a <- poly(c(0.9 * exp(-1i * pi * 0.6), 1 / 0.9 * exp(-1i * pi * 0.2)))
#' gpd <- grpdelay(b, a, 512, whole = TRUE, fs = 1)
#' print(gpd)
#' plot(gpd)
#'
#' @author Paul Kienzle, \email{pkienzle@@users.sf.net},\cr
#'  Julius O. Smith III, \email{jos@@ccrma.stanford.edu}.\cr
#'  Conversion to R by Tom Short,\cr
#'  adapted by Geert van Boxtel, \email{gjmvanboxtel@@gmail.com}
#'
#' @references
#' \url{https://ccrma.stanford.edu/~jos/filters/Numerical_Computation_Group_Delay.html}\cr
#' \url{https://en.wikipedia.org/wiki/Group_delay}
#'
#' @rdname grpdelay
#' @export

grpdelay <- function(filt, ...) UseMethod("grpdelay")

#' @rdname grpdelay
#' @export

print.grpdelay <- function(x, ...) {

  cat("- Group delay (gd) calculated at", x$ns, "points.\n")
  cat("- Frequencies (w) given in",
      if (x$HzFlag) "*Hz*." else "*radians*.", "\n")
  temp <- data.frame(do.call("cbind", x[c("gd", "w")]))
  if (nrow(temp) > 8L) {
    print(utils::head(temp, n = 4L), row.names = FALSE, ...)
    cat(" .......  .......\n")
  } else print(temp, row.names = FALSE, ...)
  invisible(x)
}

#' @rdname grpdelay
#' @export

plot.grpdelay <- function(x,
                          xlab = if (x$HzFlag) "Frequency (Hz)"
                          else "Frequency (rad/sample)",
                          ylab = "Group delay (samples)",
                          type = "l", ...) {
  graphics::plot(x$w[1:x$ns], x$gd[1:x$ns],
       xlab = xlab, ylab = ylab, type = type, ...)
}

#' @rdname grpdelay
#' @export

grpdelay.default <- function(filt, a = 1, n = 512,
                             whole = FALSE, fs = NULL, ...)   {

  b <- as.vector(filt)
  a <- as.vector(a)
  n <- as.vector(n)

  if (whole == "whole" || whole) {
    whole <- TRUE
  } else {
    whole <- FALSE
  }

  if (is.null(fs)) {
    HzFlag <- FALSE
    fs <- 1
  } else {
    HzFlag <- TRUE
  }

  if (length(n) == 1) {
    nfft <- n
    if (!whole) {
      nfft <- 2 * nfft
    }
    w <- fs * (0:(nfft - 1)) / nfft
    if (!HzFlag) {
      w <- w * 2 * pi
    }
  } else if (length(n) > 1) {
    w <- n
    nfft <- length(w) * 2
    whole <- FALSE
  } else {
    stop("n must be a vector with a length >= 1")
  }

  oa <- length(a) - 1             # order of a(z)
  if (oa < 0) {
    a <- 1
    oa <- 0
  }
  ob <- length(b) - 1             # order of b(z)
  if (ob < 0) {
    b <- 1
    ob <- 0
  }
  oc <- oa + ob                 # order of c(z)

  c <- fftconv(b, rev(Conj(a)))
  cr <- c * (0:oc)
  num <- stats::fft(postpad(cr, nfft))
  den <- stats::fft(postpad(c, nfft))
  polebins <- which(abs(den) < 2 * .Machine$double.eps)
  if (any(polebins)) {
    warning("setting group delay to 0 at singularity")
    num[polebins] <- 0
    den[polebins] <- 1
  }

  gd <- Re(num / den) - oa

  if (!whole) {
    ns <- nfft / 2     # Matlab convention ... should be nfft/2 + 1
    gd <- gd[1:ns]
    w <- w[1:ns]
  } else {
    ns <- nfft # used in plot below
  }
  res <- list(gd = gd, w = w, ns = ns, HzFlag = HzFlag)
  class(res) <- "grpdelay"
  res
}

#' @rdname grpdelay
#' @export

grpdelay.Arma <- function(filt, ...) # IIR
  grpdelay(filt$b, filt$a, ...)

#' @rdname grpdelay
#' @export

grpdelay.Ma <- function(filt, ...) # FIR
  grpdelay(as.Arma(filt), ...)

#' @rdname grpdelay
#' @export

grpdelay.Sos <- function(filt, ...) # Second-order sections
  grpdelay(as.Arma(filt), ...)

#' @rdname grpdelay
#' @export

grpdelay.Zpg <- function(filt, ...) # Zero-pole-gain ARMA
  grpdelay(as.Arma(filt), ...)
