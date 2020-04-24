# impz.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Original Octave function:
# Copyright (C) 1999 Paul Kienzle <pkienzle@users.sf.net>
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
# See also: http://www.gnu.org/licenses/gpl-2.0.txt
#
# Version history
# 20200423  GvB       setup for gsignal v0.1.0
#---------------------------------------------------------------------------------------------------------------------

#' Impulse response of digital filter
#'
#' Compute the z-plane impulse response of an ARMA model or rational IIR
#' filter. A plot of the impulse and step responses are generated.
#' 
#' @note When results of \code{impz} are printed, \code{plot} will be called to
#'   display a plot of the impulse response against frequency. As with lattice
#'   plots, automatic printing does not work inside loops and function calls, so
#'   explicit calls to print or plot are needed there.
#'
#' @param filt for the default case, the moving-average coefficients of an ARMA
#'   model or filter. Generically, filt specifies an arbitrary model or filter
#'   operation.
#' @param a the autoregressive (recursive) coefficients of an ARMA filter.
#' @param n	number of points at which to evaluate the frequency response. If
#'   \code{n} is a vector with a length greater than 1, then evaluate the
#'   frequency response at these points. For fastest computation, \code{n}
#'   should factor into a small number of small primes. Default: 512.
#' @param fs sampling frequency in Hz. If not specified (default = 2 * pi), the
#'   frequencies are in radians.
#' @param x	object to be printed or plotted.
#' @param ...	 for methods of \code{freqz}, arguments are passed to the default
#'   method. For \code{plot.impz}, additional arguments are passed through to plot.
#'
#' @return For \code{impz}, a list of class \code{'impz'} with items:
#' \describe{
#'   \item{x}{impulse response signal.}
#'   \item{t}{time.}
#' }
#'
#' @examples
#' ## elliptic low-pass filter
#' #elp <- ellip(4, 0.5, 20, 0.4)
#' elp <- Arma(b = c(0.1810488, 0.1048651, 0.3111133, 0.1048651, 0.1810488),
#'             a = c(1.0000000, -1.1463305, 1.5092587, -0.6975320,  0.2698624))
#' impz(elp)
#' 
#' ## or:
#' (xt <- impz(elp))
#' 
#'
#' @author Paul Kienzle \email{pkienzle@@users.sf.net}. Port to R
#'   by Tom Short; adapted by Geert van Boxtel \email{gjmvanboxtel@@gmail.com}
#'
#' @rdname impz
#' @export

impz <- function(filt, ...) UseMethod("impz")

#' @rdname impz
#' @export

print.impz <- plot.impz <- function(x, ...) {

  mini <- min(x$x)
  maxi <- max(x$x)
  op <- graphics::par(mfrow = c(2, 1), mar = c(4, 4, 1.5, 1))
  on.exit(graphics::par(op))
  graphics::plot(x$t, x$x, type = "l", xlab = "", ylab = "Impulse response",
                 ylim=c(min(0, mini), max(1, maxi)), main = "", yaxp = c(0, 1, 1), ...)
  graphics::abline(h = 0, col="red")
  graphics::arrows(0, 0, 0, 1, col = "red", length = 0.1)
  
  step <- cumsum(x$x)
  mini <- min(step)
  maxi <- max(step)
  graphics::plot(x$t, cumsum(x$x), type = "l", xlab = "", ylab = "Step response", 
                 ylim=c(min(0, mini), max(1, maxi)), main = "", yaxp = c(0, 1, 1), ...)
  graphics::segments(0, 0, 0, 1, col="red", lty = 2)
  graphics::segments(0, 1, x$t[length(x$t)], 1, col="red", lty = 2)
  
}

#' @rdname impz
#' @export

impz.Arma <- function(filt, ...) # IIR
  impz(filt$b, filt$a, ...)

#' @rdname impz
#' @export

impz.Ma <- function(filt, ...) # FIR
  impz(filt$b, 1, ...)

#' @rdname impz
#' @export

impz.Sos <- function(filt, ...) # second-order sections
  impz(as.Arma(filt), ...)

#' @rdname impz
#' @export

impz.default <- function(filt, a = 1, n = NULL, fs = 1, ...)  {

  b = filt

  if (length(n) == 0 && length(a) > 1) {
    precision <- 1e-6
    r <- pracma::roots(a)
    maxpole <- max(abs(r))
    if (maxpole > 1 + precision) {        # unstable -- cutoff at 120 dB
      n <- floor(6 / log10(maxpole))
    } else if (maxpole < 1 - precision) { # stable -- cutoff at -120 dB
      n <- floor(-6 / log10(maxpole))
    } else {                              # periodic -- cutoff after 5 cycles
      n <- 30
      
      ## find longest period less than infinity
      ## cutoff after 5 cycles (w=10*pi)
      rperiodic <- r[abs(r) >= 1 - precision & abs(Arg(r)) > 0]
      if (!is.null(rperiodic) && length(rperiodic) > 0) {
        n_periodic <- ceiling( 10 * pi / min(abs(Arg(rperiodic))))
        if (n_periodic > n) {
          n <- n_periodic
        }
      } 
      
      ## find most damped pole
      ## cutoff at -60 dB
      rdamped <- r[abs(r) < 1 - precision]
      if (!is.null(rdamped) && length(rdamped) > 0 ) {
        n_damped <- floor(-3 / log10(max(abs(rdamped))))
      }
      if (n_damped > n) {
        n <- n_damped
      }
    } 
    n <- n + length(b)
  } else if (is.null(n)) {
    n <- length(b)
  } else if (length(n) > 1) {
    t <- n
    n <- length(t)
  }
  if (length(a) == 1) {
    x <- fftfilt(b / a, c(1, numeric(n - 1)))
  } else {
    x <- filter(b, a, c(1, numeric(n - 1)))
  }
  if(!exists("t", mode = "numeric")) {
    t <- (0:(n-1)) / fs
  }
  
  res <- list(x = x, t = t)
  class(res) = "impz"
  res
}
