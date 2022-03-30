# zplane.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Original Octave function:
# Copyright (C) 1999, 2001 Paul Kienzle <pkienzle@users.sf.net>
# Copyright (C) 2004 Stefan van der Walt <stefan@sun.ac.za>
# Copyright (C) 2019 Mike Miller
# R signal version:
# Copyright (C) 2006 EPRI Solutions, Inc.
# by Tom Short, tshort@eprisolutions.com
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
# 20200425  GvB       setup for gsignal v0.1.0
# 20201214  GvB       changes to S3 setup: do them all via as.Zpg
#------------------------------------------------------------------------------

#' Zero-pole plot
#'
#' Plot the poles and zeros of a filter or model on the complex Z-plane
#'
#' Poles are marked with an \code{x}, and zeros are marked with an \code{o}.
#'
#' @note When results of \code{zplane} are printed, \code{plot} will be called.
#'   As with lattice plots, automatic printing does not work inside loops and
#'   function calls, so explicit calls to print or plot are needed there.
#'
#' @param filt for the default case, the moving-average coefficients of an ARMA
#'   model or filter. Generically, \code{filt} specifies an arbitrary model or
#'   filter operation.
#' @param a the autoregressive (recursive) coefficients of an ARMA filter.
#' @param ...	additional arguments are passed through to plot.
#'
#' @return No value is returned.
#'
#' @examples
#' ## elliptic low-pass filter
#' elp <- ellip(4, 0.5, 20, 0.4)
#' zplane(elp)
#'
#' @references \url{https://en.wikipedia.org/wiki/Pole-zero_plot}
#'
#' @seealso \code{\link{freqz}}
#'
#' @author Paul Kienzle, \email{pkienzle@@users.sf.net},\cr
#'  Stefan van der Walt \email{stefan@@sun.ac.za},\cr
#'  Mike Miller.\cr
#'   Conversion to R by Tom Short,\cr
#'    adapted by Geert van Boxtel, \email{gjmvanboxtel@@gmail.com}
#'
#' @rdname zplane
#' @export

zplane <- function(filt, ...) UseMethod("zplane")

#' @rdname zplane
#' @export

zplane.Arma <- function(filt, ...) # IIR
  zplane(as.Zpg(filt), ...)

#' @rdname zplane
#' @export

zplane.Ma <- function(filt, ...) # FIR
  zplane(as.Zpg(filt), ...)

#' @rdname zplane
#' @export

zplane.Sos <- function(filt, ...)
  zplane(as.Zpg(filt), ...)

#' @rdname zplane
#' @export

zplane.Zpg <- function(filt, ...) {
  x <- filt
  r <- exp(2i * pi * (0:100) / 100)
  xlim <- range(c(-1.1, 1.1, Re(x$p), Re(x$z)))
  ylim <- range(c(-1.1, 1.1, Im(x$p), Im(x$z)))
  graphics::plot(Re(r), Im(r), col = "red", xlab = "", ylab = "",
                 xlim = xlim, ylim = ylim, type = "l", asp = 1, ...)
  graphics::points(Re(x$p), Im(x$p), pch = 4)
  graphics::points(Re(x$z), Im(x$z), pch = 1)
}

#' @rdname zplane
#' @export

zplane.default <- function(filt, a, ...) {
  zplane(Zpg(pracma::roots(as.numeric(filt)),
             pracma::roots(as.numeric(a)), 1), ...)
}
