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
# 20200425  GvB       setup for gsignal v0.1.0
#---------------------------------------------------------------------------------------------------------------------

#' Zero-pole plot
#'
#' Plot the poles and zeros of a filter or model on the complex Z-plane
#' 
#' Poles are marked with an ‘x’, and zeros are marked with an ‘o’.
#' 
#' @note When results of \code{zplane} are printed, \code{plot} will be called.
#'   As with lattice plots, automatic printing does not work inside loops and
#'   function calls, so explicit calls to print or plot are needed there.
#'
#' @param filt for the default case, the moving-average coefficients of an ARMA
#'   model or filter. Generically, filt specifies an arbitrary model or filter
#'   operation.
#' @param a the autoregressive (recursive) coefficients of an ARMA filter.
#' @param ...	additional arguments are passed through to plot.
#'
#' @return No value is returned.
#'
#' @examples
#' ## elliptic low-pass filter
#' #elp <- ellip(4, 0.5, 20, 0.4)
#' elp <- Arma(b = c(0.1810488, 0.1048651, 0.3111133, 0.1048651, 0.1810488),
#'             a = c(1.0000000, -1.1463305, 1.5092587, -0.6975320,  0.2698624))
#' zplane(elp)
#' 
#' @references \url{http://en.wikipedia.org/wiki/Pole-zero_plot}
#' 
#' @seealso \code{\link{freqz}}
#' 
#' @author Paul Kienzle \email{pkienzle@@users.sf.net}, Stefan van der Walt
#'   \email{stefan@@sun.ac.za}, Mike Miller. Port to R by Tom Short; adapted by
#'   Geert van Boxtel \email{gjmvanboxtel@@gmail.com}
#'
#' @rdname zplane
#' @export

zplane <- function(filt, ...) UseMethod("zplane")

#' @rdname zplane
#' @export

zplane.Arma <- function(filt, ...) # IIR
  zplane(filt$b, filt$a, ...)

#' @rdname zplane
#' @export

zplane.Ma <- function(filt, ...) # FIR
  zplane(filt, 1, ...)

#' @rdname zplane
#' @export

zplane.Sos <- function(filt, ...)
  zplane(as.Zpg(filt, ...))

#' @rdname zplane
#' @export

zplane.Zpg <- function(filt, ...) {
  x <- filt
  r <- exp(2i * pi * (0:100) / 100)
  xlim <- range(c(-1.1, 1.1, Re(x$p), Re(x$z)))
  ylim <- range(c(-1.1, 1.1, Im(x$p), Im(x$z)))
  graphics::plot(Re(r), Im(r), col = "red", xlab = "", ylab = "", xlim = xlim, ylim = ylim, type = "l", asp = 1, ...)
  graphics::points(Re(x$p), Im(x$p), pch = 4)
  graphics::points(Re(x$z), Im(x$z), pch = 1)
}

#' @rdname zplane
#' @export

zplane.default <- function(filt, a, ...) {
  zplane(Zpg(pracma::roots(filt), pracma::roots(a), 1), ...)
}
