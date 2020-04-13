# filter.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
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
# 20200208  GvB       setup for gsignal v0.1.0
# 20200413  GvB       added S3 method for Sos
#---------------------------------------------------------------------------------------------------------------------

#' Filter a signal
#' 
#' Apply a 1-D digital filter to the data in \code{x}, compatible with
#' Matlab/Octave.
#' 
#' \code{filter(b, a, x)} returns the solution to the following linear,
#' time-invariant difference equation:
#' \deqn{\sum_{k=0}^{N} a(k+1) y(n-k) + \sum_{k=0}^{M} b(k+1) x(n-k) = 0, 1 \le n \le length(x)}{\sum k=0:N a(k+1) y(n-k) + \sum k=0:M b(k+1) x(n-k) = 0, 1 \le n \le length(x)}
#' where \code{N = length(a) - 1} and \code{M = length(b) - 1}.
#' 
#' An equivalent form of this equation is:
#' \deqn{y(n) = \sum_{k=1}^{N} c(k+1) y(n-k) + \sum_{k=0}^{M} d(k+1) x(n-k), 1 \le n \le length(x)}{y(n) = \sum k=1:N c(k+1) y(n-k) + \sum k=0:M d(k+1) x(n-k), 1 \le n \le length(x)}
#' where \code{c = a/a(1) and d = b/a(1)}.
#' 
#' Specifying \code{init.x}, \code{init.y}, or \code{init} (alias of
#' \code{init.y}) allows filtering very large time series in pieces in case of
#' limited computer memory. This works by passing the last \code{M} data points
#' of the truncated input series \code{x} for the convolution part to
#' \code{init.x}, or the last \code{N} data points of the truncated output
#' series \code{y} to \code{init.y}. See the examples.
#' 
#' The default filter calls \code{\link[stats]{filter}}, which returns an object
#' of class \code{\link{ts}}. However, the output of \code{filter} is converted
#' to a vector before the result is returned.
#' 
#' @param filt For the default case, the moving-average coefficients of an ARMA
#'   filter (normally called ‘b’). Generically, \code{filt} specifies an arbitrary
#'   filter operation.
#' @param a the autoregressive (recursive) coefficients of an ARMA filter.
#' @param x the input signal to be filtered.
#' @param init.x initial data for the convolution part of the filter (FIR).
#' @param init.y,init initial data for the recursive part of the filter (IIR).
#' @param ... additional arguments (ignored).
#' 
#' @return The filtered signal, normally of the same length as the input signal
#'   \code{x}, returned as a vector
#' 
#' @examples
#' #bf <- butter(3, 0.1)                                # 10 Hz low-pass filter
#' bf <- Arma(c(0.002898195, 0.008694584, 0.008694584, 0.002898195),
#'  c(1.0000000, -2.3740947,  1.9293557, -0.5320754))   # change later
#' t <- seq(0, 1, len = 100)                            # 1 second sample
#' x <- sin(2* pi * t * 2.3) + 0.25 * rnorm(length(t))  # 2.3 Hz sinusoid+noise
#' z <- filter(bf, x)                                   # apply filter
#' plot(t, x, type = "l")
#' lines(t, z, col = "red")
#' 
#' ## example of filtering data in pieces (FIR)
#' b <- c(1, 1); a <- 1
#' x <- 1:100
#' y <- filter(b, a, x)
#' x1 <- x[1:50]
#' x2 <- x[51:100]
#' y1 <- filter(b, a, x1)
#' y2 <- filter(b, a, x2, init.x = x1[(length(x1) - (length(b) - 1) + 1):length(x1)])
#' all.equal(as.numeric(y), c(y1,y2))
#' 
#' ## example of filtering data in pieces (IIR)
#' b <- 1; a <- c(1,1)
#' x <- 1:100
#' y <- filter(b, a, x)
#' x1 <- x[1:50]
#' x2 <- x[51:100]
#' y1 <- filter(b, a, x1)
#' y2 <- filter(b, a, x2, init = y1[(length(y1) - (length(a) - 1) + 1):length(y1)])
#' all.equal(as.numeric(y), c(y1,y2))
#' 
#' @seealso \code{\link[stats]{filter}} in the \code{stats} package
#' 
#' @author Tom Short, EPRI Solutions, Inc., (\email{tshort@@eprisolutions.com}), slightly
#' adapted by Geert van Boxtel \email{G.J.M.vanBoxtel@@gmail.com}.
#' 
#' @rdname filter
#' @export

filter <- function(filt, ...) UseMethod("filter")

#' @rdname filter
#' @method filter default
#' @export

filter.default <- function(filt, a, x, init, init.x, init.y, ...) {
  
  if(missing(init.x)) {
    init.x <- c(rep(0, length(filt) - 1))
  }
  
  if(length(init.x) != length(filt) - 1) {
    stop("length of init.x should match filter length-1 = ", length(filt)-1)
  }
  
  if(missing(init) && !missing(init.y)) {
    init <- rev(init.y)
  }
  
  if(all(is.na(x))) {
    return(x)
  }
  
  if (length(filt)) {
    x1 <- stats::filter(c(init.x, x), filt / a[1], sides = 1)
    if(all(is.na(x1))) {
      return(x)
    }
    x <- stats::na.omit(x1)
  }
  
  if (length(a) >= 2) {
    x <- stats::filter(x, -a[-1] / a[1], method = "recursive", init = init)
  }
  
  as.numeric(x)
}

#' @rdname filter
#' @method filter Arma
#' @export
filter.Arma <- function(filt, x, ...) # IIR
  filter(filt$b, filt$a, x, ...)

#' @rdname filter
#' @method filter Ma
#' @export
filter.Ma <- function(filt, x, ...) # FIR
  filter(unclass(filt), 1, x, ...)

#' @rdname filter
#' @method filter Sos
#' @export
filter.Sos <- function(filt, x, ...) { # Second-order sections
  if (filt$g != 1) {
    filt$sos[1, 1:3] <- filt$sos[1, 1:3] * filt$g
  }
  sosfilt(filt$sos, x, ...)
}

#' @rdname filter
#' @method filter Zpg
#' @export
filter.Zpg <- function(filt, x, ...) # zero-pole-gain form
  filter(as.Arma(filt), x)
