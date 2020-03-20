# medfilt1.R
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
# 20200320    GvB       setup for gsignal v0.1.0
#---------------------------------------------------------------------------------------------------------------------

#' 1-D median filtering
#' 
#' Apply a running median of odd span to the input \code{x}
#' 
#' This function computes a running median over the input \code{x}, using the
#' \code{\link[stats]{runmed}} function. Because of that, it works a little
#' differently than the Matlab/Octave version (i.e., it does not produce exactly
#' the same values). 
#' 
#' \describe{
#' \item{missing values}{The Mablab/Octave function has a \code{'nanflag'}
#'   option that allows to include or remove missing values. If inclusion is
#'   specifies, then the function returns a signal so that the median of any
#'   segment containing NAs is also Na. Because the \code{runmed} function does
#'   not include an \code{na.omit} option, implementing this functionality would
#'   lead to a considerable speed loss. Instead, an \code{na.omit} parameter was
#'   implemented that allows either omitting NAs or interpolating them with a
#'   spline function.}
#' \item{endpoint filtering}{Instead of the \code{'zeropad'} and
#'   \code{'truncate'} options to the \code{'padding'} argument in the
#'   Matlab/Octave function, the present version uses the standard
#'   \code{'endrule'} parameter of the \code{runmed} function, with options
#'   \code{'keep'}, \code{'constant'}, or \code{'median'}.}
#'}
#' 
#' @param x Input signal, specified as a numeric vector, matrix or array.
#' @param n positive integer width of the median window; must be odd. Default: 3
#' @param dim Dimension to filter along in case \code{x} is a matrix or array.
#'   Default: 2 (columns)
#' @param na.omit logical indicating whether to omit missing values,
#'   or interpolate then using a cubic spline function
#'   (\code{\link[stats]{splinefun}}). Default: FALSE
#' @param ... other arguments passed to \code{runmed}
#' 
#' @return Filtered signal, returned as a numeric vector, matrix, or array, of
#'   the same size as \code{x}.
#' 
#' @examples
#' ## noise suppression
#' fs <- 100
#' t <- seq(0, 1, 1/fs)
#' x <- sin(2 * pi * t * 3) + 0.25 * sin(2 * pi * t * 40)
#' plot(t, x, type = "l", xlab = "", ylab = "")
#' y <- medfilt1(x, 11)
#' lines (t, y, col = "red")
#' legend("topright", c("Original", "Filtered"), lty = 1, col = 1:2)
#' 
#' @seealso \code{\link[stats]{runmed}}, \code{\link[stats]{splinefun}}
#' 
#' @author Geert van Boxtel \email{G.J.M.vanBoxtel@@gmail.com}.
#' 
#' @export

medfilt1 <- function(x, n = 3, dim = 2, na.omit = FALSE, ...) {
  
  mf <- function (x, n, na.omit, ...) {
    if (n%%2 != 1 || n > length(x)) {
      stop ('n must be odd and smaller than the length of x')
    }
    if (any(is.na(x))) {
      if (na.omit) {
        x <- na.omit(x)
      } else {
        spl <- stats::splinefun(seq_along(x), x)
        x <- spl(seq_along(x))
      }
    }
    y <- stats::runmed(x, n, ...)
    as.vector(y)
  }
  
  if (is.vector(x)) {
    y <- mf(x, n, na.omit, ...)
  } else {
    d <- dim(x)
    if (dim > length(d)) {
      stop("invalid 'dim' argument")
    }
    y <- apply(x, dim, mf, n, na.omit, ...)
  }
  y
}
