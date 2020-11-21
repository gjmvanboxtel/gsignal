# downsample.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Original Octave code:
# Copyright (C) 2007 Paul Kienzle <pkienzle@users.sf.net>
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
# 20201129  GvB       setup for gsignal v0.1.0
#---------------------------------------------------------------------------------------------------------------------

#' Decrease sample rate
#' 
#' Downsample a signal by an integer factor.
#' 
#' For most signals you will want to use \code{\link{decimate}} instead since it
#' prefilters the high frequency components of the signal and avoids aliasing
#' effects.
#'   
#' @param x input data, specified as a numeric vector or matrix. In case of a
#'   vector it represents a single signal; in case of a matrix each column is a
#'   signal.
#' @param n downsampling factor, specified as a positive integer.
#' @param phase offset, specified as a positive integer from \code{0} to \code{n
#'   - 1}. Default: 0.
#' 
#' @return downsampled signal, returned as a vector or matrix.
#' 
#' @examples
#' x <- seq_len(10)
#' downsample(x, 3)
#' downsample(x, 3, 2)
#' 
#' x <- matrix(seq_len(12), 4, 3, byrow = TRUE)
#' downsample(x, 3)
#' 
#' @seealso \code{\link{decimate}}, \code{\link{resample}}
#'
#' @author Paul Kienzle, \email{pkienzle@@users.sf.net}; port to R by Geert van
#'   Boxtel \email{G.J.M.vanBoxtel@@gmail.com}.
#
#' @export

downsample <- function (x, n, phase = 0) {

  if (is.vector(x)) {
    ns <- 1
    lx <- length(x)
    x <- matrix(x, ncol = 1)
    vec <- TRUE
  } else if (is.matrix(x)) {
    ns <- ncol(x)
    lx <- nrow(x)
    vec <- FALSE
  } else {
    stop ('x must be a numeric vector or matrix')
  }
  
  if (!(isPosscal(n) && isWhole(n))) {
    stop('n must be a positive integer')
  }
  if (!(isPosscal(phase) && isWhole(phase)) || phase > n - 1) {
    stop('phase must be a positive integer between 0 and n - 1')
  }
  
  y <- x[seq(phase + 1, lx, n), ]
  
  if (vec) {
    y <- as.vector(y)
  }
  y  
}
