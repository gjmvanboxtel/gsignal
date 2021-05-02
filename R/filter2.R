# filter2.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Original Octave code:
# Copyright (C) 2001-2019 Paul Kienzle
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
# 20201130  GvB       setup for gsignal v0.1.0
#------------------------------------------------------------------------------

#' 2-D digital filter
#'
#' Apply a 2-D digital filter to the data in \code{x}.
#'
#' The \code{filter2} function filters data by taking the 2-D convolution of the
#' input \code{x} and the coefficient matrix \code{h} rotated 180 degrees. More
#' specifically, \code{filter2(h, x, shape)} is equivalent to \code{conv2(x,
#' rot90(h, 2), shape)}.
#'
#' @param h transfer function, specified as a matrix.
#' @param x numeric matrix containing the input signal to be filtered.
#' @param shape Subsection of convolution, partially matched to:
#' \describe{
#'   \item{"same"}{Return the central part of the filtered data; same size as
#'   \code{x} (Default)}
#'   \item{"full"}{Return the full 2-D filtered data, with zero-padding on all
#'   sides before filtering}
#'   \item{"valid"}{Return only the parts which do not include zero-padded
#'   edges.}
#' }
#'
#' @return The filtered signal, returned as a matrix
#'
#' @examples
#' op <- par(mfcol = c(1, 2))
#' x <- seq(-10, 10, length.out = 30)
#' y <- x
#' f <- function(x, y) { r <- sqrt(x^2+y^2); 10 * sin(r)/r }
#' z <- outer(x, y, f)
#' z[is.na(z)] <- 1
#' persp(x, y, z, theta = 30, phi = 30, expand = 0.5, col = "lightblue")
#' title( main = "Original")
#'
#' h <- matrix(c(1, -2, 1, -2, 3, -2, 1, -2, 1), 3, 3)
#' zf <-filter2(h, z, 'same')
#' persp(x, y, zf, theta = 30, phi = 30, expand = 0.5, col = "lightgreen")
#' title( main = "Filtered")
#' par(op)
#'
#' @seealso \code{\link{conv2}}
#'
#' @author Paul Kienzle.
#' Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#'
#' @export

filter2 <- function(h, x, shape = c("same", "full", "valid")) {

  if (!is.numeric(h) || !is.matrix(h)) {
    stop("h must be a numeric matrix")
  }
  if (!is.numeric(x) || !is.matrix(x)) {
    stop("x must be a numeric matrix")
  }
  shape <- match.arg(shape)

  nr <- nrow(h)
  nc <- ncol(h)

  y <- conv2(x, h[seq(nr, 1, -1), seq(nc, 1, -1)], shape)
  y
}
