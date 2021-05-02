# bitrevorder.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Original Octave code:
# Copyright (C) 2013-2019 Mike Miller
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
# 20200821  GvB       setup for gsignal v0.1.0
#------------------------------------------------------------------------------

#' Permute input to bit-reversed order
#'
#' Reorder the elements of the input vector in bit-reversed order.
#'
#' This function is equivalent to calling \code{digitrevorder(x, 2)}, and is
#' useful for prearranging filter coefficients so that bit-reversed ordering
#' does not have to be performed as part of an fft or ifft computation.
#'
#' @param x input data, specified as a vector. The length of \code{x} must be an
#'   integer power of 2.
#' @param index.return logical indicating if the ordering index vector should be
#'   returned as well. Default: \code{FALSE}.
#'
#' @return The bit-reversed input vector. If \code{index.return = TRUE}, then
#'   a list containing the bit-reversed input vector (\code{y}), and the
#'   digit-reversed indices (\code{i}).
#'
#' @examples
#' x <- 0:15
#' v <- bitrevorder(x)
#' dec2bin <- function(x, l)
#'   substr(paste(as.integer(rev(intToBits(x))), collapse = ""),
#'   32 - l + 1, 32)
#' x_bin <- sapply(x, dec2bin, 4)
#' v_bin <- sapply(v, dec2bin, 4)
#' data.frame(x, x_bin, v, v_bin)
#'
#' @author Mike Miller.\cr
#'  Port to to by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#'
#' @seealso \code{\link{digitrevorder}}, \code{\link{fft}}, \code{\link{ifft}}
#'
#' @export

bitrevorder <- function(x, index.return = FALSE) {

  if (!is.vector(x)) {
    stop("x must be a vector")
  } else if (trunc(log2(length(x))) != log2(length(x))) {
    stop("x must have length equal to an integer power of 2")
  }
  if (!is.logical(index.return)) {
    stop("index.return must be TRUE or FALSE")
  }

  digitrevorder(x, 2, index.return)
}
