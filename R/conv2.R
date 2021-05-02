# conv2.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
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
# 20200227  GvB       setup for gsignal v0.1.0
# 20200228  GvB       coerce inputs a and to to matrices instead of checking
#------------------------------------------------------------------------------

#' 2-D convolution
#'
#' Compute the two-dimensional convolution of two matrices.
#'
#' @param a,b Input matrices, coerced to numeric.
#' @param shape Subsection of convolution, partially matched to:
#' \describe{
#'   \item{"full"}{Return the full convolution (default)}
#'   \item{"same"}{Return the central part of the convolution with the same size
#'   as A. The central part of the convolution begins at the indices
#'   \code{floor(c(nrow(b), ncol(b)) / 2 + 1)}}
#'   \item{"valid"}{Return only the parts which do not include zero-padded
#'   edges. The size of the result is \code{max(nrow(a) - nrow(a) + 1, 0)} by
#'   \code{max(ncol(A) - ncol(B) + 1, 0)}}
#' }
#'
#' @return Convolution of input matrices, returned as a matrix.
#'
#' @examples
#' a <- matrix(1:16, 4, 4)
#' b <- matrix(1:9, 3,3)
#' cnv <- conv2(a, b)
#' cnv <- conv2(a, b, "same")
#' cnv <- conv2(a, b, "valid")
#'
#' @seealso \code{\link{conv}}, \code{\link[stats]{convolve}}
#'
#' @author Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#'
#' @export

conv2 <- function(a, b, shape = c("full", "same", "valid")) {

  a <- as.matrix(a)
  b <- as.matrix(b)
  shape <- match.arg(shape)

  if (length(a) < length(b)) {
    x <- a
    a <- b
    b <- x
  }

  y <- switch(shape,
              "full"  = .Call("_gsignal_conv2df",
                              PACKAGE = "gsignal", Re(a), Re(b)),
              "same"  = .Call("_gsignal_conv2ds",
                              PACKAGE = "gsignal", Re(a), Re(b)),
              "valid" = .Call("_gsignal_conv2dv",
                              PACKAGE = "gsignal", Re(a), Re(b))
              )
  y
}
