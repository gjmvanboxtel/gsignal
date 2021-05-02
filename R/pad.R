# pad.R
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
# 20200121  GvB       setup for gsignal v0.1.0
#------------------------------------------------------------------------------

#' Pad data
#'
#' Pre- or postpad the data object \code{x} with the value \code{c} until it is
#' of length \code{l}.
#'
#' @param x Vector or matrix to be padded
#' @param l Length of output data along the padding dimension. If \code{length
#'   (x) > l}, elements from the beginning (\code{dimension = "pre"}) or the end
#'   (\code{direction = "post"}) of \code{x} are removed until a vector of
#'   length \code{l} is obtained. If \code{direction = "both"}, values are
#'   removed from both ends, and in case of an uneven length the smallest number
#'   of elements is removed from the beginning of vector.
#' @param c Value to be used for the padding (scalar). Must be of the same type
#'   as the elements in \code{x}. Default: 0
#' @param MARGIN A vector giving the subscripts which the function will be
#'   applied over. E.g., for a matrix 1 indicates rows, 2 indicates columns,
#'   c(1, 2) indicates rows and columns. Where \code{x} has named dimnames, it
#'   can be a character vector selecting dimension names.  If \code{MARGIN} is
#'   larger than the dimensions of \code{x}, the result will have \code{MARGIN}
#'   dimensions. Default: 2 (columns).
#' @param direction Where to pad the array along each dimension. One of the
#'   following:
#' \describe{
#'   \item{"pre"}{Before the first element}
#'   \item{"post"}{After the last element}
#'   \item{"both"}{(default) Before the first and after the last element}
#' }
#'
#' @return Padded data, returned as a vector or matrix.
#'
#' @examples
#' v <- 1:24
#' res <- postpad(v, 30)
#' res <- postpad(v, 20)
#' res <- prepad(v, 30)
#' res <- prepad(v, 20)
#'
#' m <- matrix(1:24, 4, 6)
#' res <- postpad(m, 8, 100)
#' res <- postpad(m, 8, 100, MARGIN = 1)
#' res <- prepad(m, 8, 100)
#' res <- prepad(m, 8, 100, MARGIN = 1)
#'
#' res <- postpad(m, 2)
#' res <- postpad(m, 2, MARGIN = 1)
#' res <- prepad(m, 2)
#' res <- prepad(m, 2, MARGIN = 1)
#'
#' @author Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#' @export

pad <- function(x, l, c = 0, MARGIN = 2,
                direction = c("both", "pre", "post")) {

  vec <- FALSE
  if (is.vector(x)) vec <- TRUE

  if (vec) {
    x <- as.matrix(x)
    MARGIN <- 2
  }
  xdim <- dim(x)
  ld <- length(xdim)

  if (!isPosscal(l) && !isWhole(l)) {
    stop("l must be a positive integer")
  }

  if (!isScalar(c)) {
    stop("c must be a scalar")
  }

  # handle character dimnames (adapted from apply source code)
  if (is.character(MARGIN)) {
    if (is.null(dnn <- names(dimnames(x))))
      stop("'x' must have named dimnames")
    MARGIN <- match(MARGIN, dnn)
    if (anyNA(MARGIN))
      stop("not all elements of 'MARGIN' are names of dimensions")
  } else if (!isPosscal(MARGIN) || MARGIN > ld || MARGIN < 1 || MARGIN > 2) {
    stop("'MARGIN' must be a positive integer and a valid margin (1 or 2)")
  }

  pd <- function(x, l, c, dir) {
    lx <- length(x)
    if (l == lx) {
      y <- x
    } else {
      if (dir == "pre") {
        if (l > lx) {
          y <- c(rep(c, l - lx), x)
        } else {
          y <- x[(lx - l + 1):lx]
        }
      } else if (dir == "post") {
        if (l > lx) {
          y <- c(x, rep(c, l - lx))
        } else {
          y <- x[1:l]
        }
      } else if (dir == "both") {
        r1 <- floor(abs(lx - l) / 2)
        r2 <- ceiling(abs(lx - l) / 2)
        if (l > lx) {
          y <- c(rep(c, r1), x, rep(c, r2))
        } else {
          y <- x[(r1 + 1):(lx - r2)]
        }
      }
    }
    y
  }

  direction <- match.arg(direction)
  y <- unlist(apply(X = x, MARGIN = MARGIN, FUN = pd,
                    l = l, c = c, dir = direction))

  if (MARGIN == 1) y <- t(y)
  if (vec) y <- as.vector(y)
  y
}

#' @rdname pad
#' @export
prepad <- function(x, l, c = 0, MARGIN = 2) {

  pad(x, l, c, MARGIN, direction = "pre")

}
#' @rdname pad
#' @export
postpad <- function(x, l, c = 0, MARGIN = 2) {

  pad(x, l, c, MARGIN, direction = "post")

}
