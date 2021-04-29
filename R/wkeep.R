# wkeep.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Original Octave code:
# Copyright (C) 2008 Sylvain Pelissier <sylvain.pelissier@gmail.com>
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
# 20201128  GvB       setup for gsignal v0.1.0
#------------------------------------------------------------------------------

#' Keep part of vector or matrix
#'
#' Extract elements from a vector or matrix.
#'
#' @param x input data, specified as a numeric vector or matrix.
#' @param l either a positive integer value, specifying the length to extract
#'   from the input *vector* \code{x}, or a vector of length 2, indicating the
#'   submatrix to extract from the *matrix* \code{x}. See the examples.
#' @param opt One of:
#' \describe{
#'   \item{character string}{matched against \code{c("centered", "left",
#'   "right")}, indicating the location of the *vector* \code{x} to extract}
#'   \item{positive integer}{starting index of the input *vector* \code{x}}
#'   \item{two-element vector}{starting row and columns from the *matrix*
#'   \code{x}}
#' }
#' See the examples. Default: "centered".
#'
#' @return extracted vector or matrix
#'
#' @examples
#' ## create a vector
#' x <- 1:10
#' ## Extract a vector of length 6 from the central part of x.
#' y <- wkeep(x, 6, 'c')
#'
#' ## Extract two vectors of length 6, one from the left part of x, and the
#' ## other from the right part of x.
#' y <- wkeep(x, 6, 'l')
#' y <- wkeep(x, 6, 'r')
#'
#' ## Create a 5-by-5 matrix.
#' x <- matrix(round(runif(25, 0, 25)), 5, 5)
#'
#' ## Extract a 3-by-2 matrix from the center of x
#' y <- wkeep(x, c(3, 2))
#'
#' ## Extract from x the 2-by-4 submatrix starting at x[3, 1].
#' y <- wkeep(x, c(2, 4), c(3, 1))
#'
#' @author Sylvain Pelissier, \email{sylvain.pelissier@@gmail.com}.\cr
#' Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#
#' @export

wkeep <- function(x, l, opt = "centered") {

  if (is.vector(x) && !is.character(x)) {
    lx <- length(x)
    if (!isPosscal(l) || l > length(x)) {
      stop("l must be a positive integer <= length(x)")
    }
    if (is.character(opt)) {
      opt <- match.arg(opt, c("centered", "left", "right"))
      if (opt == "centered") {
        s <- (lx - l) / 2
        y <- x[(1 + floor(s)):(lx - ceiling(s))]
      } else if (opt == "left") {
        y <- x[1:l]
      } else if (opt == "right") {
        y <- x[(lx - l + 1):lx]
      } else {
        stop('opt must be a character string ("centered", "left", or "right")')
      }
    } else if (isPosscal(opt)) {
      if (opt == 0 || opt + l - 1 > lx) {
        stop(paste("opt must be an integer value between 1 and", lx - l))
      } else {
        y <- x[opt:(opt + l - 1)]
      }
    } else {
      stop("opt must be a character string or a positive scalar value")
    }
  } else if (is.matrix(x)) {
    nr <- nrow(x)
    nc <- ncol(x)
    lx <- max(nr, nr)
    if (!is.vector(l) || length(l) != 2) {
      stop("When x is a matrix l must be a vector of length 2")
    } else {
      s1 <- (lx - l[1]) / 2
      s2 <- (lx - l[2]) / 2
    }
    if (!is.numeric(opt)) {
      y <- x[(1 + floor(s1)):(nc - ceiling(s1)),
             (1 + floor(s2)):(nc - ceiling(s2))]
    } else {
      if (length(opt) == 2) {
        firstr <- opt[1]
        firstc <- opt[2]
      } else {
        stop("When x is a matrix opt must be a vector of length 2")
      }
      y <- x[firstr:(firstr + l[1] - 1), firstc:(firstc + l[2] - 1)]
    }
  } else {
    stop("x must be a vector or a matrix")
  }
  y
}
