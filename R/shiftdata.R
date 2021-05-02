# shiftdata.R
# Copyright (C) 2019 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Matlab/Octave signal package:
# Copyright (C) 2014 Georgios Ouzounis <ouzounis_georgios@hotmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; see the file COPYING. If not, see
# <https://www.gnu.org/licenses/>.
#
# 20191202 Geert van Boxtel          First version for v0.1.0
#------------------------------------------------------------------------------

#' Shift data to operate on specified dimension
#'
#' Shift data in to permute the dimension \code{dimx} to the first column.
#'
#' \code{shiftdata(x, dimx)} shifts data \code{x} to permute dimension
#' \code{dimx} to the first column using the same permutation as the built-in
#' \code{filter} function. The vector \code{perm} in the output list returns the
#' permutation vector that is used.
#'
#' If \code{dimx} is missing or empty, then the first nonsingleton dimension is
#' shifted to the first column, and the number of shifts is returned in
#' \code{nshifts}.
#'
#' \code{shiftdata} is meant to be used in tandem with \code{unshiftdata}, which
#' shifts the data back to its original shape. These functions are useful for
#' creating functions that work along a certain dimension, like
#' \code{\link{filter}}, \code{\link{sgolayfilt}}, and \code{\link{sosfilt}}.
#'
#' @param x The data to be shifted. Can be of any type.
#' @param dimx Dimension of \code{x} to be shifted to the first column. Named
#'   "dimx" instead of "dim" to avoid confusion with R's dim() function.
#'   Default: NULL (shift the first nonsingleton dimension)
#'
#' @return A list containing 3 variables; \code{x}, the shifted data,
#'   \code{perm}, the permutation vector, and \code{nshifts}, the number of
#'   shifts
#'
#' @examples
#'
#' ## create a 3x3 magic square
#' x <- pracma::magic(3)
#' ## Shift the matrix x to work along the second dimension.
#' ## The permutation vector, perm, and the number of shifts, nshifts,
#' ## are returned along with the shifted matrix.
#' sd <- shiftdata(x, 2)
#'
#' ## Shift the matrix back to its original shape.
#' y <- unshiftdata(sd)
#'
#' ## Rearrange Array to Operate on First nonsingleton Dimension
#' x <- 1:5
#' sd <- shiftdata(x)
#' y <- unshiftdata(sd)
#'
#' @author Georgios Ouzounis, \email{ouzounis_georgios@@hotmail.com}.\cr
#' Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#'
#' @seealso \code{\link{unshiftdata}}
#'
#' @export

shiftdata <- function(x, dimx) {

  x <- as.array(x)   # needed for aperm

  if (missing(dimx) || is.null(dimx)) {
    dimx <- which((dim(x) - 1) > 0)[1]
    shift <- TRUE
  } else {
    shift <- FALSE
  }

  if (!isScalar(dimx) || !isWhole(dimx))
    stop("dimx must be an integer")
  if (dimx > length(dim(x)))
    stop(paste("dimx should be between 1 and", length(dim(x))))

  perm <- dimx
  if (dimx - 1 >= 1) {
    perm <- c(dimx, 1:(dimx - 1))
  }
  d1 <- dimx + 1
  d2 <- length(dim(x))
  if (d1 <= d2) perm <- c(perm, d1:d2)
  out <- aperm(x, perm)

  if (shift) {
    perm <- NA
    nshifts <- dimx - 1
  } else {
    nshifts <- NA
  }

  list(x = out, perm = perm, nshifts = nshifts)
}
