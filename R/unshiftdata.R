# unshiftdata.R
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
# 20200507 GvB                       Bugfix
#------------------------------------------------------------------------------

#' Inverse of shiftdata
#'
#' Reverse what has been done by \code{shiftdata()}.
#'
#' \code{unshiftdata} restores the orientation of the data that was shifted with
#' shiftdata. The permutation vector is given by \code{perm}, and \code{nshifts}
#' is the number of shifts that was returned from \code{shiftdata()}.
#'
#' \code{unshiftdata} is meant to be used in tandem with \code{shiftdata}. These
#' functions are useful for creating functions that work along a certain
#' dimension, like filter, goertzel, sgolayfilt, and sosfilt. These functions
#' are useful for creating functions that work along a certain dimension, like
#' \code{\link{filter}}, \code{\link{sgolayfilt}}, and \code{\link{sosfilt}}.
#'
#' @param sd A list of objects named \code{x}, \code{perm}, and \code{nshifts},
#'   as returned by \code{shiftdata()}
#'
#' @return Array with the same values and dimensions as passed to a previous
#'   call to \code{shiftdata}.
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
#' ## Rearrange Array to Operate on First Nonsingleton Dimension
#' x <- 1:5
#' sd <- shiftdata(x)
#' y <- unshiftdata(sd)
#'
#' @author Georgios Ouzounis, \email{ouzounis_georgios@@hotmail.com}.\cr
#'   Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#'
#' @seealso \code{\link{shiftdata}}
#'
#' @export

unshiftdata <- function(sd) {

  nm <- names(sd)
  if (!is.list(sd) || !("x" %in% nm && "perm" %in% nm && "nshifts" %in% nm))
    stop("sd must be a list with elements x, perm and nshifts")

  if (length(sd$perm) > 0 && !anyNA(sd$perm)) {
    if (!isWhole(sd$perm))
      stop(paste0(deparse(substitute(sd)),
                  "$perm must be a vector of integers"))
    dimx <- sd$perm[1]
  } else if (length(sd$nshifts) > 0) {
    if (!isWhole(sd$nshifts))
      stop(paste0(deparse(substitute(sd)), "$nshifts must be an integer"))
    dimx <- sd$nshifts + 1
  } else {
    stop(paste0("Either perm or nshifts must not be empty"))
  }

  perm <- dimx
  if (dimx - 1 >= 1) {
    perm <- c(dimx, 1:(dimx - 1))
  }
  d1 <- dimx + 1
  d2 <- (length(dim(sd$x)))
  if (d1 <= d2) perm <- c(perm, d1:d2)

  iaperm <- function(x, p) {
    p[p] <- seq_along(dim(x))
    aperm(x, p)
  }
  out <- iaperm(sd$x, perm)
  out
}
