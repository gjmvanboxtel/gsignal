# fwhm.R
# Copyright (C) 2019 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Octave signal package:
# Author: Petr Mikulik (2009)
# This program is granted to the public domain.
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
# 20200415 Geert van Boxtel          First version for v0.1.0
#------------------------------------------------------------------------------

#' Full width at half maximum
#'
#' Compute peak full-width at half maximum or at another level of peak maximum
#' for a vector or matrix.
#'
#' @param x samples at which \code{y} is measured, specified as a vector. I.e.,
#'   \code{y} is sampled as \code{y[x]}. Default: \code{seq_len(length(y))}.
#' @param y signal to find the width of. If \code{y} is a matrix, widths of all
#'   columns are computed.
#' @param ref reference. Compute the width with reference to:
#' \describe{
#'   \item{\code{"max" | "zero"}}{\code{max(y)}}
#'   \item{\code{"middle" | "min"}}{\code{min(y) + max(y)}}
#'   \item{\code{"absolute"}}{an absolute level of \code{y}}
#' }
#' @param level the level at which to compute the width. Default: 0.5.
#'
#' @return Full width at half maximum, returned as a vector with a length equal
#'   to the number of columns in \code{y}, or 1 in case of a vector.
#'
#' @examples
#' x <- seq(-pi, pi, 0.001)
#' y <- cos(x)
#' w <- fwhm(x, y)
#' m <- x[which.max(y)]
#' f <- m - w/2
#' t <- m + w/2
#' plot(x, y, type="l",
#'      panel.first = {
#'        usr <- par('usr')
#'        rect(f, usr[3], t, usr[4], col = rgb(0, 1, 0, 0.4), border = NA)
#'      })
#' abline(h = max(y) / 2, lty = 2, col = "gray")
#'
#' @author Petr Mikulik.\cr
#'  Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#
#' @export

fwhm <- function(x = seq_len(length(y)), y,
                 ref = c("max", "zero", "middle", "min", "absolute"),
                 level = 0.5) {

  if (!is.vector(x)) {
    stop("x must be a vector")
  }
  if (length(x) != NROW(y)) {
    stop("length of x must match length or number of rows in y")
  }
  if (is.vector(y)) {
    y <- as.matrix(y)
  }
  if (is.array(y) && length(dim(y)) > 2) {
    stop("y must be a vector or a matrix")
  }
  nc <- ncol(y)
  ref <- match.arg(ref)
  if (!is.numeric(level)) {
    stop("level must be numeric")
  }

  w <- rep(0L, nc)
  for (icol in seq_len(nc)) {
    yy <- y[, icol]
    ly <- length(yy)
    if (ref == "absolute") {
      yy <- yy - level
    } else if (ref == "max" || ref == "zero") {
      yy <- yy - level * max(yy)
    } else if (ref == "middle" || ref == "min") {
      yy <- yy - level * (max(yy) + min(yy))
    }
    ind <- which(yy[1:(ly - 1)] * yy[2:ly] <= 0)
    if (length(ind) >= 2 && yy[ind[1]] > 0) {
      ind <- ind[2:length(ind)]
    }
    imax <- which.max(yy)[1]
    li <- length(ind)
    if (li >= 2 && imax >= ind[1] && imax <= ind[li]) {
      ind1 <- ind[1]
      ind2 <- ind1 + 1
      xx1 <- x[ind1] - yy[ind1] * (x[ind2] - x[ind1]) / (yy[ind2] - yy[ind1])
      ind1 <- ind[li]
      ind2 <- ind1 + 1
      xx2 <- x[ind1] - yy[ind1] * (x[ind2] - x[ind1]) / (yy[ind2] - yy[ind1])
      w[icol] <- xx2 - xx1
    }
  }
  w
}
