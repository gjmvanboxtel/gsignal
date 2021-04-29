# schtrig.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Original Octave code:
# Copyright (C) 2012 Juan Pablo Carbajal <carbajal@ifi.uzh.ch>
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
# 20201127  GvB       setup for gsignal v0.1.0
#------------------------------------------------------------------------------

#' Schmitt Trigger
#'
#' Multisignal Schmitt trigger with levels.
#'
#' The trigger works compares each column in \code{x} to the levels in
#' \code{lvl}, when the value is higher than \code{max(lvl)}, the output
#' \code{v} is high (i.e. 1); when the value is below \code{min(lvl)} the output
#' is low (i.e. 0); and when the value is between the two levels the output
#' retains its value.
#'
#' @param x input data, specified as a numeric vector or matrix. In case of a
#'   vector it represents a single signal; in case of a matrix each column is a
#'   signal.
#' @param lvl threshold levels against which \code{x} is compared, specified as
#'   a vector. If this is a scalar, the thresholds are symmetric around 0, i.e.
#'   \code{c(-lvl, lvl)}.
#' @param st trigger state, specified as a vector of length \code{ncol(x}. The
#'   trigger state is returned in the output list and may be passed again to a
#'   subsequent call to \code{schtrig}. Default: NULL.
#'
#' @return a \code{\link{list}} containing the following variables:
#' \describe{
#'   \item{v}{vector or matrix of 0's and 1's, according to whether \code{x} is
#'   above or below \code{lvl}, or the value of \code{x} if indeterminate}
#'   \item{rng}{ranges in which the output is high, so the indexes
#'   \code{rng[1,i]:rng[2,i]} point to the i-th segment of 1s in \code{v}. See
#'   \code{\link{clustersegment}} for a detailed explanation.}
#'   \item{st}{trigger state, returned as a vector with a length of the number
#'   of columns in \code{x}.}
#' }
#'
#' @examples
#' t <- seq(0, 1, length.out = 100)
#' x <- sin(2 * pi * 2 * t) + sin(2 * pi * 5 * t) %*% matrix(c(0.8, 0.3), 1, 2)
#' lvl <- c(0.8, 0.25)
#' trig  <- schtrig (x, lvl)
#'
#' op <- par(mfrow = c(2, 1))
#' plot(t, x[, 1], type = "l", xlab = "", ylab = "")
#' abline(h = lvl, col = "blue")
#' lines(t, trig$v[, 1], col = "red", lwd = 2)
#' plot(t, x[, 2], type = "l", xlab = "", ylab = "")
#' abline(h = lvl, col = "blue")
#' lines(t, trig$v[, 2], col = "red", lwd = 2)
#' par(op)
#'
#' @seealso \code{\link{clustersegment}}
#'
#' @author Juan Pablo Carbajal, \email{carbajal@@ifi.uzh.ch}.\cr
#' Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#
#' @export

schtrig <- function(x, lvl, st = NULL) {

  if (!is.numeric(x)) {
    stop("x must be a numeric vector or matrix")
  }

  if (is.vector(x)) {
    x <- matrix(x, ncol = 1)
    vec <- TRUE
  } else if (is.matrix(x)) {
    vec <- FALSE
  } else {
    stop("x must be a numeric vector or matrix")
  }
  nc <- ncol(x)
  nr <- nrow(x)

  if (!is.numeric(lvl) || !is.vector(lvl)) {
    stop("lvl must be a numeric vector")
  }
  if (length(lvl) == 1) {
    lvl <- abs(lvl) * c(1, -1)
  } else {
    lvl <- sort(lvl, decreasing = TRUE)
  }

  if (is.null(st)) {
    st <- rep(0, nc)
  } else if (!is.vector(st) || length(st) != nc) {
    stop(paste("st must be NULL a vector of length", nc))
  }

  v      <- matrix(NA, nr, nc)
  v[1, ] <- st

  ## Signal is above up level
  up    <- x > lvl[1]
  v[up] <- 1

  ## Signal is below down level
  dw    <- x < lvl[2]
  v[dw] <- 0

  ## Resolve intermediate states
  ## Find data between the levels
  idx    <- is.na(v)
  rng <- clustersegment(t(idx))
  if (!is.list(rng)) {
    crng <- list(rng)
  } else {
    crng <- rng
  }

  if (length(crng) > 0 && nrow(crng[[1]]) == 2) {

    for (i in seq_len(nc)) {
      ## Record the state at the beginning of the interval between levels
      prev           <- crng[[i]][1, ] - 1
      prev[prev < 1] <- 1
      st             <- v[prev, i]

      ## Copy the initial state to the interval
      ini_idx <- crng[[i]][1, ]
      end_idx <- crng[[i]][2, ]
      for (j in seq_along(ini_idx)) {
        v[ini_idx[j]:end_idx[j], i] <- st[j]
      }
    }
  }
  st <- v[nr, ]

  if (vec) {
    v <- as.vector(v)
  }

  list(v = v, rng = rng, st = st)
}
