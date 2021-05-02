# sigmoid_train.R
# Copyright (C) 2019 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Matlab/Octave signal package:
# Copyright (C) 2011-2013 Juan Pablo Carbajal <carbajal@ifi.uzh.ch>
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
# 20191204 Geert van Boxtel          First version for v0.1.0
# 20200322 Geert van Boxtel          used NROW and NCOL; expand rc
#------------------------------------------------------------------------------

#' Sigmoid Train
#'
#' Evaluate a train of sigmoid functions at \code{t}.
#'
#' The number and duration of each sigmoid is determined from ranges. Each row
#' of \code{ranges} represents a real interval, e.g. if sigmoid \code{i} starts
#' at \code{t = 0.1} and ends at \code{t = 0.5}, then \code{ranges[i, ] = c(0.1,
#' 0.5)}. The input \code{rc} is an array that defines the rising and falling
#' time constants of each sigmoid. Its size must equal the size of ranges.
#'
#' The individual sigmoids are returned in \code{s}. The combined sigmoid train
#' is returned in the vector \code{y} of length equal to \code{t}, and such that
#' \code{y = max(s)}.
#'
#' @param t Vector (or coerced to a vector) of time values at which the sigmoids
#'   are calculated.
#' @param ranges Matrix or array with 2 columns containing the time values
#'   within \code{t} at which each sigmoid is evaluated. The number of sigmoids
#'   is determined by the number of rows in \code{ranges}.
#' @param rc Time constant. Either a scalar or a matrix or array with 2 columns
#'   containing the rising and falling time constants of each sigmoid. If a
#'   matrix or array is passed in \code{rc}, its size must equal the size of
#'   \code{ranges}. If a single scalar is passed in \code{rc}, then all sigmoids
#'   have the same time constant and are symmetrical.
#'
#' @return A list consisting two variables; \code{y} the combined sigmoid train
#'   (length identical to \code{t}), and \code{s}, the individual sigmoids
#'   (number of rows equal to number of rows in \code{ranges} and \code{rc}.
#'
#' @examples
#'
#' t <- seq(0, 2, length.out = 500)
#' ranges <- rbind(c(0.1, 0.4), c(0.6, 0.8), c(1, 2))
#' rc <- rbind(c(1e-2, 1e-3), c(1e-3, 2e-2), c(2e-2, 1e-2))
#' st <- sigmoid_train (t, ranges, rc)
#' plot(t, st$y[1,], type="n", xlab = "Time(s)", ylab = "S(t)",
#'      main = "Vectorized use of sigmoid train")
#' for (i in 1:3) rect(ranges[i, 1], 0, ranges[i, 2], 1,
#'                     border = NA, col="pink")
#' for (i in 1:3) lines(t, st$y[i,])
#' # The colored regions show the limits defined in range.
#'
#' t <- seq(0, 2, length.out = 500)
#' ranges <- rbind(c(0.1, 0.4), c(0.6, 0.8), c(1, 2))
#' rc <- rbind(c(1e-2, 1e-3), c(1e-3, 2e-2), c(2e-2, 1e-2))
#' amp <- c(4, 2, 3)
#' st <- sigmoid_train (t, ranges, rc)
#' y <- amp %*% st$y
#' plot(t, y[1,], type="l", xlab = 'time', ylab = 'signal',
#'      main = 'Varying amplitude sigmoid train', col="blue")
#' lines(t, st$s, col = "orange")
#' legend("topright", legend = c("Sigmoid train", "Components"),
#'        lty = 1, col = c("blue", "orange"))
#'
#' @author Juan Pablo Carbajal, \email{carbajal@@ifi.uzh.ch}.\cr
#' Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#
#' @export

sigmoid_train <- function(t, ranges, rc) {

  t <- as.vector(t)

  ## number of sigmoids
  if (is.vector(ranges)) {
    ranges <- as.matrix(t(ranges))
  }
  nr <- NROW(ranges)
  nc <- NCOL(ranges)
  if (is.null(nr) || nr <= 0)
    stop("ranges must be a vector, or an array or matrix with at least 1 row")
  if (is.null(nc) || nc != 2)
    stop("ranges must be a vector, or an array or matrix with 2 columns")

  ## Parse time constants
  if (isScalar(rc)) {
    # All sigmoids have the same time constant and are symmetric
    rc <- rc * matrix(1L, nr, 2)
  } else if (is.vector(rc)) {
    rc <- as.matrix(t(rc))
  } else if ((nrow(rc) == 1 || ncol(rc) == 1) && nr > 1) {
    # All sigmoids have different time constants but are symmetric
    if (nrow(rc) == 1) {
      rc <- t(rc)
    }
    if (nrow(rc) != nr)
      stop("length of time constant must equal number of ranges")
    rc <- cbind(rc, rc)
  }

  a_up <- apply(t(apply(t(ranges[, 1]), 2,
                        function(x) x - t)), 2, function(x) x / rc[, 1])
  a_dw <- apply(t(apply(t(ranges[, 2]), 2,
                        function(x) x - t)), 2, function(x) x / rc[, 2])

  ## Evaluate the sigmoids and mix them
  y <- 1 / (1 + exp(a_up)) * (1 - 1 / (1 + exp(a_dw)))
  if (nr == 1) {
    s <- y
  } else {
    s <- apply(y, 2, max)
  }

  list(y = y, s = s)
}
