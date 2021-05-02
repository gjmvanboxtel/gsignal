# clustersegment.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Original Octave code:
# Copyright (C) 2010 Juan Pablo Carbajal <carbajal@ifi.uzh.ch>
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
# 20201122  GvB       setup for gsignal v0.1.0
# 20201127  GvB       also accept numeric values other than  0 and 1
#------------------------------------------------------------------------------

#' Cluster Segments
#'
#' Calculate boundary indexes of clusters of 1’s.
#'
#' The function calculates the initial index and end index of sequences of 1's
#' rising and falling phases of the signal in \code{x}. The clusters are sought
#' in the rows of the array \code{x}. The function works by finding the indexes
#' of jumps between consecutive values in the rows of \code{x}.
#'
#' @param x input data, specified as a numeric vector or matrix, coerced to
#'   contain only 0's and 1's, i.e., every nonzero element in \code{x} will
#'   be replaced by 1.
#'
#' @return A list of size \code{nr}, where \code{nr} is the number
#'   of rows in \code{x}. Each element of the list contains a matrix with two
#'   rows. The first row is the initial index of a sequence of 1’s and the
#'   second row is the end index of that sequence.
#'
#' @examples
#' (x <- c(0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1))
#' (ranges <- clustersegment(x))
#' # The first sequence of 1's in x lies in the interval
#' (r <- ranges[1,1]:ranges[2,1])
#'
#' x <- matrix(as.numeric(runif(30) > 0.4), 3, 10)
#' ranges <- clustersegment(x)
#'
#' x <- c(0, 1.2, 3, -8, 0)
#' ranges <- clustersegment(x)
#'
#' @author Juan Pablo Carbajal, \email{carbajal@@ifi.uzh.ch}.\cr
#'  Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#
#' @export

clustersegment <- function(x) {

  if (!(is.numeric(x) || is.logical(x) || is.complex(x))) {
    stop("x must be numeric, logical or complex")
  }

  if (is.vector(x)) {
    x <- matrix(x, nrow = 1)
  } else if (!is.matrix(x)) {
    stop("x must be a vector or matrix")
  }
  nc <- ncol(x)
  nr <- nrow(x)

  # coerce to 0's and 1's
  x <- apply(x, c(1, 2), function(x) as.integer(as.logical(x)))

  y <- list()
  for (i in seq_len(nr)) {
    bool_discon <- diff(x[i, ])
    idxUp <- which(bool_discon > 0) + 1L
    idxDwn <- which(bool_discon < 0)
    tLen <- length(idxUp) + length(idxDwn)

    if (tLen <= 0) {
      y[[i]] <- matrix(NA, 1, 1)
    } else  {
      contRange <- rep(0, tLen)
      if (x[i, 1] == 1) {
        ## first event was down
        contRange <- c(contRange, 0)
        contRange[1] <- 1
        if (tLen >= 1) contRange[seq(2, tLen + 1, 2)] <- idxDwn
        if (tLen >= 2) contRange[seq(3, tLen + 1, 2)] <- idxUp
      } else {
        ## first event was up
        contRange[seq(1, tLen, 2)] <- idxUp
        if (tLen >= 2) contRange[seq(2, tLen, 2)] <- idxDwn
      }

      if (x[i, nc] == 1) {
        ## last event was up
        contRange <- c(contRange, nc)
      }

      tLen <- length(contRange)
      if (tLen != 0) {
        dim(contRange) <- c(2, tLen / 2)
        y[[i]] <- contRange
      }
    }
  }

  if (nr == 1 && length(y) > 0) {
    y <- as.matrix(y[[1]])
  }
  y
}
