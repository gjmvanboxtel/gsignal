# clustersegment.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Original Octave code:
# Copyright (C) 2010 Juan Pablo Carbajal <carbajal@ifi.uzh.ch>
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
# See also: http://www.gnu.org/licenses/gpl-2.0.txt
#
# Version history
# 20201122  GvB       setup for gsignal v0.1.0
#---------------------------------------------------------------------------------------------------------------------

#' Cluster Segments
#' 
#' Calculate boundary indexes of clusters of 1’s.
#' 
#' The function calculates the initial index and end index of the sequences of
#' 1’s in the rows of \code{x}. The clusters are sought in the rows of the array
#' \code{x}. The function works by finding the indexes of jumps between
#' consecutive values in the rows of \code{x}.
#'   
#' @param x input data, specified as a numeric vector or matrix.
#' 
#' @return a \code{\link{list}} of matrices size \code{nr}, where \code{nr} is
#'   the number of rows in \code{x}. Each element of the list contains a matrix
#'   with two rows. The first row is the initial index of a sequence of 1’s and
#'   the second row is the end index of that sequence. If \code{nr == 1} the
#'   output is a matrix with two rows.
#' 
#' @examples
#' (x <- c(0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1))
#' # [1] 0 0 1 1 1 0 0 1 0 0 0 1 1
#' (ranges <- clustersegment(x))
#' #       [,1] [,2] [,3]
#' # [1,]    3    8   12
#' # [2,]    5    8   13
#' # The first sequence of 1's in x lies in the interval
#' ranges[1,1]:ranges[2,1]
#' # [1] 3 4 5
#' 
#' x <- matrix(as.numeric(runif(30) > 0.4), 3, 10)
#' ranges <- clustersegment(x)
#'
#' @author Juan Pablo Carbajal \email{carbajal@@ifi.uzh.ch}; port to R by Geert
#'   van Boxtel \email{G.J.M.vanBoxtel@@gmail.com}.
#
#' @export

clustersegment <- function (x) {

  if (is.vector(x)) {
    nr <- 1
    nc <- length(x)
    x <- matrix(x, nrow = 1)
    vec <- TRUE
  } else if (is.matrix(x)) {
    nc <- ncol(x)
    nr <- nrow(x)
    vec <- FALSE
  } else {
    stop ('x must be a numeric vector or matrix')
  }
  
  bool_discon <- t(apply(x, 1, diff))
  y <- list()
  for (i in seq_len(nr)) {
    idxUp <- which(bool_discon[i, ] > 0) + 1
    idxDwn <- which(bool_discon[i, ] < 0)
    tLen <- length(idxUp) + length(idxDwn)
    
    contRange <- matrix(0, nrow = 1, ncol = tLen)
    if (x[i, 1] == 1) {
      ## first event was down
      contRange <- c(contRange, 0)
      contRange[1] <- 1
      contRange[seq(2, tLen + 1, 2)] <- idxDwn
      contRange[seq(3, tLen + 1, 2)] <- idxUp
    } else {
      ## first event was up
      contRange[seq(1, tLen, 2)] <- idxUp
      contRange[seq(2, tLen, 2)] <- idxDwn
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
  
  if (nr == 1) {
    y <- as.matrix(y[[1]])
  }
  y  
}
