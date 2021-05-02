# fwht.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Original Octave code:
# Copyright (C) 2013-2019 Mike Miller
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
# 20201023  GvB       setup for gsignal v0.1.0
#------------------------------------------------------------------------------

#' Fast Walsh-Hadamard Transform
#'
#' Compute the (inverse) Fast Walsh-Hadamard transform of a signal.
#'
#' @param x input data, specified as a numeric vector or matrix. In case of a
#'   vector it represents a single signal; in case of a matrix each column is a
#'   signal. \code{fwht} operates only on signals with length equal to a power
#'   of 2. If the length of \code{x} is less than a power of 2, its length is
#'   padded with zeros to the next greater power of two before processing.
#' @param n transform length, specified as a positive integer scalar. Default:
#'   \code{NROW(x)}.
#' @param ordering order of the Walsh-Hadamard transform coefficients, one of:
#' \describe{
#'   \item{"sequency"}{(Default) Coefficients in order of increasing sequency
#'   value, where each row has an additional zero crossing.}
#'   \item{"hadamard"}{Coefficients in normal Hadamard order}
#'   \item{"dyadic"}{Coefficients in Gray code order, where a single bit change
#'   occurs from one coefficient to the next}
#' }
#'
#' @return (Inverse) Fast Walsh Hadamard transform, returned as a vector or
#'   matrix.
#'
#' @examples
#' x <- c(19, -1, 11, -9, -7, 13, -15, 5)
#' X <- fwht(x)
#' all.equal(x, ifwht(X))
#'
#' @author Mike Miller.\cr
#' Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#'
#' @references \url{https://en.wikipedia.org/wiki/Hadamard_transform}
#' @references \url{https://en.wikipedia.org/wiki/Fast_Walsh-Hadamard_transform}
#'
#' @rdname fwht
#' @export

ifwht <- function(x, n = NROW(x),
                  ordering = c("sequency", "hadamard", "dyadic")) {

  # check parameters
  if (!(is.vector(x) || is.matrix(x)) || !is.numeric(x)) {
    stop("x must be a numeric or vector or matrix")
  }

  if (is.vector(x)) {
    vec <- TRUE
    x <- as.matrix(x, ncol = 1)
  } else {
    vec <- FALSE
  }
  nr <- nrow(x)

  if (!isPosscal(n) || !isWhole(n)) {
    stop("n must be a positive integer")
  }
  # force n to be a power of 2
  n <- nextpow2(n)
  if (n != nr) {
    x <- postpad(x, n)
  }

  ordering <- match.arg(ordering)

  # Zero-based index for normal Hadamard ordering
  idx <- seq(0, n - 1)
  nbits <- floor(log2(max(idx))) + 1   # number of significant bits

  # Gray code permutation of index for alternate orderings
  idx_bin <- matrix(0, n, nbits)
  if (ordering == "dyadic") {
    for (i in seq_len(n)) {
      idx_bin[i, ] <-  as.integer(intToBits(idx[i]))[1:nbits]
      idx[i] <- bin2dec(paste(idx_bin[i, ], collapse = "")) + 1
    }
  } else if (ordering == "sequency") {
    for (i in seq_len(n)) {
      idx_bin[i, ] <- as.integer(rev(intToBits(idx[i])))[(32 - nbits + 1):32]
      idx_bin_a <- idx_bin[i, 1:(nbits - 1)]
      idx_bin_b <- idx_bin[i, 2:nbits]
      idx_bin[i, 2:nbits] <- (idx_bin_a + idx_bin_b) %% 2
      idx[i] <- bin2dec(paste(rev(idx_bin[i, ]), collapse = "")) + 1
    }
  } else {
    idx <- idx + 1
  }

  # do the transform
  if (n < 2) {
    y <- x
  } else {
    y <- .Call("_gsignal_fwht", PACKAGE = "gsignal", x)
  }

  # apply ordering
  y <- y[idx, ]

  # cleanup and exit
  if (vec) {
    y <- as.vector(y)
  }
  y
}

#' @rdname fwht
#' @export

 fwht <- function(x, n = NROW(x),
                  ordering = c("sequency", "hadamard", "dyadic")) {

   if (!isPosscal(n) || !isWhole(n)) {
     stop("n must be a positive integer")
   }
   n <- nextpow2(n)
   y <- ifwht(x, n, ordering)
   y <- y / n
   y
 }
