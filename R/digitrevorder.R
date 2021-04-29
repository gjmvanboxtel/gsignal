# digitrevorder.R
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
# 20200821  GvB       setup for gsignal v0.1.0
#------------------------------------------------------------------------------

#' Permute input to digit-reversed order
#'
#' Reorder the elements of the input vector in digit-reversed order.
#'
#' This function is useful for pre-ordering a vector of filter coefficients for
#' use in frequency-domain filtering algorithms, in which the fft and ifft
#' transforms are computed without digit-reversed ordering for improved run-time
#' efficiency.
#'
#' @param x input data, specified as a vector. The length of \code{x} must be an
#'   integer power of \code{r}.
#' @param r radix base used for the number conversion, which can be any integer
#'   from 2 to 36. The elements of \code{x} are converted to radix \code{r} and
#'   reversed.
#' @param index.return logical indicating if the ordering index vector should be
#'   returned as well. Default \code{FALSE}.
#'
#' @return The digit-reversed input vector. If \code{index.return = TRUE}, then
#'   a list containing the digit-reversed input vector (\code{y}, and the
#'   digit-reversed indices (\code{i}).
#'
#' @examples
#'
#' res <- digitrevorder(0:8, 3)
#'
#' @author Mike Miller.\cr
#'   Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#'
#' @seealso \code{\link{bitrevorder}}, \code{\link{fft}}, \code{\link{ifft}}
#'
#' @export

digitrevorder <- function(x, r, index.return = FALSE) {

  if (!is.vector(x)) {
    stop("x must be a vector")
  } else if (!isPosscal(r) || !isWhole(r) || !(r >= 2 && r <= 36)) {
    stop("r must be an integer between 2 and 36")
  } else {
    tmp <- log(length(x)) / log(r)
    if (trunc(tmp) != tmp) {
      stop(paste("x must have length equal to an integer power of", r))
    }
  }
  if (!is.logical(index.return)) {
    stop("index.return must be TRUE or FALSE")
  }

  old_ind <- seq(0, length(x) - 1)
  new_ind <- gdec2base(old_ind, r)

  i <- new_ind + 1
  y <- rep(0L, length(x))
  y[old_ind + 1] <- x[i]

  if (index.return) {
    retval <- list(y, i)
  } else {
    retval <- y
  }
  retval
}
# R version of Octave function dec2base, simplified and adapted for use with
# digitrevorder, not exported to the namespace)

## Original Octave dec2base:
## Author: Daniel Calvelo <dcalvelo@yahoo.com>
## Adapted-by: Paul Kienzle <pkienzle@kienzle.powernet.co.uk>

gdec2base <- function(d, base, len = 0) {

  # better safe than sorry
  d <- as.vector(round(abs(as.numeric(d))))

  symbols <- c(as.character(0:9), LETTERS)
  if (is.character(base)) {
    symbols <- unique(unlist(strsplit(gsub("[[:space:]]", "", base), "")))
    base <- length(symbols)
  } else if (!isScalar(base)) {
    base <- base[1]
  } else if (base < 2 || base > length(symbols)) {
    base <- max(min(base, length(symbols)), 2)
  }

  ## determine number of digits required to handle all numbers, can overflow
  ## by 1 digit
  max_len <- round(log(max(max(d), 1)) / log(base)) + 1
  max_len <- max(max_len, len)

  ## determine digits for each number
  digits <- matrix(0L, length(d), max_len)
  for (k in seq(max_len, 1, -1)) {
    digits[, k] <- d %% base
    d <- round((d - digits[, k]) / base)
  }

  ## convert digits to symbols
  retval <- matrix(symbols[digits + 1],
                   nrow = NCOL(digits), ncol = NROW(digits),
                   byrow = TRUE)

  ## Check if the first element is the zero symbol.  It seems possible
  ## that LEN is provided, and is less than the computed MAX_LEN and
  ## MAX_LEN is computed to be one larger than necessary, so we would
  ## have a leading zero to remove.  But if LEN >= MAX_LEN, we should
  ## not remove any leading zeros.
  if ((len == 0 || (len != 0 && max_len > len))
      && NROW(retval) != 1 && !any(retval[1, ] != symbols[1])) {
    retval <- retval[-1, ]
  }

  # GvB: flip and convert back to numeric if retval is a matrix
  nc <- NCOL(retval); nr <- NROW(retval)
  if (nc > 1) {
    if (nr > 1) {
      retval <- pracma::flipud(retval)
    }
    tmp <- array("", nc)
    for (k in seq_len(nc)) {
      tmp[k] <- paste0(retval[, k], collapse = "")
    }
    retval <- tmp
  }
  strtoi(retval, base)
}
