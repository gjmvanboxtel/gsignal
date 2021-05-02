# cplxpair.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
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
# 20200327  GvB       setup for gsignal v0.1.0
# 20210405  GvB       changed 'dim' argument to MARGIN
# 20210506  GvB       bugfix in Check if real parts occur in pairs.
#------------------------------------------------------------------------------

#' Complex conjugate pairs
#'
#' Sort complex numbers into complex conjugate pairs ordered by increasing real
#' part.
#'
#' The negative imaginary complex numbers are placed first within each pair. All
#' real numbers (those with \code{abs(Im (z) / z) < tol)} are placed after the
#' complex pairs.
#'
#' An error is signaled if some complex numbers could not be paired and if all
#' complex numbers are not exact conjugates (to within \code{tol}).
#'
#' @note There is no defined order for pairs with identical real parts but
#'   differing imaginary parts.
#'
#' @param z Vector, matrix, or array of complex numbers.
#' @param tol Weighting factor \code{0 < tol < 1}, which determines the
#'   tolerance of matching. Default: \code{100 * .Machine$double.eps}. (This
#'   definition differs from the 'Octave' usage).
#' @param MARGIN Vector giving the subscripts which the function will be applied
#'   over. E.g., for a matrix 1 indicates rows, 2 indicates columns, c(1, 2)
#'   indicates rows and columns. Where X has named dimnames, it can be a
#'   character vector selecting dimension names. Default: 2 (columns).
#'
#' @return Vector, matrix or array containing ordered complex conjugate pairs by
#'   increasing real parts.
#'
#' @examples
#' r <- rbind(t(cplxpair(exp(2i * pi * 0:4 / 5))),
#'            t(exp(2i * pi *c(3, 2, 4, 1, 0) / 5)))
#'
#' @author Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#'
#' @seealso \code{\link{cplxreal}}
#'
#' @export

cplxpair <- function(z, tol = 100 * .Machine$double.eps, MARGIN = 2) {

  vec <- FALSE
  if (is.vector(z)) {
    vec <- TRUE
    z <- as.matrix(z)
  }

  if (!isPosscal(tol) || tol > 1) {
    stop("tol must be a positive scalar between 0 and 1")
  }

  d <- dim(z)
  if (any(d <= 0)) {
    y <- array(0L, d)
    if (vec) y <- as.vector(y)
    return(y)
  }

  sort_vec <- function(v) {

    v <- as.vector(v)
    l <- length(v)
    y <- rep(0L, l)

    # Find real values and put them (sorted) at the end of y
    idx <- which(abs(Im(v)) <= tol * abs(v))
    n <- length(idx)
    if (n > 0) {
      y[(l - n + 1):l] <- sort(Re(v[idx]))
      v <- v[-idx]
    }

    # remaining values are all complex, if any
    nv <- length(v)
    if (nv > 0) {

      if (nv %% 2 == 1) {
        stop("Could not pair all complex numbers")
      }

      # Sort v based on real part
      s <- sort(Re(v), index.return = TRUE)
      v <- v[s$ix]

      # Check if real parts occur in pairs. If not: error
      a <- matrix(s$x, ncol = 2, byrow = TRUE)
      if (any(abs(a[, 1] - a[, 2]) > tol * abs(a[, 1]))) {
        stop("Could not pair all complex numbers")
      }

      # Check if imag part of real part pairs are conjugates
      yix <- 1
      while (length(v) > 0) {

        # Find all real parts equal to to real(v[1])
        idx <- which(abs(Re(v) - Re(v[1])) <= tol * abs(Re(v)))
        nn <- length(idx)
        if (nn <= 1) {
          stop("Could not pair all complex numbers")
        }

        # Sort the imag parts of those values
        si <- sort(Im(v[idx]), index.return = TRUE)
        q <- v[si$ix]    # Get values with identical real parts,
        lq <- length(q)  # now sorted by imaginary parts

        # Verify conjugate-pairing of imag parts
        if (any(abs(si$x + rev(si$x)) > tol * abs(q))) {
          stop("Could not pair all complex numbers")
        }

        # Keep value with positive imag part, and compute conjugate
        # Value with smallest neg imag part first, then its conj
        y[yix:(yix + nn - 1)] <-
          as.vector(t(cbind(Conj(q[seq(lq, (nn / 2 + 1), -1)]),
                            q[seq(lq, (nn / 2 + 1), -1)])))

        yix <- yix + nn # update y index
        v <- v[-idx]    # Remove entries from v
      }
    }
    y
  }

  if (vec) {
    y <- as.vector(sort_vec(z))
  } else {
    y <- apply(z, MARGIN, sort_vec)
  }
  y
}
