# filter_zi.R
# Copyright (C) 2021 Geert van Boxtel <gjmvanboxtel@gmail.com>
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
# 20210319  GvB       setup for gsignal v0.3.0
#------------------------------------------------------------------------------

#' Filter initial conditions
#'
#' Construct initial conditions for a filter
#'
#' This function computes an initial state for the filter function that
#' corresponds to the steady state of the step response. In other words, it
#' finds the initial condition for which the response to an input of all ones is
#' a constant. Therefore, the results returned by this function can also be
#' obtained using the function \code{\link{filtic}} by setting \code{x} and
#' \code{y} to all 1's (see the examples).
#'
#' A typical use of this function is to set the initial state so that the output
#' of the filter starts at the same value as the first element of the signal to
#' be filtered.
#'
#' @param filt For the default case, the moving-average coefficients of an ARMA
#'   filter (normally called ‘b’), specified as a vector.
#' @param a the autoregressive (recursive) coefficients of an ARMA filter,
#'   specified as a vector.
#' @param ... additional arguments (ignored).
#'
#' @return The initial state for the filter, returned as a vector.
#'
#' @examples
#' ## taken from Python scipy.signal.lfilter_zi documentation
#'
#' h <- butter(5, 0.25)
#' zi <- filter_zi(h)
#' y <- filter(h, rep(1, 10), zi)
#' ## output is all 1, as expected.
#' y2 <- filter(h, rep(1, 10))
#' ## if the zi argument is not given, the output
#' ## does not return the final conditions
#'
#' x <- c(0.5, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0)
#' y <- filter(h, x, zi = zi*x[1])
#' ## Note that the zi argument to filter was computed using
#' ## filter_zi and scaled by x[1]. Then the output y has no
#' ## transient until the input drops from 0.5 to 0.0.
#'
#' ## obtain the same results with filtic
#' lab <- max(length(h$b), length(h$a)) - 1
#' ic <- filtic(h, rep(1, lab), rep(1, lab))
#' all.equal(zi, ic)
#'
#' @author Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com},
#'  converted to R from Python scipy.signal.lfilter_zi.
#'
#' @seealso \code{\link{filtic}}
#'
#' @references Gustafsson, F. (1996). Determining the initial states in
#'   forward-backward filtering. IEEE Transactions on Signal Processing, 44(4),
#'   988 - 992.
#'
#' @rdname filter_zi
#' @export

filter_zi <- function(filt, ...) UseMethod("filter_zi")

#' @rdname filter_zi
#' @method filter_zi default
#' @export

filter_zi.default <- function(filt, a, ...) {

  if (!is.vector(filt) || ! is.vector(a) ||
      !is.numeric(filt) || !is.numeric(a)) {
    stop("b and a must be numeric vectors")
  }

  while (length(a) > 1 && a[1] == 0) {
    a <- a[2:length(a)]
  }
  if (length(a) < 1) {
    stop("There must be at least one nonzero element in the vector a")
  }

  if (a[1] != 1) {
    # Normalize the coefficients so a[1] == 1.
    filt <- filt / a[1]
    a <- a / a[1]
  }

  n <- max(length(a), length(filt))

  # Pad a or b with zeros so they are the same length.
  filt <- postpad(filt, n)
  a <- postpad(a, n)

  B <- filt[2:length(filt)] - a[2:length(a)] * filt[1]
  IminusA <- diag(rep(1, n - 1)) - t(pracma::compan(a))
  # Solve zi = A*zi + B
  zi <- solve(IminusA, B)
  zi
}

#' @rdname filter_zi
#' @method filter_zi Arma
#' @export
filter_zi.Arma <- function(filt, ...) # IIR
  filter_zi(filt$b, filt$a, ...)

#' @rdname filter_zi
#' @method filter_zi Ma
#' @export
filter_zi.Ma <- function(filt, ...) # FIR
  filter_zi(unclass(filt), 1, ...)

#' @rdname filter_zi
#' @method filter_zi Sos
#' @export
filter_zi.Sos <- function(filt, ...) { # Second-order sections
  if (filt$g != 1) {
    filt$sos[1, 1:3] <- filt$sos[1, 1:3] * filt$g
  }
  L <- NROW(filt$sos)
  zi <- matrix(0, L, 2)
  scale <- 1.0
  for (l in seq_len(L)) {
    b <- filt$sos[l, 1:3]
    a <- filt$sos[l, 4:6]
    zi[l, ] <- scale * filter_zi.default(b, a)
    # If H(z) = B(z)/A(z) is this section's transfer function, then
    # b.sum()/a.sum() is H(1), the gain at omega=0.  That's the steady
    # state value of this section's step response.
    scale <- scale *  sum(b) / sum(a)
  }
  zi
}

#' @rdname filter
#' @method filter Zpg
#' @export
filter_zi.Zpg <- function(filt, x, ...) # zero-pole-gain form
  filter_zi(as.Arma(filt), ...)
