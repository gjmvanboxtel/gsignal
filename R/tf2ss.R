# tf2ss.R
# Copyright (C) 2026 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Original Octave version:
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
# 20261031 GvB       setup for gsignal v0.4.0
#------------------------------------------------------------------------------

#' Transfer function to state-space form
#'
#' Convert digital filter transfer function data to state-space form. The
#' state space system:
#' \deqn{x = Ax + Bu}
#' \deqn{y = Cx + Du}
#' is obtained from a transfer function:
#' \deqn{G(s) = \frac{b(s)}{a(s)}}
#' 
#' The state space system matrices obtained from this function will be in
#' observable companion form as Wolovich’s Observable Structure Theorem is used.
#' 
#' @note This function will generally result in output that is different from
#'   that in Octave. This is because the Octave implementation uses the
#'   `control` package, which in turn calls SLICOT routine TD04AD to compute the
#'   realization, which is unavailable in R. However, the resulting state-space
#'   equations  describe the same input–output system. This function does
#'   generally return the same results as Matlab's `tf2ss` though.
#'
#' @param num moving average polynomial coefficients (b)
#' @param den autoregressive polynomial coefficients (a)
#'
#' @return A list with the following list elements:
#' \describe{
#'   \item{A}{State matrix}
#'   \item{B}{Input-to-state matrix}
#'   \item{C}{State-to-output matrix}
#'   \item{D}{Feedthrough matrix}
#' }
#'
#' @seealso See also \code{\link{ss2tf}}, \code{\link{tf2sos}},
#'   \code{\link{tf2zp}}
#'
#' @examples
#' b <- c(1, 0, 0, 0.3^3)
#' a <- c(1, 0, 0, 0, 0, 0.9^5)
#' tf2ss(b, a)
#' 
#' # example with multiple outputs
#' num <- rbind(
#'   c(0, 2, 3),
#'   c(1, 2, 1)
#' )
#' den <- c(1, 0.4, 1)
#' ss <- tf2ss(num, den)
#'
#' @author Geert van Boxtel, \email{gjmvanboxtel@@gmail.com}.
#'
#' @export

tf2ss <- function(num, den) {

  # Convert denominator to numeric vector
  den <- as.numeric(den)
  
  # Convert numerator to matrix: one row per output
  if (is.null(dim(num))) {
    num <- matrix(as.numeric(num), nrow = 1L)
  } else {
    num <- as.matrix(num)
    storage.mode(num) <- "double"
  }
  
  # Validate inputs
  if (length(den) == 0L || nrow(num) == 0L || ncol(num) == 0L)
    stop("Numerator and denominator must not be empty.")
  
  if (any(!is.finite(den)) || any(!is.finite(num)))
    stop("All coefficients must be finite.")
  
  # Store original leading denominator coefficient
  gain <- den[1L]
  
  if (gain == 0)
    stop("The leading denominator coefficient must be nonzero.")
  
  # Normalize denominator
  den <- den / gain
  
  # Number of outputs and denominator coefficients
  nout <- nrow(num)
  nden <- length(den)
  
  # Pad numerator rows with leading zeros
  if (ncol(num) > nden)
    stop("Numerator degree must not exceed denominator degree.")
  
  num_pad <- matrix(0, nrow = nout, ncol = nden)
  
  for (i in seq_len(nout)) {
    idx <- (nden - ncol(num) + 1L):nden
    num_pad[i, idx] <- num[i, ]
  }
  
  # Normalize numerator using original leading denominator coefficient
  num_pad <- num_pad / gain
  
  # Number of states
  n <- nden - 1L
  
  # Direct feedthrough matrix
  D <- matrix(num_pad[, 1L], ncol = 1L)
  
  # Zero-order system
  if (n == 0L) {
    return(list(
      A = matrix(numeric(0), nrow = 0L, ncol = 0L),
      B = matrix(numeric(0), nrow = 0L, ncol = 1L),
      C = matrix(numeric(0), nrow = nout, ncol = 0L),
      D = D
    ))
  }
  
  # State matrix A: companion form
  A <- matrix(0, nrow = n, ncol = n)
  A[1L, ] <- -den[-1L]
  
  if (n > 1L) {
    A[cbind(2L:n, 1L:(n - 1L))] <- 1
  }
  
  # Input matrix B
  B <- matrix(0, nrow = n, ncol = 1L)
  B[1L, 1L] <- 1
  
  # Output matrix C
  C <- num_pad[, -1L, drop = FALSE] -
    D %*% matrix(den[-1L], nrow = 1L)
  
  # Return state-space realization
  list(
    A = A,
    B = B,
    C = C,
    D = D
  )
}
