# ss2tf.R
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
# 20261102 GvB       setup for gsignal v0.4.0
#------------------------------------------------------------------------------

#' Convert state-space representation to transfer function
#'
#' Convert digital filter transfer function data to state-space form. The
#' state space system:
#' \deqn{x = Ax + Bu}
#' \deqn{y = Cx + Du}
#' is converted to a transfer function:
#' \deqn{G(s) = \frac{b(s)}{a(s)}}
#'
#' @param A State matrix (must be square)
#' @param B Input-to-state matrix (must have dimensions nrow(A) x 1)
#' @param C State-to-output matrix (must have nrow(A) columns)
#' @param D Feedthrough matrix (must have dimensions nrow(C) x 1)
#' @param tol tolerance to remove insignificant numerical noise (default: 1e-10)
#' 
#' @return A list of class "Arma" with the following list elements:
#' \describe{
#'   \item{b}{moving average (MA) polynomial coefficients}
#'   \item{a}{autoregressive (AR) polynomial coefficients}
#' }
#' 
#' 
#' @seealso See also \code{\link{tf2ss}} 
#'
#' @examples
#' A <- rbind(
#'   c(0, 0, 0, 0, -0.9^5),
#'   c(1, 0, 0, 0, 0),
#'   c(0, 1, 0, 0, 0),
#'   c(0, 0, 1, 0, 0),
#'   c(0, 0, 0, 1, 0)
#' )
#' B <- matrix(c(1, 0, 0, 0, 0), ncol = 1)
#' C <- matrix(c(0, 1, 0, 0, 0.3^3), nrow = 1)
#' D <- 0
#' ss2tf(A, B, C, D)
#' 
#' # example with multiple outputs
#' A <- rbind(
#'   c(-0.4,-1),
#'   c(1, 0)
#' )
#' B <- matrix(c(1, 0), ncol = 1)
#' C <- rbind(
#'   c(2, 3),
#'   c(1.6, 0)
#' )
#' D <- rbind(0, 1)
#' ss2tf(A, B, C, D)
#' 
#' @author Geert van Boxtel, \email{gjmvanboxtel@@gmail.com}.
#'
#' @export

ss2tf <- function(A, B, C, D, tol = 1e-10) {
  
  A <- as.matrix(A)
  B <- as.matrix(B)
  C <- as.matrix(C)
  D <- as.matrix(D)
  
  storage.mode(A) <- "double"
  storage.mode(B) <- "double"
  storage.mode(C) <- "double"
  storage.mode(D) <- "double"
  
  n <- nrow(A)
  nout <- nrow(C)
  
  # Validate dimensions
  if (n != ncol(A)) {
    stop("A must be square")
  }
  if (!all(dim(B) == c(n, 1L))) {
    stop(paste("B must have dimensions", n, "x 1"))
  }
  if (ncol(C) != n) {
    stop(paste("C must have", n, "columns"))
  }
  if (!all(dim(D) == c(nout, 1L))) {
    stop(paste("D must have dimensions", nout, "x 1"))
  }
  
  if (any(!is.finite(A)) ||
      any(!is.finite(B)) ||
      any(!is.finite(C)) ||
      any(!is.finite(D)))
    stop("All matrix elements must be finite.")
  
  # Zero-order system
  if (n == 0L) {
    return(list(
      num = matrix(as.numeric(D), nrow = nout),
      den = 1
    ))
  }
  
  # --------------------------------------------------
  # Denominator: characteristic polynomial det(zI - A)
  # --------------------------------------------------
  
  den <- numeric(n + 1L)
  den[1L] <- 1
  
  traces <- numeric(n)
  Apower <- diag(n)
  
  for (k in seq_len(n)) {
    Apower <- Apower %*% A
    traces[k] <- sum(diag(Apower))
    
    den[k + 1L] <- -sum(
      den[k:1L] * traces[seq_len(k)]
    ) / k
  }
  
  # --------------------------------------------------
  # Markov parameters:
  # H[k] = C A^(k-1) B
  # One row per output
  # --------------------------------------------------
  
  h <- matrix(0, nrow = nout, ncol = n)
  
  v <- B
  
  for (k in seq_len(n)) {
    h[, k] <- as.numeric(C %*% v)
    v <- A %*% v
  }
  
  # --------------------------------------------------
  # Numerator for each output
  # --------------------------------------------------
  
  num <- matrix(0, nrow = nout, ncol = n + 1L)
  
  for (i in seq_len(nout)) {
    
    num[i, ] <- as.numeric(D[i, 1L]) * den
    
    for (k in seq_len(n)) {
      
      num[i, k + 1L] <- num[i, k + 1L] +
        sum(
          den[seq_len(k)] *
            rev(h[i, seq_len(k)])
        )
    }
  }
  
  # Remove insignificant numerical noise
  num[abs(num) < tol] <- 0
  den[abs(den) < tol] <- 0
  
  Arma(
    b = as.vector(num),
    a = as.vector(as.numeric(den))
  )
}

