# sgolay.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Original Matlab/Octave version:
# Copyright (C) 2001 Paul Kienzle <pkienzle@users.sf.net>
# Copyright (C) 2004 Pascal Dupuis <Pascal.Dupuis@esat.kuleuven.ac.be>
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
# 20200322    GvB       setup for gsignal v0.1.0
#------------------------------------------------------------------------------

#' Savitzky-Golay filter design
#'
#' Compute the filter coefficients for all Savitzky-Golay FIR smoothing filters.
#'
#' The early rows of the resulting filter smooth based on future values and
#' later rows smooth based on past values, with the middle row using half future
#' and half past.  In particular, you can use row \code{i} to estimate
#' \code{x(k)} based on the \code{i-1} preceding values and the \code{n-i}
#' following values of \code{x} values as \code{y(k) = F[i, ] *
#' x[(k - i + 1):(k + n -i)]}.
#'
#' Normally, you would apply the first \code{(n-1)/2} rows to the first \code{k}
#' points of the vector, the last \code{k} rows to the last \code{k} points of
#' the vector and middle row to the remainder, but for example if you were
#' running on a real-time system where you wanted to smooth based on the all the
#' data collected up to the current time, with a lag of five samples, you could
#' apply just the filter on row \code{n - 5} to your window of length \code{n}
#' each time you added a new sample.
#'
#' @param p Polynomial filter order; must be smaller than \code{n}.
#' @param n Filter length; must a an odd positive integer.
#' @param m Return the m-th derivative of the filter coefficients. Default: 0
#' @param ts Scaling factor. Default: 1
#'
#' @return An square matrix with dimensions \code{length(n)} that is of class
#'   \code{"sgolayFilter"}, so it can be used with \code{filter}.
#'
#' @examples
#' ## Generate a signal that consists of a 0.2 Hz sinusoid embedded
#' ## in white Gaussian noise and sampled five times a second for 200 seconds.
#' dt <- 1 / 5 
#' t <- seq(0, 200 - dt, dt)
#' x <- 5 * sin(2 * pi * 0.2 * t) + rnorm(length(t))
#' ## Use sgolay to smooth the signal.
#' ## Use 21-sample frames and fourth order polynomials.
#' p <- 4
#' n <- 21
#' sg <- sgolay(p, n)
#' ## Compute the steady-state portion of the signal by convolving it
#' ## with the center row of b.
#' ycenter <- conv(x, sg[(n + 1)/2, ], 'valid')
#' ## Compute the transients. Use the last rows of b for the startup
#' ## and the first rows of b for the terminal.
#' ybegin <- sg[seq(nrow(sg), (n + 3) / 2, -1), ] %*% x[seq(n, 1, -1)]
#' yend <- sg[seq((n - 1)/2, 1, -1), ] %*%
#'         x[seq(length(x), (length(x) - (n - 1)), -1)]
#' ## Concatenate the transients and the steady-state portion to
#' ## generate the complete smoothed signal.
#' ## Plot the original signal and the Savitzky-Golay estimate.
#' y = c(ybegin, ycenter, yend)
#' plot(t, x, type = "l", xlab = "", ylab = "", ylim = c(-8, 10))
#' lines(t, y, col = 2)
#' legend("topright", c('Noisy Sinusoid','S-G smoothed sinusoid'),
#'   lty = 1, col = c(1,2))
#'
#' @seealso \code{\link{sgolayfilt}}
#'
#' @author Paul Kienzle \email{pkienzle@@users.sf.net},\cr
#'  Pascal Dupuis, \email{Pascal.Dupuis@@esat.kuleuven.ac.be}.\cr
#'  Conversion to R Tom Short,\cr
#'  adapted by Geert van Boxtel \email{G.J.M.vanBoxtel@@gmail.com}.
#'
#' @export

sgolay <- function(p, n, m = 0, ts = 1) {

  if (!(isPosscal(p) && isWhole(p))) {
    stop("p must be a positive integer")
  }
  if (!(isPosscal(n) && isWhole(n) && n %% 2 == 1 && n > p)) {
    stop("n must be an odd positive integer > p")
  }

  Fm <- matrix(0L, n, n)
  k <- floor(n / 2)
  for (row in 1:(k + 1)) {
    ## Construct a matrix of weights Cij = xi ^ j.  The points xi are
    ## equally spaced on the unit grid, with past points using negative
    ## values and future points using positive values.
    Ce <- (((1:n) - row) %*% matrix(1, 1, p + 1)) ^ (matrix(1, n) %*% (0:p))
    ## A = pseudo-inverse (C), so C*A = I; this is constructed from the SVD
    A <- pracma::pinv(Ce, tol = .Machine$double.eps)
    ## Take the row of the matrix corresponding to the derivative
    ## you want to compute.
    Fm[row, ] <- A[(1 + m), ]
  }
  ## The filters shifted to the right are symmetric with those to the left.
  Fm[((k + 2):n), ] <- (-1)^m * Fm[k:1, n:1]
  if (m > 0) {
    Fm <- Fm * prod(1:m) / (ts^m)
  }
  class(Fm) <- "sgolayFilter"
  Fm
}
