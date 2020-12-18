# movingrms.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Original Matlab/Octave version:
# Copyright (C) 2012 Juan Pablo Carbajal <carbajal@ifi.uzh.ch>
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
#---------------------------------------------------------------------------------------------------------------------

#' Moving Root Mean Square
#'
#' Compute the moving root mean square (RMS) of the input signal.
#'
#' The signal is convoluted against a sigmoid window of width \code{w} and
#' risetime \code{rc}. The units of these parameters are relative to the value
#' of the sampling frequency given in \code{fs}.
#'
#' @param x Input signal, specified as a numeric vector or matrix. In case of a
#'   matrix, the function operates along the columns
#' @param width width of the sigmoid window, in units relative to \code{fs}.
#'   Default: 0.1
#' @param rc Rise time (time constant) of the sigmoid window, in units relative
#'   to \code{fs}. Default: 1e-3
#' @param fs Sampling frequency. Deafult: 1
#'
#' @return A \code{\link{list}} containg 2 variables:
#' \describe{
#'   \item{rmsx}{Output signal with the same dimensions as \code{x}}
#'   \item{w}{Window, returned as a vector}
#' }
#'
#' @examples
#' N <- 128
#' t <- seq(0, 1, length.out = N)
#' x <- sigmoid_train(t, c(0.4, Inf), 1e-2)$y * (2 * runif(length(t)) - 1)
#' fs <- 1 / diff(t[1:2])
#' width <- 0.05
#' rc <- 5e-3
#' ret <- movingrms(as.numeric(scale(x)), width, rc, fs)
#' plot(t, x, type = "l", col = "red", xlab = "", ylab = "")
#' lines(t, ret$rmsx, lwd = 4, col = "black")
#' polygon(c(0, t, length(t)), c(0, ret$rmsx, 0), col = "blue")
#' lines (t, ret$w, lwd = 2, col = "green")
#' legend("topleft", c("data", "window", "movingrms"), lty = 1,
#' col = c("red", "green", "blue"))
#'
#' @seealso \code{\link{sigmoid_train}}
#'
#' @author Juan Pablo Carbajal, \email{carbajal@@ifi.uzh.ch}.\cr
#' Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#'
#' @export

movingrms <- function(x, width = 0.1, rc = 1e-3, fs = 1) {

  if (is.vector(x)) {
    N <- length(x)
  } else if (is.matrix(x)) {
    N <- nrow(x)
  } else {
    stop('x must be a vector or a matrix')
  }
  if(!isPosscal(width)) {
    stop('width must be a positive scalar')
  }
  if(!isPosscal(rc)) {
    stop('rc must be a positive scalar')
  }
  if(!isPosscal(fs)) {
    stop('fs must be a positive scalar')
  }

  if (width * fs > N / 2) {
    idx <- c(1, N)
    w <- rep(1L, N)
  } else {
    idx <- round((N + width * fs * c(-1, 1)) / 2)
    w <- sigmoid_train((1:N), idx, rc * fs)$y
  }

  rmsx_singleChan <- function (x, w, N, idx) {
    fx    <- stats::fft(as.vector(x)^2)
    fw    <- stats::fft(as.vector(w)^2)
    out  <- ifft(fx * fw) / (N-1)
    out[out < .Machine$double.eps * max(out)] <- 0
    out <- pracma::circshift (sqrt(out), round(mean(idx)))
  }

  if (is.vector(x)) {
    rmsx <- rmsx_singleChan(x, w, N, idx)
  } else {
    rmsx <- apply(x, 2, rmsx_singleChan, w, N, idx)
  }

  list(rmsx = rmsx, w = w)
}
