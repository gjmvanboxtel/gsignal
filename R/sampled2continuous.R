# sampled2continuous.R
# Copyright (C) 2020 Geert van Bxtel <G.J.M.vanBoxtel@gmail.com>
# Original Octave code:
# Copyright (C) 2009 Muthiah Annamalai <muthiah.annamalai@uta.edu>
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
# 20201124  GvB       setup for gsignal v0.1.0
#------------------------------------------------------------------------------

#' Signal reconstruction
#'
#' Analog signal reconstruction from discrete samples.
#'
#' Given a discrete signal x[n] sampled with a frequency of \code{fs} Hz, this
#' function reconstruct the original analog signal x(t) at time points \code{t}.
#' The function can be used, for instance, to calculate sampling rate effects on
#' aliasing.
#'
#' @param xn the sampled input signal, specified as a vector
#' @param fs sampling frequency in Hz used in collecting \code{x}, specified as
#'   a positive scalar value. Default: 1
#' @param t time points at which data is to be reconstructed, specified as a
#'   vector relative to \code{x[0]} (not real time).
#'
#' @return Reconstructed signal x(t), returned as a vector.
#'
#' @examples
#' # 'analog' signal: 3 Hz cosine
#' t <- seq(0, 1, length.out = 100)
#' xt <- cos(3 * 2 * pi * t)
#' plot(t, xt, type = "l", xlab = "", ylab = "", ylim = c(-1, 1.2))
#'
#' # 'sample' it at 4 Hz to simulate aliasing
#' fs <- 4
#' n <- ceiling(length(t) / fs)
#' xn <- xt[seq(ceiling(n / 2), length(t), n)]
#' s4 <- sampled2continuous(xn, fs, t)
#' lines(t, s4, col = "red")
#'
#' # 'sample' it > 6 Hz to avoid aliasing
#' fs <- 7
#' n <- ceiling(length(t) / fs)
#' xn <- xt[seq(ceiling(n / 2), length(t), n)]
#' s7 <- sampled2continuous(xn, fs, t)
#' lines(t, s7, col = "green")
#' legend("topright", legend = c("original", "aliased", "non-aliased"),
#'   lty = 1, col = c("black", "red", "green"))
#'
#'
#' @author Muthiah Annamalai, \email{muthiah.annamalai@@uta.edu}.
#' Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#'
#' @export

sampled2continuous <- function(xn, fs, t) {

  if (!is.vector(xn)) {
    stop("xn must be a vector")
  }
  if (!isPosscal(fs) || fs <= 0) {
    stop("fs must be a positive scalar")
  }
  if (!is.vector(t)) {
    stop("t must be a vector")
  }

  T <- 1 / fs
  N <- length(xn)
  mg <- pracma::meshgrid(T * seq(0, N - 1), t)
  S <- sinc((mg$Y - mg$X) / T)
  xt <- S %*% xn
  xt
}
