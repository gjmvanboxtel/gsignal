# ultrwin.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Octave signal package:
# Copyright (C) 2013 Rob Sykes <robs@users.sourceforge.net>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; see the file COPYING. If not, see
# <https://www.gnu.org/licenses/>.
#
# 20201116 Geert van Boxtel          First version for v0.1.0
#---------------------------------------------------------------------------------------------------------------------------------

#' Ultraspherical window
#' 
#' Return the coefficients of an ultraspherical window
#' 
#' Windows can be characterized as fixed or adjustable. Fixed windows have only
#' one independent parameter, namely, the window length which controls the
#' main-lobe width. Adjustable windows have one or more additional parameters
#' that can control other window characteristics, such as the relative side-lobe
#' amplitude (e.g., \code{kaiser}. \code{chebwin}). The ultraspherical window
#' has an additional parameter that can be used to set the rate at which
#' side-lobes decrease (or increase) in amplitude as a function of frequency.
#' This can be useful, for instance, if a very sharp attenuation immediately
#' next to the transition band is desired.
#' 
#' Ref [1] details methods to set the parameters that can be used for designing
#' specific spectral characteristics of the window.
#' 
#' @param n Window length, specified as a positive integer (controls main lobe
#'   width)
#' @param mu parameter that controls the side-lobe roll-off ratio of the Fourier
#'   transform of the window.
#' @param pvalue parameter value used for controlling the transform’s main-lobe
#'   width/side-lobe-ratio.
#' @param ptype type of parameter passed in \code{pvalue}. Can be one of:
#'   \describe{
#'     \item{xmu (default)}{controls the ripple ratio by setting the
#'     (un-normalized) window’s Fourier transform according to its canonical
#'     definition:
#'       \preformatted{
#'        (MU)
#' W(k) = C   [ xmu cos(pi k/M) ],  k = 0, 1, ..., M-1,
#'        M-1
#'      }
#'      where C is the Ultraspherical (a.k.a. Gegenbauer) polynomial, which can
#'      be defined using the recurrence relationship:
#'        \preformatted{
#'      (l)    1                  (l)                    (l)
#'     C (x) = - [ 2x(m + l - 1) C   (x) - (m + 2l - 2) C   (x) ]
#'      m      m                  m-1                    m-2
#'                                 (l)        (l)
#'      for m an integer > 1, and C (x) = 1, C (x) = 2lx.
#'                                 0          1
#'      }
#'      Note that when not giving \code{xmu}, stability issues may occur with
#'      \code{mu <= -1.5}.
#'    }
#'     \item{beta}{sets the main lobe width to \code{beta} times that of a
#'     rectangular window.}
#'     \item{att}{sets the ripple ratio at the first }
#'     \item{latt}{sets the ripple ratio at the last side-lobe}
#'   }
#' 
#' @return ultraspherical window, returned as a vector.
#' 
#' @examples
#' ## Window with sharp attenuation next to transition band
#' w <- ultrwin(120, -1, 40, "l")
#' fz <- freqz(w)
#' op <- par(mfrow = c(2, 1))
#' plot (seq(0, length(w) - 1), w, type = "l", xlab = "", ylab ="")
#' plot (fz$w / pi, 20 * log10(abs(fz$h) / abs(fz$h[1])), type = "l", 
#'       xlab = "", ylab = "")
#' par(op)
#'       
#' ## Varying beta with fixed mu
#' x <- y <- NULL
#' for (beta in 2:5) {
#'   w <- ultrwin(80, -.5, beta, "beta"); fz  <- freqz(w)
#'   x <- cbind(x, fz$w / pi)
#'   y <- cbind(y, 20 * log10(abs(fz$h) / abs(fz$h[1])))
#' }
#' matplot(x, y, type = "l", lty = 1, xlab = "", ylab = "", ylim = c(-150,0),
#'   main = expression(paste("Varying ", beta, " with ", mu, " = 0.5")))
#' legend("topright", legend = 2:5, title = "Beta", lty = 1, col = 1:4)
#' 
#' ## Varying n with fixed mu and beta
#' x <- y <- NULL
#' for (n in 2:10) {
#'   w <- ultrwin(n * 20, 1, 3, "beta"); fz  <- freqz(w, 1, 2^11)
#'   x <- cbind(x, fz$w / pi)
#'   y <- cbind(y, 20 * log10(abs(fz$h) / abs(fz$h[1])))
#' }
#' matplot(x, y, type = "l", lty = 1, xlab = "", ylab = "", 
#'   xlim = c(0, 0.25), ylim = c(-100, 0), col = 1:9,
#'   main = expression(paste("Varying n with ", mu, " = 1 and ", beta, " = 3")))
#' legend("topright", legend = 20*(2:10), title = "n", lty = 1, col = 1:9)
#' 
#' ## Varying mu with fixed m and att
#' x <- y <- NULL
#' for (j in 0:4) {
#'   w <- ultrwin(80, j * .6 - 1.2, 50, "att"); fz  <- freqz(w)
#'   x <- cbind(x, fz$w / pi)
#'   y <- cbind(y, 20 * log10(abs(fz$h) / abs(fz$h[1])))
#' }
#' matplot(x, y, type = "l", lty = 1, xlab = "", ylab = "", 
#'   xlim = c(0, 1), ylim = c(-100, 0), col = 1:5,
#'   main = expression(paste("Varying ", mu, " with n = 80 and att = 50")))
#' legend("topright", legend = (0:4) * .6 - 1.2, title = expression(mu),
#'   lty = 1, col = 1:5)
#' 
#' ## Varying mu with fixed m and latt
#' x <- y <- NULL
#' for (j in rev(0:4)) {
#'   w <- ultrwin(80, j * .75 - 1.5, 50, "latt"); fz  <- freqz(w)
#'   x <- cbind(x, fz$w / pi)
#'   y <- cbind(y, 20 * log10(abs(fz$h) / abs(fz$h[1])))
#' }
#' matplot(x, y, type = "l", lty = 1, xlab = "", ylab = "", 
#'   xlim = c(0, 1), ylim = c(-100, 0), col = 1:5,
#'   main = expression(paste("Varying ", mu, " with n = 80 and latt = 50")))
#' legend("topright", legend = (0:4) * .6 - 1.2, title = expression(mu),
#'   lty = 1, col = 1:5)
#'   
#' ## Compare ultraspherical, Dolph-Chebyshev and Kaiser windows
#' x <- y <- NULL
#' for (i in 1:3) {
#'   w <- switch(i, ultrwin(153, 0.5, 2.6, "beta"), ultrwin(165, 0, 2.73, "beta"), kaiser(159, 7.91))
#'   fz <- freqz(w, fs = 1.2 * pi)
#'   x <- cbind(x, fz$w)
#'   y <- cbind(y, 20 * log10(abs(fz$h) / abs(fz$h[1])))
#' }
#' matplot(x, y, type = "l", lty = 1, xlab = "", ylab = "", col = 1:3,
#'   xlim = c(0, 1), ylim = c(-130, 0))
#' legend("topright", lty = 1, col = 1:3,
#'   legend = c("Ultraspherical", "Dolph-Chebyshev", "Kaiser"))
#'   
#' @note The Dolph-Chebyshev and Saramaki windows are special cases of the
#'   Ultraspherical window, with mu set to 0 and 1 respectively.
#' 
#' @author Rob Sykes, \email{robs@@users.sourceforge.net}; port to R by Geert
#'   van Boxtel \email{G.J.M.vanBoxtel@@gmail.com}.
#'
#' @references [1] Bergen, S.W.A., Antoniou, A. (2004). Design of Ultraspherical Window Functions
#' with Prescribed Spectral Characteristics. EURASIP J. Appl. Sign. Proc. 13, 2053-2065\cr
#' [2] \url{https://en.wikipedia.org/wiki/Window_function#Ultraspherical_window}
#
#' @export

ultrwin <- function (n, mu, pvalue, ptype = "beta") {
  
  if (!isPosscal(n) || !isWhole(n) || n <= 0) stop ("n must be a positive integer")
  mu <- as.double(mu)
  if (!isScalar(mu)) stop ("mu must be a real scalar")
  pvalue <- as.double(pvalue)
  if (!isScalar(pvalue)) stop ("paramater value must be a real scalar")
  types <- c("xmu", "beta", "att", "latt")
  ptype = match.arg(ptype, types)
  
  
  w <- .Call("_gsignal_ultrwin", PACKAGE = "gsignal", n, mu, pvalue, which(types == ptype) - 1, 0L)
  if (is.null(w)) {
    stop("Parameter(s) out of range")
  }
  w
}