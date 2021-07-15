# ultrwin.R
# Copyright (C) 2019 Geert van Boxtel <gjmvanboxtel@gmail.com>
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
# 20191215 Geert van Boxtel          First version for v0.1.0
# 20210715 Geert van Boxtel          Rewritten from scratch v0.3-3
#---------------------------------------------------------------------------------------------------------------------------------

#' Ultraspherical window
#' 
#' Return the coefficients of an ultraspherical window
#' 
#' @param n Window length, specified as a positive integer.
#' @param mu parameter that controls the side-lobe roll-off ratio. Default: 3.
#' @param xmu parameters that provides a trade-off between the ripple ratio and
#'   a width characteristic. Default: 1
#' 
#' @return ultraspherical window, returned as a vector.
#' 
#' @examples
#' 
#' w <- ultrwin(101, 3, 1)
#' plot (w, type = "l", xlab = "Samples", ylab =" Amplitude")
#' freqz(w)
#' 
#' w2 <- ultrwin(101, 2, 1)
#' f2 <- freqz(w2)
#' w3 <- ultrwin(101, 3, 1)
#' f3 <- freqz(w3)
#' w4 <- ultrwin(101, 4, 1)
#' f4 <- freqz(w4)
#' op <- par(mfrow = c(2, 1))
#' plot(w2, type = "l", col = "black", xlab = "", ylab = "")
#' lines(w3, col = "red")
#' lines(w4, col = "blue")
#' legend("topright", legend = 2:4, col = c("black", "red", "blue"), lty = 1)
#' plot (f2$w, 20 * log10(abs(f2$h)), type = "l", col = "black",
#'       xlab = "", ylab = "", ylim = c(-100, 50))
#' lines(f3$w, 20 * log10(abs(f3$h)), col = "red")
#' lines(f4$w, 20 * log10(abs(f4$h)), col = "blue")
#' legend("topright", legend = 2:4, col = c("black", "red", "blue"), lty = 1)
#' par(op)
#' title(main = "Effect of increasing the values of mu (xmu = 1)")
#'
#' w1 <- ultrwin(101, 2, 1)
#' f1 <- freqz(w1)
#' w2 <- ultrwin(101, 2, 1.001)
#' f2 <- freqz(w2)
#' w3 <- ultrwin(101, 2, 1.002)
#' f3 <- freqz(w3)
#' op <- par(mfrow = c(2, 1))
#' plot(w1, type = "l", col = "black", xlab = "", ylab = "")
#' lines(w2, col = "red")
#' lines(w3, col = "blue")
#' legend("topright", legend = 2:4, col = c("black", "red", "blue"), lty = 1)
#' plot (f1$w, 20 * log10(abs(f1$h)), type = "l", col = "black",
#'       xlab = "", ylab = "", ylim = c(-100, 50))
#' lines(f2$w, 20 * log10(abs(f2$h)), col = "red")
#' lines(f3$w, 20 * log10(abs(f3$h)), col = "blue")
#' legend("topright", legend = c(1, 1.001, 1.002),
#'        col = c("black", "red", "blue"), lty = 1)
#' par(op)
#' title(main = "Effect of increasing the values of xmu (mu = 2)")
#'
#' @note The Dolph-Chebyshev and Saramaki windows are special cases of the
#'   Ultraspherical window, with mu set to 0 and 1, respectively.
#' 
#' @author Geert van Boxtel \email{G.J.M.vanBoxtel@@gmail.com}.
#' 
#' @references [1] Bergen, S.W.A., and Antoniou, A. Design of Ultraspherical
#'   Window Functions with Prescribed Spectral Characteristics. EURASIP Journal
#'   on Applied Signal Processing 2004:13, 2053–2065.
#
#' @export

ultrwin <- function (n, mu = 3, xmu = 1) {
  
  if (!isPosscal(n) || !isWhole(n) || n <= 0) stop ("n must be a positive integer")
  if (!isScalar(mu) || !is.double(mu)) stop ("mu must be a real scalar")
  if (!isScalar(xmu) || !is.double(xmu)) stop ("xmu must be a real scalar")
  
  if (n == 1) {
    w <- 1
  } else {
    m <- (n - 1) / 2
    if (n%%2) {
      a <- 0
    } else {
      a <- 0.5
    }
    b <- 1 - xmu^(-2)
    w <- rep(0, n)
    for (k in seq(0, m)) {
      nn = k - m
      #cat("k+1=",k+1,"nn=",nn,"nn-k=",n-k,"\n")
      w[k + 1] <- usc(nn, m, b, mu, xmu)
      w[n - k] <- w[k + 1]
    }
  }
  w <- w / usc(a, m, b, mu, xmu)
  w
}

# function to calculate the ultraspherical coefficients
# not exported to the namespace
usc <- function(n, m, b, mu, xmu) {
  sum <- 0
  for (mm  in seq(0, (m - abs(n)))) {
    sum <- sum + 
      choose(mu + m - abs(n) - 1, m - abs(n) - mm) *
      choose(m + abs(n), mm) * b^(mm)
  }
  if (mu == 0) {
    cf <- xmu^(2 * m) / (m + abs(n)) * 
      choose(mu + m + abs(n) - 1, m + abs(n) - 1) * sum
  } else {
    cf <- mu * xmu^(2 * m) / (m + abs(n)) * 
      choose(mu + m + abs(n) - 1, m + abs(n) - 1) * sum
  }
  cf
}

