# residue.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Original Octave function:
# Copyright (C) 1994-2017 John W. Eaton
# Copyright (C) 2007 Ben Abbott
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
# 20200612  GvB       setup for gsignal v0.1.0
#------------------------------------------------------------------------------

#' Partial fraction expansion
#'
#' Finds the residues, poles, and direct term of a Partial Fraction Expansion of
#' the ratio of two polynomials.
#'
#' The call \code{res <- residue(b, a)} computes the partial fraction expansion
#' for the quotient of the polynomials, \code{b} and \code{a}.
#'
#' The call \code{res <- rresidue(r, p, k)} performs the inverse operation and
#' computes the reconstituted quotient of polynomials, b(s) / a(s), from the
#' partial fraction expansion; represented by the residues, poles, and a direct
#' polynomial specified by \code{r}, \code{p} and \code{k}, and the pole
#' multiplicity \code{e}.
#'
#' @param b coefficients of numerator polynomial
#' @param a coefficients of denominator polynomial
#' @param r residues of partial fraction expansion
#' @param p poles of partial fraction expansion
#' @param k direct term
#' @param tol tolerance. Default: 0.001
#'
#' @return For \code{residue}, a list containing \code{r}, \code{p} and
#'   \code{k}. For \code{rresidue}, a list containing \code{b} and \code{a}.
#'
#' @examples
#' b <- c(-4, 8)
#' a <- c(1, 6, 8)
#' rpk <- residue(b, a)
#' ba <- rresidue(rpk$r, rpk$p, rpk$k)
#'
#' @author Tony Richardson, \email{arichard@@stark.cc.oh.us},\cr
#'  Ben Abbott, \email{bpabbott@@mac.com},\cr
#'  adapted by John W. Eaton.\cr
#'   Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}
#'
#' @rdname residue
#' @export

residue <- function(b, a, tol = 0.001) {

  if (!is.vector(b) || !is.vector(a)) {
    stop("b and a must be vectors")
  }
  tol <- abs(tol[1])

  ## Make sure both polynomials are in reduced form.
  a <- polyreduce(a)
  b <- polyreduce(b)

  b <- b / a[1]
  a <- a / a[1]

  la <- length(a)
  lb <- length(b)

  ## Handle special cases here.
  if (la == 0 || lb == 0) {
    return(list(r = NULL, p = NULL, k = NULL))
  } else if (la == 1) {
    k <- b / a
    return(list(r = NULL, p = NULL, k = k))
  }

  ## Find the poles.
  p <- pracma::roots(a)
  lp <- length(p)

  ## Sort poles so that multiplicity loop will work.
  mn <- mpoles(p, tol = tol, reorder = TRUE, index.return = TRUE)
  p <- p[mn$n]

  ## For each group of pole multiplicity, set the value of each
  ## pole to the average of the group.  This reduces the error in
  ## the resulting poles.
  p_group <- cumsum(mn$m == 1)
  for (ng in seq_len(length(p_group))) {
    m <- which(p_group == ng)
    p[m] <- mean(p[m])
  }

  ## Find the direct term if there is one.
  if (lb >= la) {
    ## Also return the reduced numerator.
    qr <- pracma::deconv(b, a)
    k <- qr$q
    b <- qr$r
    lb <- length(b)
  } else {
    k <- NULL
  }

  ## Determine if the poles are (effectively) zero.
  small <- max(abs(p))
  small <- max(small, 1) * .Machine$double.eps * 1e4 * (1 + length(p))^2
  p[abs(p) < small] <- 0

  ## Determine if the poles are (effectively) real, or imaginary.
  index <- (abs(Im(p)) < small)
  if (any(index)) {
    p[index] <- Re(p[index])
  }
  index <- (abs(Re(p)) < small)
  if (any(index)) {
    p[index] <- 1i * Im(p[index])
  }

  ## The remainder determines the residues.  The case of one pole
  ## is trivial.
  if (lp == 1) {
    r <- pracma::polyval(b, p)
  } else {

    ## Determine the order of the denominator and remaining numerator.
    ## With the direct term removed the potential order of the numerator
    ## is one less than the order of the denominator.
    aorder <- length(a) - 1
    border <- aorder - 1

    ## Construct a system of equations relating the individual
    ## contributions from each residue to the complete numerator.
    A <- matrix(0L, nrow = border + 1, ncol = border + 1)
    B <- prepad(b, border + 1, 0)
    for (ip in seq_along(p)) {
      ri <- rep(0L, length(p))
      ri[ip] <- 1
      A[, ip] <- prepad(rresidue(ri, p, NULL, tol)$b, border + 1, 0)
    }

    ## Solve for the residues.
    r <- as.vector(pracma::mldivide(A, B, pinv = FALSE))
  }

  r <- zapIm(r)
  p <- zapIm(p)
  list(r = r, p = p, k = k)
}

#' @rdname residue
#' @export

rresidue <- function(r, p, k, tol = 0.001) {

  if (!is.vector(r) || !is.vector(p)) {
    stop("r and p must be vectors")
  }
  tol <- abs(tol[1])

  mn <- mpoles(p, tol, reorder = FALSE, index.return = TRUE)
  indx <- mn$n
  p <- p[indx]
  r <- r[indx]

  indx <- seq_along(p)
  for (n in indx) {
    pn <- c(1, -p[n])
    if (n == 1) {
      pden <- pn
    } else {
      pden <- fftconv(pden, pn)
    }
  }

  ## D is the order of the denominator
  ## K is the order of the direct polynomial
  ## N is the order of the resulting numerator
  ## pnum(1:(N+1)) is the numerator's polynomial
  ## pden(1:(D+1)) is the denominator's polynomial
  ## pm is the multible pole for the nth residue
  ## pn is the numerator contribution for the nth residue

  D <- length(pden) - 1
  K <- length(k) - 1
  N <- K + D
  pnum <- rep(0L, N + 1)
  for (n in indx[abs(r) > 0]) {
    p1 <- c(1, -p[n])
    pn <- 1
    if (n > 1) {
      for (j in 1:(n - 1)) {
        pn <- fftconv(pn, c(1, -p[j]))
      }
    }
    if (n + 1 <= length(p)) {
      for (j in (n + 1):length(p)) {
        pn <- fftconv(pn, c(1, -p[j]))
      }
    }
    if (mn$m[n] > 1) {
      for (j in 1:(mn$m[n] - 1)) {
        pn <- pracma::deconv(pn, p1)$q
      }
    }
    pn <- r[n] * pn
    pnum <- pnum + prepad(pn, N + 1, 0, 2)
  }

  ## Add the direct term.
  if (length(k)) {
    pnum <- pnum + fftconv(pden, k)
  }

  pnum <- polyreduce(pnum)
  pden <- polyreduce(pden)

  list(b = zapIm(pnum), a = zapIm(pden))
}
