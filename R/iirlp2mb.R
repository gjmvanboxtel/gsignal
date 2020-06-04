# iirlp2mb.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Octave signal package:
# Copyright (C) 2011 Alan J. Greenberger <alanjg@ptd.net>
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
# 20200604  Geert van Boxtel              First version for v0.1.0
#---------------------------------------------------------------------------------------------------------------------------------

#' IIR lowpass filter to IIR multiband
#' 
#' Transform an IIR lowpass filter prototype to an IIR multiband filter.
#' 
#' The utility of a prototype filter comes from the property that all other
#' filters can be derived from it by applying a scaling factor to the components
#' of the prototype. The filter design need thus only be carried out once in
#' full, with other filters being obtained by simply applying a scaling factor.
#' Especially useful is the ability to transform from one bandform to another.
#' In this case, the transform is more than a simple scale factor. Bandform here
#' is meant to indicate the category of passband that the filter possesses. The
#' usual bandforms are lowpass, highpass, bandpass and bandstop, but others are
#' possible. In particular, it is possible for a filter to have multiple
#' passbands. In fact, in some treatments, the bandstop filter is considered to
#' be a type of multiple passband filter having two passbands. Most commonly,
#' the prototype filter is expressed as a lowpass filter, but other techniques
#' are possible[1].
#' 
#' Filters with multiple passbands may be obtained by applying the general
#' transformation described in [2].
#' 
#' Because \code{iirlp2mb} is generic, it can be extended to accept other inputs.
#' 
#' @param b numerator polynomial of prototype low pass filter
#' @param a denominator polynomial of prototype low pass filter
#' @param Wo (normalized angular frequency)/pi to be transformed
#' @param Wt vector of (norm. angular frequency)/pi transform targets
#' @param type one of "pass" or "stop". Specifies to filter to produce: bandpass
#'   (default) or bandstop.
#' @param ... additional arguments (not used)
#' 
#' @return list of class \code{'\link{Arma}'} numerator and denominator
#'   polynomials of the resulting filter
#' 
#' @examples
#' ## Design a prototype real IIR lowpass elliptic filter with a gain of about
#' ## –3 dB at 0.5pi rad/sample.
#' el <- ellip(3, 0.1, 30, 0.409)
#' ## Create a real multiband filter with two passbands.
#' mb1 <- iirlp2mb(el, 0.5, c(.2, .4, .6, .8), 'pass')
#' ## Create a real multiband filter with two stopbands.
#' mb2 <- iirlp2mb(el, 0.5, c(.2, .4, .6, .8), 'stop')
#' ## Compare the magnitude responses of the filters.
#' hfl <- freqz(el)
#' hf1 <- freqz(mb1)
#' hf2 <- freqz(mb2)
#' plot(hfl$w, 20 * log10(abs(hfl$h)), type = "l",
#'     xlab = "Normalized frequency (* pi rad/sample)",
#'     ylab = "Magnitude (dB)")
#' lines(hf1$w, 20 * log10(abs(hf1$h)), col="red")
#' lines(hf2$w, 20 * log10(abs(hf2$h)), col="blue")
#' legend('bottomleft', legend = c('Prototype', 'Two passbands', 'Two Stopbands'),
#'        col=c("black", "red", "blue"), lty = 1)
#'
#' @references [1] \url{https://en.wikipedia.org/wiki/Prototype_filter}\cr
#' [2] \url{https://en.wikipedia.org/wiki/Prototype_filter#Lowpass_to_multi-band}
#' 
#' @author Original Octave code by Alan J. Greenberger \email{alanjg@@ptd.net}.
#'   Port to R by Geert van Boxtel \email{G.J.M.vanBoxtel@@gmail.com}.
#'   
#' @rdname iirlp2mb
#' @export

iirlp2mb <- function(b, ...) UseMethod("iirlp2mb")

#' @rdname iirlp2mb
#' @export

iirlp2mb.Arma <- function(b, Wo, Wt, type, ...) {
  iirlp2mb (b$b, b$a, Wo, Wt, type, ...)
}

#' @rdname iirlp2mb
#' @export

iirlp2mb.Zpg <- function (b, Wo, Wt, type, ...) {
  ba <- as.Arma(b)
  iirlp2mb (ba$b, ba$a, Wo, Wt, type, ...)
}

#' @rdname iirlp2mb
#' @export

iirlp2mb.Sos <- function (b, Wo, Wt, type, ...) {
  ba <- as.Arma(b)
  iirlp2mb (ba$b, ba$a, Wo, Wt, type, ...)
}

#' @rdname iirlp2mb
#' @export

iirlp2mb.default <- function (b, a, Wo, Wt, type = c("pass", "stop"), ...) {

  # input validation
  type <- match.arg(type)
  if(type == "pass") {
    pass_stop <- -1
  } else if (type == "stop") {
    pass_stop <- 1
  }
  if (!isPosscal(Wo) || Wo > 1) {
    stop('Frequency value Wo of prototype filter must be a scalar between 0 and 1')
  }
  if (any(Wt < 0) || any(Wt > 1)) {
    stop('Frequency values Wt of target filter must be between 0 and 1')
  }
  oWt <- unique(sort(Wt))

  ##                                                             B(z)
  ## Inputs B,A specify the low pass IIR prototype filter G(z) = ---- .
  ##                                                             A(z)
  ## This module transforms G(z) into a multiband filter using the iterative
  ## algorithm from:
  ## [FFM] G. Feyh, J. Franchitti, and C. Mullis, "All-Pass Filter
  ## Interpolation and Frequency Transformation Problem", Proceedings 20th
  ## Asilomar Conference on Signals, Systems and Computers, Nov. 1986, pp.
  ## 164-168, IEEE.
  ## [FFM] moves the prototype filter position at normalized angular frequency
  ## .5*pi to the places specified in the Wt vector times pi.  In this module,
  ## a generalization allows the position to be moved on the prototype filter
  ## to be specified as Wo*pi instead of being fixed at .5*pi.  This is
  ## implemented using two successive allpass transformations.
  ##                                         KK(z)
  ## In the first stage, find allpass J(z) = ----  such that
  ##                                         K(z)
  ##    jWo*pi     -j.5*pi
  ## J(e      ) = e                    (low pass to low pass transformation)
  ##
  ##                                          PP(z)
  ## In the second stage, find allpass H(z) = ----  such that
  ##                                          P(z)
  ##    jWt(k)*pi     -j(2k - 1)*.5*pi
  ## H(e         ) = e                 (low pass to multiband transformation)
  ##
  ##                                          ^
  ## The variable PP used here corresponds to P in [FFM].
  ## len = length(P(z)) == length(PP(z)), the number of polynomial coefficients
  ##
  ##        len      1-i           len       1-i
  ## P(z) = SUM P(i)z   ;  PP(z) = SUM PP(i)z   ; PP(i) == P(len + 1 - i)
  ##        i=1                    i=1              (allpass condition)
  ## Note: (len - 1) == n in [FFM] eq. 3
  ##
  ## The first stage computes the denominator of an allpass for translating
  ## from a prototype with position .5 to one with a position of Wo. It has the
  ## form:
  ##          -1
  ## K(2)  - z
  ## -----------
  ##          -1
  ## 1 - K(2)z
  ##
  ## From the low pass to low pass transformation in Table 7.1 p. 529 of A.
  ## Oppenheim and R. Schafer, Discrete-Time Signal Processing 3rd edition,
  ## Prentice Hall 2010, one can see that the denominator of an allpass for
  ## going in the opposite direction can be obtained by a sign reversal of the
  ## second coefficient, K(2), of the vector K (the index 2 not to be confused
  ## with a value of z, which is implicit).
  
  ## The first stage allpass denominator computation
  K <- apd(pi * Wo)
  
  ## The second stage allpass computation
  phi <- pi * Wt                                # vector of normalized angular frequencies between 0 and pi
  P <- apd(phi)                                 # calculate denominator of allpass for this target vector
  PP <- rev(P)                                  # numerator of allpass has reversed coefficients of P
  
  ## The total allpass filter from the two consecutive stages can be written as
  ##          PP
  ## K(2) -  ---
  ##          P         P
  ## -----------   *   ---
  ##          PP        P
  ## 1 - K(2)---
  ##          P
  AllpassDen <- P - (K[2] * PP)
  AllpassDen <- AllpassDen / AllpassDen[1]      # normalize
  AllpassNum <- pass_stop * rev(AllpassDen)
  ba <- transform(b, a, AllpassNum, AllpassDen, pass_stop)
  ba
}


######################################################################################
# Helper functions for iirlp2mb, not exported from the namespace

# all pass denominator
apd <- function (phi) {
  ## Given phi, a vector of normalized angular frequency transformation targets,
  ## return the denominator of an allpass H(z)
  Pkm1 <- 1                 # P0 initial condition from [FFM] eq. 22
  for (k in 1:length(phi)) {
    P <- pk(Pkm1, k, phi[k])
    Pkm1 <- P
  }
  P
}

# kth iteration of P(z)
pk <- function (Pkm1, k, phik){

  ## Given Pkminus1, k, and phi(k) in radians , return Pk
  ##
  ## From [FFM] eq. 19 :                     k
  ## Pk =     (z+1  )sin(phi(k)/2)Pkm1 - (-1) (z-1  )cos(phi(k)/2)PPkm1
  ## Factoring out z
  ##              -1                         k    -1
  ##    =   z((1+z  )sin(phi(k)/2)Pkm1 - (-1) (1-z  )cos(phi(k)/2)PPkm1)
  ## PPk can also have z factored out.  In H=PP/P, z in PPk will cancel z in Pk,
  ## so just leave out.  Use
  ##              -1                         k    -1
  ## PK =     (1+z  )sin(phi(k)/2)Pkm1 - (-1) (1-z  )cos(phi(k)/2)PPkm1
  ## (expand)                                k
  ##    =            sin(phi(k)/2)Pkm1 - (-1)        cos(phi(k)/2)PPkm1
  ##
  ##              -1                         k   -1
  ##           + z   sin(phi(k)/2)Pkm1 + (-1)    z   cos(phi(k)/2)PPkm1

  Pk <- rep(0L, k + 1)  # there are k+1 coefficients in Pk
  sin_k <- sin(phik / 2)
  cos_k <- cos(phik / 2)
  for (i in 1:k) {
    Pk[i] <- Pk[i] +  sin_k * Pkm1[i] - ((-1)^k * cos_k * Pkm1[k + 1 - i]);
    ##
    ##                    -1
    ## Multiplication by z   just shifts by one coefficient
    Pk[i + 1] <- Pk[i + 1] + sin_k * Pkm1[i] + ((-1)^k * cos_k * Pkm1[k + 1 - i]);
  }
  ## now normalize to Pk(1) = 1 (again will cancel with same factor in PPk)
  Pk <- Pk / Pk[1]
  Pk
}


# Regenerate ith power of P from stored PPower
ppower <- function (Ppower, i, powcols) {

  if(i == 0) {
    p  <- 1
  } else {
    p  <- NULL
    for (j in 1:powcols) {
      if (is.na(Ppower[i, j])) break
      p =  cbind(p, Ppower[i,j])
    }
  }
  p
}

# add polynomials of possibly different length
polysum <- function (p1, p2) {

  n1 <- length(p1)
  n2 <- length(p2)
  if (n1 > n2) {
    ## pad p2
    p2 <- c(p2, rep(0L, n1 - n2))
  } else if (n2 > n1) {
    ## pad p1
    p1 = c(p1, rep(0L, n2 - n1))
  }
  poly <- p1 + p2
  poly
}

transform <- function (B, A, PP, P, pass_stop) {

  ## Given G(Z) = B(Z)/A(Z) and allpass H(z) = PP(z)/P(z), compute G(H(z))
  ## For Pass = 'pass', transformed filter is:
  ##                          2                   nb-1
  ## B1 + B2(PP/P) + B3(PP/P)^  + ... + Bnb(PP/P)^
  ## -------------------------------------------------
  ##                          2                   na-1
  ## A1 + A2(PP/P) + A3(PP/P)^  + ... + Ana(PP/P)^
  ## For Pass = 'stop', use powers of (-PP/P)
  ##
  na <- length(A)  # the number of coefficients in A
  nb <- length(B)  # the number of coefficients in B
  ## common low pass iir filters have na == nb but in general might not
  n  = max(na, nb) # the greater of the number of coefficients
  ##                              n-1
  ## Multiply top and bottom by P^   yields:
  ##
  ##      n-1             n-2          2    n-3                 nb-1    n-nb
  ## B1(P^   ) + B2(PP)(P^   ) + B3(PP^ )(P^   ) + ... + Bnb(PP^    )(P^    )
  ## ---------------------------------------------------------------------
  ##      n-1             n-2          2    n-3                 na-1    n-na
  ## A1(P^   ) + A2(PP)(P^   ) + A3(PP^ )(P^   ) + ... + Ana(PP^    )(P^    )

  ## Compute and store powers of P as a matrix of coefficients because we will
  ## need to use them in descending power order
  np <- length(P)
  powcols <- np + (np - 1) * (n - 2) # number of coefficients in P^(n-1)
  ## initialize to "Not Available" with n-1 rows for powers 1 to (n-1) and
  ## the number of columns needed to hold coefficients for P^(n-1)
  Ppower <- matrix(NA, nrow = n - 1, ncol = powcols) # to hold coefficients of powers of P
  Ptemp <- P                      # start with P to the 1st power
  for (i in 1:(n - 1)) {          # i is the power
    for (j in 1:length(Ptemp)) {  # j is the coefficient index for this power
      Ppower[i, j]  <- Ptemp[j]
    }
    Ptemp <- conv(Ptemp, P)       # increase power of P by one
  }

  ## Compute numerator and denominator of transformed filter
  Num <- Den <- NULL
  for (i in 1:n) {
    ##              n-i
    ## Regenerate P^    (p_pownmi)
    if((n - i) == 0) {
      p_pownmi <- 1
    } else {
      p_pownmi <- ppower(Ppower, n - i, powcols)
    }
    ##               i-1
    ## Regenerate PP^   (pp_powim1)
    if(i == 1) {
      pp_powim1 <- 1
    } else{
      pp_powim1 <- rev(ppower(Ppower, i-1, powcols))
    } 
    if(i <= nb) {
      Bterm <- (pass_stop^(i - 1)) * B[i] * conv(pp_powim1, p_pownmi)
      Num <- polysum(Num, Bterm)
    }
    if(i <= na) {
      Aterm <- (pass_stop^(i - 1)) * A[i] * conv(pp_powim1, p_pownmi)
      Den <- polysum(Den, Aterm)
    }
  }
  ## Scale both numerator and denominator to have Den(1) = 1
  Den <- Den / Den[1]
  Num <- Num / Den[1]

  Arma(Num, Den)
}

