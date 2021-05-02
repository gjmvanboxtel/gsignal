# sftrans.R
# Copyright (C) 2019 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Octave signal package:
# Copyright (C) 1999-2001 Paul Kienzle <pkienzle@users.sf.net>
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
# 20200501 Geert van Boxtel          First version for v0.1.0
#------------------------------------------------------------------------------

#' Transform filter band edges
#'
#' Transform band edges of a generic lowpass filter to a filter with different
#' band edges and to other filter types (high pass, band pass, or band stop).
#'
#' Given a low pass filter represented by poles and zeros in the splane, you can
#' convert it to a low pass, high pass, band pass or band stop by transforming
#' each of the poles and zeros individually. The following summarizes the
#' transformations:
#' \if{latex}{
#'   \tabular{lll}{
#'      \strong{Transform}               \tab \strong{Zero at x}                   \tab \strong{Pole at x}  \cr
#'      -------------------------        \tab -------------------------            \tab ------------------------- \cr
#'      \strong{Low-Pass}                \tab  zero: \eqn{Fc x/C}                  \tab  pole: \eqn{Fc x/C} \cr
#'      \eqn{S \rightarrow C S/Fc}       \tab  gain: \eqn{C/Fc}                    \tab  gain: \eqn{Fc/C}   \cr
#'      -------------------------        \tab -------------------------            \tab ------------------------- \cr
#'      \strong{High Pass}               \tab  zero: \eqn{Fc C/x}                  \tab  pole: \eqn{Fc C/x} \cr
#'      \eqn{S \rightarrow C Fc/S}       \tab  pole: \eqn{0}                       \tab  zero: \eqn{0} \cr
#'                                       \tab  gain: \eqn{-x}                      \tab  gain: \eqn{-1/x} \cr
#'      -------------------------        \tab -------------------------            \tab ------------------------- \cr
#'      \strong{Band Pass}               \tab  zero: \eqn{b +- \sqrt{(b^2-FhFl)}}  \tab  pole: \eqn{b \pm \sqrt{(b^2-FhFl)}} \cr
#'                                       \tab  pole: \eqn{0}                       \tab  zero: \eqn{0} \cr
#'      S -> \eqn{C \frac{S^2+FhFl}{S(Fh-Fl)}}
#'                                       \tab gain: \eqn{C/(Fh-Fl)}                \tab  gain: \eqn{(Fh-Fl)/C} \cr
#'                                       \tab  \eqn{b=x/C (Fh-Fl)/2}               \tab  \eqn{b=x/C (Fh-Fl)/2} \cr
#'      -------------------------        \tab -------------------------            \tab ------------------------- \cr
#'      \strong{Band Stop}               \tab  zero: \eqn{b \pm \sqrt{(b^2-FhFl)}} \tab  pole: \eqn{b +- \sqrt{(b^2-FhFl)}} \cr
#'                                       \tab  pole: \eqn{\pm \sqrt{(-FhFl)}}      \tab  zero: \eqn{\pm \sqrt{(-FhFl)}}      \cr
#'      S -> \eqn{C \frac{S(Fh-Fl)}{S^2+FhFl}}
#'                                       \tab  gain: \eqn{-x}                      \tab  gain: \eqn{-1/x} \cr
#'                                       \tab  \eqn{b=C/x (Fh-Fl)/2}               \tab  \eqn{b=C/x (Fh-Fl)/2} \cr
#'      -------------------------        \tab -------------------------            \tab ------------------------- \cr
#'      \strong{Bilinear}                \tab zero: \eqn{(2+xT)/(2-xT)}            \tab  pole: \eqn{(2+xT)/(2-xT)} \cr
#'                                       \tab  pole: \eqn{-1}                      \tab  zero: \eqn{-1} \cr
#'      \eqn{S \rightarrow \frac{2 z-1}{T z+1}}
#'                                       \tab  gain: \eqn{(2-xT)/T}                \tab  gain: \eqn{(2-xT)/T} \cr
#'      -------------------------        \tab -------------------------            \tab ------------------------- \cr
#' }}
#' \if{html}{\preformatted{
#'
#'   Transform         Zero at x                  Pole at x
#'   ----------------  -------------------------  --------------------------
#'   Low-Pass          zero: Fc x/C               pole: Fc x/C
#'   S -> C S/Fc       gain: C/Fc                 gain: Fc/C
#'   ----------------  -------------------------  --------------------------
#'   High Pass         zero: Fc C/x               pole: Fc C/x
#'   S -> C Fc/S       pole: 0                    zero: 0
#'                     gain: -x                   gain: -1/x
#'   ----------------  -------------------------  --------------------------
#'   Band Pass         zero: b +- sqrt(b^2-FhFl)  pole: b +- sqrt(b^2-FhFl)
#'          S^2+FhFl   pole: 0                    zero: 0
#'   S -> C --------   gain: C/(Fh-Fl)            gain: (Fh-Fl)/C
#'          S(Fh-Fl)   b=x/C (Fh-Fl)/2            b=x/C (Fh-Fl)/2
#'   ----------------  -------------------------  --------------------------
#'   Band Stop         zero: b +- sqrt(b^2-FhFl)  pole: b +- sqrt(b^2-FhFl)
#'          S(Fh-Fl)   pole: +-sqrt(-FhFl)        zero: +-sqrt(-FhFl)
#'   S -> C --------   gain: -x                   gain: -1/x
#'          S^2+FhFl   b=C/x (Fh-Fl)/2            b=C/x (Fh-Fl)/2
#'   ----------------  -------------------------  --------------------------
#'   Bilinear          zero: (2+xT)/(2-xT)        pole: (2+xT)/(2-xT)
#'        2 z-1        pole: -1                   zero: -1
#'   S -> -----        gain: (2-xT)/T             gain: (2-xT)/T
#'        T z+1
#'   ----------------  -------------------------  --------------------------
#' }}
#'
#' where C is the cutoff frequency of the initial lowpass filter, F_c is the
#' edge of the target low/high pass filter and [F_l,F_h] are the edges of the
#' target band pass/stop filter. With abundant tedious algebra, you can derive
#' the above formulae yourself by substituting the transform for S into
#' \eqn{H(S)=S-x} for a zero at x or \eqn{H(S)=1/(S-x)} for a pole at x, and
#' converting the result into the form:
#'
#' \deqn{g prod(S-Xi) / prod(S-Xj)}
#'
#' Please note that a pole and a zero at the same place exactly cancel. This is
#' significant for High Pass, Band Pass and Band Stop filters which create
#' numerous extra poles and zeros, most of which cancel. Those which do not
#' cancel have a ‘fill-in’ effect, extending the shorter of the sets to have the
#' same number of as the longer of the sets of poles and zeros (or at least
#' split the difference in the case of the band pass filter). There may be other
#' opportunistic cancellations, but it does not check for them.
#'
#' Also note that any pole on the unit circle or beyond will result in an
#' unstable filter. Because of cancellation, this will only happen if the number
#' of poles is smaller than the number of zeros and the filter is high pass or
#' band pass. The analytic design methods all yield more poles than zeros, so
#' this will not be a problem.
#'
#' @param Sz In the generic case, a model to be transformed. In the default case,
#'   a vector containing the zeros in a pole-zero-gain model.
#' @param Sp a vector containing the poles in a pole-zero-gain model.
#' @param Sg a vector containing the gain in a pole-zero-gain model.
#' @param w critical frequencies of the target filter specified in radians.
#'   \code{w} must be a scalar for low-pass and high-pass filters, and \code{w}
#'   must be a two-element vector c(low, high) specifying the lower and upper
#'   bands in radians.
#' @param stop FALSE for a low-pass or band-pass filter, TRUE for a high-pass or
#'   band-stop filter.
#' @param ...	arguments passed to the generic function.
#'
#' @return For the default case or for sftrans.Zpg, an object of class "Zpg",
#'   containing the list elements:
#' \describe{
#'   \item{z}{complex vector of the zeros of the transformed model}
#'   \item{p}{complex vector of the poles of the transformed model}
#'   \item{g}{gain of the transformed model}
#' }
#' For sftrans.Arma, an object of class "Arma", containing the list elements:
#' \describe{
#'   \item{b}{moving average (MA) polynomial coefficients}
#'   \item{a}{autoregressive (AR) polynomial coefficients}
#' }
#'
#' @examples
#' ## 6th order Bessel bandpass
#' zpg <- besselap(6)
#' bp <- sftrans(zpg, c(2, 3), stop = TRUE)
#' freqs(bp, seq(0, 4, length.out = 128))
#' bp <- sftrans(zpg, c(0.1,0.3), stop = FALSE)
#' freqs(bp, seq(0, 4, length.out = 128))
#'
#' @references Proakis & Manolakis (1992). \emph{Digital Signal Processing}. New
#'   York: Macmillan Publishing Company.
#'
#' @author Paul Kienzle, \email{pkienzle@@users.sf.net}.\cr
#'   Conversion to R by Tom Short,\cr
#'   adapted by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#'
#' @rdname sftrans
#' @export

sftrans <- function(Sz, ...) UseMethod("sftrans")

#' @rdname sftrans
#' @export

sftrans.Zpg <- function(Sz, w, stop = FALSE, ...)
  sftrans.default(Sz$z, Sz$p, Sz$g, w, stop)

#' @rdname sftrans
#' @export

sftrans.Arma <- function(Sz, w, stop = FALSE, ...)
  as.Arma(sftrans(as.Zpg(Sz), w, stop))

#' @rdname sftrans
#' @export

sftrans.default <- function(Sz, Sp, Sg, w, stop = FALSE, ...)  {

  if (is.null(Sz)) Sz <- 0   #GvB 20200428
  C <- 1
  p <- length(Sp)
  z <- length(Sz)
  if (z > p || p == 0) {
    stop("must have at least as many poles as zeros in s-plane")
  }

  if (length(w) == 2) {
    Fl <- w[1]
    Fh <- w[2]
    if (stop) {
      ## ----------------  -------------------------  ------------------------
      ## Band Stop         zero: b +- sqrt(b^2-FhFl)   pole: b +- sqrt(b^2-FhFl)
      ##        S(Fh-Fl)   pole: +-sqrt(-FhFl)         zero: +-sqrt(-FhFl)
      ## S -> C --------   gain: -x                   gain: -1/x
      ##        S^2+FhFl   b=C/x (Fh-Fl)/2            b=C/x (Fh-Fl)/2
      ## ----------------  -------------------------  ------------------------
      Sg <- Sg * Re(prod(-Sz) / prod(-Sp))
      b <- (C * (Fh - Fl) / 2) / Sp
      Sp <- c(b + sqrt(0i + b^2 - Fh * Fl), b - sqrt(0i + b^2 - Fh * Fl))
      extend <- c(sqrt(0i + -Fh * Fl), -sqrt(0i + -Fh * Fl))
      if (is.null(Sz) || length(Sz) == 0) {
        Sz <- extend[1 + (1:(2 * p)) %% 2]
      } else {
        b <- (C * (Fh - Fl) / 2) / Sz
        Sz <- c(b + sqrt(0i + b^2 - Fh * Fl), b - sqrt(0i + b^2 - Fh * Fl))
        if (p > z) {
          Sz <- c(Sz, extend[1 + ((1:2) * (p - z)) %% 2])
        }
      }
    } else {
      ## ----------------  -------------------------  ------------------------
      ## Band Pass         zero: b +- sqrt(b^2-FhFl)   pole: b +- sqrt(b^2-FhFl)
      ##        S^2+FhFl   pole: 0                    zero: 0
      ## S -> C --------   gain: C/(Fh-Fl)            gain: (Fh-Fl)/C
      ##        S(Fh-Fl)   b=x/C (Fh-Fl)/2            b=x/C (Fh-Fl)/2
      ## ----------------  -------------------------  ------------------------
      Sg <- Sg * (C / (Fh - Fl)) ^ (z - p)
      b <- Sp * (Fh - Fl) / (2 * C)
      Sp <- c(b + sqrt(0i + b^2 - Fh * Fl), b - sqrt(0i + b^2 - Fh * Fl))
      if (is.null(Sz) || length(Sz) == 0) {
        Sz <- numeric(p)
      } else {
        b <- Sz * (Fh - Fl) / (2 * C)
        Sz <- c(b + sqrt(0i + b^2 - Fh * Fl), b - sqrt(0i + b^2 - Fh * Fl))
        if (p > z) {
          Sz <- c(Sz, numeric(p - z))
        }
      }
    }
  } else {
    Fc <- w
    if (stop) {
      ## ----------------  -------------------------  ------------------------
      ## High Pass         zero: Fc C/x               pole: Fc C/x
      ## S -> C Fc/S       pole: 0                    zero: 0
      ##                   gain: -x                   gain: -1/x
      ## ----------------  -------------------------  ------------------------
      Sg <- Sg * Re(prod(-Sz) / prod(-Sp))
      Sp <- C * Fc / Sp
      if (is.null(Sz) || length(Sz) == 0) {
        Sz <- numeric(p)
      } else {
        Sz <- C * Fc / Sz
        if (p > z) {
          Sz <- c(Sz, numeric(p - z))
        }
      }
    } else {
      ## ----------------  -------------------------  ------------------------
      ## Low Pass          zero: Fc x/C               pole: Fc x/C
      ## S -> C S/Fc       gain: C/Fc                 gain: Fc/C
      ## ----------------  -------------------------  ------------------------
      Sg <- Sg * (C / Fc) ^ (z - p)
      Sp <- Fc * Sp / C
      Sz <- Fc * Sz / C
    }
  }
  Zpg(z = Sz, p = Sp, g = Sg)
}
