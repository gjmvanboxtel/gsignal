# filternorm.R
# Copyright (C) 2026 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Original Octave function:
# Copyright (C) 2023 Leonardo Araujo <leolca@gmail.com>
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
# 20260909  GvB       setup for gsignal v0.4.0
#------------------------------------------------------------------------------

#' 2-norm or infinity-norm of digital filter
#'
#' Compute the 2-norm or the infinity-norm of a stable digital filter. This is
#' useful for quantifying the gain or energy amplification of a digital filter,
#' especially when designing and implementing FIR or IIR filters. For example,
#' if a filter has a large infinity-norm, it may produce a large peak output for
#' a bounded input. If it has a large 2-norm, it may strongly amplify signal
#' energy or noise. Thus, \code{filternorm} is mainly a design and
#' implementation diagnostic—not usually a replacement for examining the full
#' frequency response with tools such as \code{\link{freqz}}.
#'
#' @param filt For the default case, the moving-average coefficients of an ARMA
#'   filter (normally called \code{b}), specified as a numeric or complex
#'   vector. Generically, \code{filt} specifies an arbitrary filter operation.
#' @param a the autoregressive (recursive) coefficients of an ARMA filter,
#'   specified as a numeric or complex vector. If \code{a[1]} is not equal to 1,
#'   then filter normalizes the filter coefficients by \code{a[1]}. Therefore,
#'   \code{a[1]} must be nonzero.
#' @param pnorm Type of filter norm to compute for the digital filter, 
#'   specified as one of these:
#'   \describe{
#'    \item{2}{compute the 2-norm of the filter (default).}
#'    \item{Inf}{computes the infinity-norm of the filter.}
#'   }
#' @param tol Tolerance allowed for filter norm estimation, specified as a
#'   positive scalar. Default: 1e-8, currently ignored.
#' @param ... additional arguments (ignored).
#
#' @return Filter norm, returned as a real positive scalar.
#'
#' @examples
#' 
#' ## 1. Comparing two low-pass filters
#' fs <- 1000
#' fc <- 100
#' 
#' # 4th-order Butterworth IIR filter
#' iir <- butter(4, fc / (fs / 2))
#' 
#' # 80th-order FIR filter
#' fir <- fir1(80, fc / (fs / 2))
#' 
#' n1 <- filternorm(iir)
#' n2 = filternorm(fir)
#' cat("IIR norm:", n1)
#' cat("FIR norm:", n2)
#' 
#' # This lets you compare the filters’ overall amplification or energy
#' # behavior. A filter with a larger norm may amplify numerical errors or
#' #internal signals more strongly, even if both filters have similar passband
#' #responses.
#' 
#' 
#' ## 2. Normalizing a filter
#' 
#' # A filter may have excessive gain because of its coefficient values or
#' # design requirements:
#' 
#' n <- filternorm(ba)
#' b_normalized <- ba$b / n
#' a_normalized <- ba$a
#' ba_normalized <- Arma(b_normalized, a_normalized)
#' 
#' # The normalized filter has a reduced overall gain. This can be useful before
#' # exporting coefficients to DSP or embedded hardware.
#' # However, normalization changes the filter’s amplitude response, so you may
#' # need to restore a desired passband gain afterward:
#' desired_gain <- 2
#' b_normalized <- desired_gain * ba$b / n
#' ba_normalized <- Arma(b_normalized, a_normalized)
#' 
#' 
#' ## 3. Comparing filters with respect to noise amplification
#' 
#' n_butter <- filternorm(ba)
#' 
#' # A simple moving-average FIR filter
#' ma <- Ma(rep(1, 8) / 8)
#' 
#' n_ma <- filternorm(ma)
#' cat("Butterworth norm:", n_butter)
#' cat("Moving-average norm:", n_ma)
#' 
#' # If one filter has a substantially larger norm, it may amplify round-off
#' # noise or quantization effects more strongly. This is particularly relevant
#' # when implementing filters with short word lengths.
#' 
#' 
#' ## 4. Scaling sections in a cascaded filter
#' 
#' # Suppose a filter is implemented as two cascaded sections:
#' 
#' baA <- butter(2, 0.2)
#' baB <- butter(2, 0.5)
#' 
#' nA <- filternorm(baA)
#' nB <- filternorm(baB)
#' 
#' cat("Section A norm:", nA)
#' cat("Section B norm:", nB)
#' 
#' # If section A has much greater gain than section B, the internal signal
#' # between the sections may become unnecessarily large. You can redistribute
#' # gain:
#' bA_scaled <- baA$b / nA
#' baA_scaled <- Arma(bA_scaled, baA$a)
#' bB_scaled <- baB$b * nA
#' baB_scaled <- Arma(bB_scaled, baB$a)
#' 
#' # The cascade’s overall gain is approximately preserved, while the
#' # intermediate signal produced by section A is reduced.
#' 
#' ### In practice, filternorm should be used alongside other tools such as
#' ### freqz, rather than as the only measure of filter quality.
#' 
#'
#' @author Leonardo Araujo \email{leolca@@gmail.com>}.\cr Conversion to R by
#'   Geert van Boxtel, \email{gjmvanboxtel@@gmail.com}.
#'
#' @rdname filternorm
#' @export
#' 
filternorm <- function(filt, ...) UseMethod("filternorm")

#' @rdname filternorm
#' @method filternorm default
#' @export

filternorm.default <- function(filt, a, pnorm = 2, tol = 1e-8, ...) {
  
  if (!is.vector(filt) || ! is.vector(a)) {
    stop("b and a must be numeric vectors")
  }

  if (pnorm != 2 && pnorm != Inf) {
    stop("pnorm must be 2 or Inf")
  }
  
  if (!isPosscal(tol)) {
    stop("tol must be a positive scalar")
  }
  
  if (pnorm == 2) {
    ## Parseval's theorem states that the L2-norm of a filter with frequency
    ## response H(e^{j\omega}) is the square-root of the sum of the squares
    ## of its filter impulse response (the energy of the impulse response).
    h <- impz(filt, a)$x
    L <- sqrt(sum(abs(h)^2))
  }
  else if (pnorm == Inf) {
    ## the norm in L-infinity is simply the maximum of the frequency response:
    ## ||H||_{\infty} = \max_{0 \leq \omega \leq \pi} { |H(e^{j\omega})|}
    h <- freqz(filt, a, 1024)$h
    L <- max(abs(h))
  }
  L
}

#' @rdname filternorm
#' @method filternorm Arma
#' @export
filternorm.Arma <- function(filt, pnorm = 2, tol = 1e-8, ...) # IIR
  filternorm(filt$b, filt$a, pnorm, tol, ...)

#' @rdname filternorm
#' @method filternorm Ma
#' @export
filternorm.Ma <- function(filt, pnorm = 2, tol = 1e-8, ...) # FIR
  filternorm(unclass(filt), 1, pnorm, tol, ...)

#' @rdname filternorm
#' @method filternorm Sos
#' @export
filternorm.Sos <- function(filt, pnorm = 2, tol = 1e-8, ...) { # Second-order sections
  filternorm(as.Arma(filt), pnorm, tol, ...)
}

#' @rdname filternorm
#' @method filternorm Zpg
#' @export
filternorm.Zpg <- function(filt, pnorm = 2, tol = 1e-8, ...) # zero-pole-gain form
  filternorm(as.Arma(filt), pnorm, tol, ...)
