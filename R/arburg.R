# arburg.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Original Octave function:
# Copyright (C) 2006 Peter V. Lanspeary <pvl@mecheng.adelaide.edu.au>
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
# 20201104  GvB       setup for gsignal v0.1.0
# 20210507  GvB       bugfix in inner product to compute v
#------------------------------------------------------------------------------

#' Autoregressive model coefficients - Burg's method
#'
#' Calculate the coefficients of an autoregressive model using the whitening
#' lattice-filter method of Burg (1968)[1].
#'
#' The inverse of the autoregressive model is a moving-average filter which
#' reduces \code{x} to white noise. The power spectrum of the AR model is an
#' estimate of the maximum entropy power spectrum of the data. The function
#' \code{ar_psd} calculates the power spectrum of the AR model.
#'
#' For data input \code{x(n)} and white noise \code{e(n)}, the autoregressive
#' model is
#' \if{latex}{
#' \deqn{x(n) = \sqrt{v} \cdot e(n) + \sum_{k=1}^{p+1} a(k) \cdot x(n-k)}
#' }
#' \if{html}{\preformatted{
#'                           p+1
#'     x(n) = sqrt(v).e(n) + SUM a(k).x(n-k)
#'                           k=1
#'  }}
#'
#' \code{arburg} does not remove the mean from the data. You should remove the
#' mean from the data if you want a power spectrum. A non-zero mean can produce
#' large errors in a power-spectrum estimate.  See \code{\link{detrend}}
#'
#' @note AIC, AICc, KIC and AKICc are based on information theory. They  attempt
#'   to balance the complexity (or length) of the model against how well the
#'   model fits the data.  AIC and KIC are biased estimates of the asymmetric
#'   and the symmetric Kullback-Leibler divergence, respectively. AICc and AKICc
#'   attempt to correct the bias. See reference [2].
#'
#' @param x input data, specified as a numeric or complex vector or matrix. In
#'   case of a vector it represents a single signal; in case of a matrix each
#'   column is a signal.
#' @param p model order; number of poles in the AR model or limit to the number
#'   of poles if a valid criterion is provided. Must be < length(x) - 2.
#' @param criterion model-selection criterion. Limits the number of poles so
#'   that spurious poles are not added when the whitened data has no more
#'   information in it. Recognized values are:
#'   \describe{
#'     \item{AKICc}{approximate corrected Kullback information criterion
#'     (recommended)}
#'     \item{KIC}{Kullback information criterion}
#'     \item{AICc}{corrected Akaike information criterion}
#'     \item{AIC}{Akaike information criterion}
#'     \item{FPE}{final prediction error}
#'   }
#'   The default is to NOT use a model-selection criterion (NULL)
#'
#' @return A \code{list} containing the following elements:
#'   \describe{
#'     \item{a}{vector or matrix containing \code{(p+1)} autoregression
#'     coefficients. If \code{x} is a matrix, then each row of a corresponds to
#'     a column of \code{x}. \code{a} has \code{p + 1} columns.}
#'     \item{e}{white noise input variance, returned as a vector. If \code{x} is
#'     a matrix, then each element of e corresponds to a column of \code{x}.}
#'     \item{k}{Reflection coefficients defining the lattice-filter embodiment
#'     of the model returned as vector or a matrix. If \code{x} is a matrix,
#'     then each column of \code{k} corresponds to a column of \code{x}.
#'     \code{k} has \code{p} rows.}
#'   }
#'
#' @examples
#' A <- Arma(1, c(1, -2.7607, 3.8106, -2.6535, 0.9238))
#' y <- filter(A, 0.2 * rnorm(1024))
#' coefs <- arburg(y, 4)
#'
#' @author Peter V. Lanspeary, \email{pvl@@mecheng.adelaide.edu.au}.\cr
#'  Conversion to R by Geert van Boxtel, \email{gjmvanboxtel@@gmail.com}.
#'
#' @references [1] Burg, J.P. (1968) A new analysis technique for time series
#'   data, NATO advanced study Institute on Signal Processing with Emphasis on
#'   Underwater Acoustics, Enschede, Netherlands, Aug. 12-23, 1968.\cr
#'   [2] Seghouane, A. and Bekara, M. (2004). A small sample model selection
#'   criterion based on Kullback’s symmetric divergence. IEEE Trans. Sign.
#'   Proc., 52(12), pp 3314-3323,
#'
#' @seealso \code{\link{ar_psd}}
#'
#' @export

arburg <- function(x, p, criterion = NULL) {

  # check parameters
  if (!(is.vector(x) || is.matrix(x)) || !is.numeric(x)) {
    stop("x must be a numeric or vector or matrix")
  }

  if (is.vector(x)) {
    vec <- TRUE
    x <- as.matrix(x, ncol = 1)
  } else {
    vec <- FALSE
  }
  nr <- nrow(x)
  nc <- ncol(x)

  if (!isScalar(p) || !isWhole(p) || !is.numeric(p) || p <= 0.5) {
    stop("p must be a positive integer")
  }
  if (p >= nr - 2) {
    stop(paste0("p must be less than the length of x (", nr, ") - 2"))
  }

  if (!is.null(criterion)) {
    criterion <- match.arg(criterion, c("AKICc", "KIC", "AICc", "AIC", "FPE"))
    #  Set the model-selection-criterion flags.
    #  is_akicc, isa_kic and is_corrected are short-circuit flags
    is_akicc <- criterion == "AKICc"                 # AKICc
    isa_kic  <- is_akicc || criterion == "KIC"       # KIC or AKICc
    is_corrected <- is_akicc || criterion == "AICc"  # AKICc or AICc
    use_inf_crit <- is_corrected || isa_kic || criterion == "AIC"
    use_fpe <- criterion == "FPE"
  } else {
    use_inf_crit <- FALSE
    use_fpe <- FALSE
  }
  # end of parameter checking

  # loop over columns
  aggr_a <- aggr_v <- aggr_k <- NULL
  for (icol in seq_len(nc)) {

    # f(n) = forward prediction error
    # b(n) = backward prediction error
    # Storage of f(n) and b(n) is a little tricky. Because f(n) is always
    # combined with b(n-1), f(1) and b(N) are never used, and therefore are
    # not stored.  Not storing unused data makes the calculation of the
    # reflection coefficient look much cleaner :)
    # N.B. {initial v} = {error for zero-order model} =
    #      {zero-lag autocorrelation} =  E(x*conj(x)) = x*x'/N
    #      E = expectation operator
    f <- x[2:nr, icol]
    b <- x[1:(nr - 1), icol]
    v <- Re(x[, icol] %*% x[, icol]) / nr
    # new_crit/old_crit is the mode-selection criterion
    new_crit <- abs(v)
    old_crit <- 2 * new_crit

    for (ip in seq_len(p)) {

      # new reflection coeff = -2* E(f.conj(b)) / ( E(f^2)+E(b(^2) )
      last_k <- as.vector(-2 * (b %*% f) / (f %*% f + b %*% b))
      ##  Levinson-Durbin recursion for residual
      new_v <- v * (1.0 - Re(last_k * Conj(last_k)))
      if (ip > 1) {

        # Apply the model-selection criterion and break out of loop if it
        # increases (rather than decreases).
        # Do it before we update the old model "a" and "v".

        # * Information Criterion (AKICc, KIC, AICc, AIC)
        if (use_inf_crit) {
          old_crit <- new_crit
          new_crit <- log(new_v) + as.integer(is_akicc) * ip / nr / (nr - ip) +
            (2 + as.integer(isa_kic) - as.integer(is_akicc) * (ip + 2) / nr) *
            (ip + 1) / (nr - as.integer(is_corrected) * (ip + 2))
          if (new_crit > old_crit) {
            break
          }
        # (FPE) Final prediction error
        } else if (use_fpe) {
          old_crit <- new_crit
          new_crit <- new_v * (nr + ip + 1) / (nr - ip - 1)
          if (new_crit > old_crit) {
            break
          }
        }
        ## Update model "a" and "v".
        ## Use Levinson-Durbin recursion formula (for complex data).
        a <- c(prev_a + last_k * Conj(prev_a[seq(ip - 1, 1, -1)]), last_k)
      } else {  # if(ip==1 )
        a <- last_k
        k <- NULL
      }
      k <- c(k, last_k)
      v <- new_v
      if (ip < p) {
        prev_a <- a
        #  calculate new prediction errors (by recursion):
        #  f(p,n) = f(p-1,n)   + k * b(p-1,n-1)        n=2,3,...n
        #  b(p,n) = b(p-1,n-1) + conj(k) * f(p-1,n)    n=2,3,...n
        #  remember f(p,1) is not stored, so don't calculate it; make f(p,2)
        #  the first element in f.  b(p,n) isn't calculated either.
        nn <- nr - ip
        new_f <- f[2:nn] + last_k * b[2:nn]
        b <- b[1:(nn - 1)] + Conj(last_k) * f[1:(nn - 1)]
        f <- new_f
      }
    }  # loop over p
    aggr_a <- rbind(aggr_a, c(1, a))
    aggr_v <- c(aggr_v, v)
    aggr_k <- rbind(aggr_k, k)
  }    # loop over signals

  if (vec) {
    rv <- list(a = as.vector(aggr_a), e = aggr_v, k = as.vector(aggr_k))
  } else {
    rv <- list(a = aggr_a, e = aggr_v, k = t(aggr_k))
  }
  rv
}
