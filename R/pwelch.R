# pwelch.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Original Octave code:
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
# 20201001  GvB       setup for gsignal v0.1.0
# 20201112  GvB       bug in assigning colnames when ns > 1
# 20210302  GvB       bug when x is matrix: xx[1:seg_len] -> xx[1:seg_len, ]
# 20210408  GvB       added range parameter ('half' o 'whole')
# 20211021  GvB       corrected bug in col ptr into coh, phase, trans matrices
# 20220405  GvB       remove padding to nearest power of 2 (Github Disc #6)
#                     bug in returning vector when ncol is 1 (Github Issue #5)
# 20220511  GvB       use inherits() instead of direct comparison of class name
# 20220512  GvB       plot method for class 'pwelch'
#------------------------------------------------------------------------------

#' Welch’s power spectral density estimate
#'
#' Compute power spectral density (PSD) using Welch's method.
#'
#' The Welch method [1] reduces the variance of the periodogram estimate to the
#' PSD by splitting the signal into (usually) overlapping segments and windowing
#' each segment, for instance by a Hamming window. The periodogram is then
#' computed for each segment, and the squared magnitude is computed, which is
#' then averaged for all segments. See also [2].
#'
#' The spectral density is the mean of the modified periodograms, scaled so that
#' area under the spectrum is the same as the mean square of the data.  This
#' equivalence is supposed to be exact, but in practice there is a mismatch of
#' up to 0.5% when comparing area under a periodogram with the mean square of
#' the data.
#'
#' In case of multivariate signals, Cross-spectral density, phase, and coherence
#' are also returned. The input data can be demeaned or detrended, overall or
#' for each segment separately.
#'
#' @param x input data, specified as a numeric vector or matrix. In case of a
#'   vector it represents a single signal; in case of a matrix each column is a
#'   signal.
#' @param window If \code{window} is a vector, each segment has the same length
#'    as \code{window} and is multiplied by \code{window} before (optional)
#'    zero-padding and calculation of its periodogram. If \code{window} is a
#'    scalar, each segment has a length of \code{window} and a Hamming window is
#'    used. Default: \code{nextpow2(sqrt(length(x)))} (the square root of the
#'    length of \code{x} rounded up to the next power of two). The window length
#'    must be larger than 3.
#' @param overlap segment overlap, specified as a numeric value expressed as a
#'   multiple of window or segment length. 0 <= overlap < 1. Default: 0.5.
#' @param nfft Length of FFT, specified as an integer scalar. The default is the
#'   length of the \code{window} vector or has the same value as the scalar
#'   \code{window} argument.  If \code{nfft} is larger than the segment length,
#'   (seg_len), the data segment is padded \code{nfft - seg_len} zeros. The
#'   default is no padding. Nfft values smaller than the length of the data
#'   segment (or window) are ignored. Note that the use of padding to increase
#'   the frequency resolution of the spectral estimate is controversial.
#' @param fs sampling frequency (Hertz), specified as a positive scalar.
#'   Default: 1.
#' @param detrend character string specifying detrending option; one of:
#'   \describe{
#'     \item{\code{long-mean}}{remove the mean from the data before
#'     splitting into segments (default)}
#'     \item{\code{short-mean}}{remove the mean value of each segment}
#'     \item{\code{long-linear}}{remove linear trend from the data before
#'     splitting into segments}
#'     \item{\code{short-linear}}{remove linear trend from each segment}
#'     \item{\code{none}}{no detrending}
#'  }
#' @param range character string. one of:
#' \describe{
#'   \item{\code{"half"} or \code{"onesided"}}{frequency range of the spectrum
#'   is from zero up to but not including \code{fs / 2}. Power from negative
#'   frequencies is added to the positive side of the spectrum.}
#'   \item{\code{"whole"} or \code{"twosided"}}{frequency range of the spectrum
#'   is \code{-fs / 2} to \code{fs / 2}, with negative frequencies stored in
#'   "wrap around order" after the positive frequencies; e.g. frequencies for a
#'   10-point \code{"twosided"} spectrum are 0 0.1 0.2 0.3 0.4 0.5 -0.4 -0.3
#'   -0.2. -0.1.}
#'   \item{\code{"shift"} or \code{"centerdc"}}{same as \code{"whole"} but with
#'   the first half of the spectrum swapped with second half to put the
#'   zero-frequency value in the middle.}
#' }
#'   Default: If \code{x} are real, the default range is \code{"half"},
#'   otherwise the default range is \code{"whole"}.
#' @param plot.type character string specifying which plot to produce; one of
#'   \code{"spectrum"}, \code{"cross-spectrum"}, \code{"phase"},
#'   \code{"coherence"}, \code{"transfer"}
#' @param yscale character string specifying scaling of Y-axis; one of
#'   \code{"linear"}, \code{"log"}, \code{"dB"}
#' @param xlab,ylab,main labels passed to plotting function. Default: NULL
#' @param ... additional arguments passed to functions
#'
#' @return An object of class \code{"pwelch"}, which is a list containing the
#'   following elements:
#'   \describe{
#'     \item{\code{freq}}{vector of frequencies at which the spectral variables
#'     are estimated. If \code{x} is numeric, power from negative frequencies is
#'     added to the positive side of the spectrum, but not at zero or Nyquist
#'     (fs/2) frequencies. This keeps power equal in time and spectral domains.
#'     If \code{x} is complex, then the whole frequency range is returned.}
#'     \item{\code{spec}}{Vector (for univariate series) or matrix (for
#'     multivariate series) of estimates of the spectral density at frequencies
#'     corresponding to freq.}
#'     \item{\code{cross}}{NULL for univariate series. For multivariateseries, a
#'     matrix containing the cross-spectral density estimates between different
#'     series. Column \eqn{i + (j - 1) * (j - 2)/2 } of contains the
#'     cross-spectral estimates between columns \eqn{i} and \eqn{j} of \eqn{x},
#'     where \eqn{i < j}.}
#'     \item{\code{phase}}{NULL for univariate series. For multivariate series,
#'     a matrix containing the cross-spectrum phase between different series.
#'     The format is the same as \code{cross}.}
#'     \item{\code{coh}}{NULL for univariate series. For multivariate series, a
#'     matrix containing the squared coherence between different series. The
#'     format is the same as \code{cross}.}
#'     \item{\code{trans}}{NULL for univariate series. For multivariate series,
#'     a matrix containing estimates of the transfer function between different
#'     series. The format is the same as \code{cross}.}
#'     \item{\code{x_len}}{The length of the input series.}
#'     \item{\code{seg_len}}{The length of each segment making up the averages.}
#'     \item{\code{psd_len}}{The number of frequencies. See \code{freq}}
#'     \item{\code{nseries}}{The number of series}
#'     \item{\code{series}}{The name of the series}
#'     \item{\code{snames}}{For multivariate input, the names of the individual
#'     series}
#'     \item{\code{window}}{The window used to compute the modified periodogram}
#'     \item{\code{fs}}{The sampling frequency}
#'     \item{\code{detrend}}{Character string specifying detrending option}
#'   }
#'
#' @examples
#' fs <- 256
#' secs <- 10
#' freq <- 30
#' ampl <- 1
#' t <- seq(0, secs, length.out = fs * secs)
#'
#' x <- ampl * cos(freq * 2 * pi * t) + runif(length(t))
#' Pxx <- pwelch(x, fs = fs)              # no plot
#' pwelch(x, fs = fs)                     # plot
#'
#' # 90 degrees phase shift with with respect to x
#' y <- ampl * sin(freq * 2 * pi * t) + runif(length(t))
#' Pxy <- pwelch(cbind(x, y), fs = fs)
#' plot(Pxy, yscale = "dB")
#' plot(Pxy, plot.type = "phase")
#' # note the phase shift around 30 Hz is pi/2
#' plot(Pxy, plot.type = "coherence")
#'
#' # Transfer function estimate example
#' fs <- 1000                 # Sampling frequency
#' t <- (0:fs) / fs           # One second worth of samples
#' A <- c(1, 2)               # Sinusoid amplitudes
#' f <- c(150, 140)           # Sinusoid frequencies
#' xn <- A[1] * sin(2 * pi * f[1] * t) +
#'       A[2] * sin(2 * pi * f[2] * t) +  0.1 * runif(length(t))
#' h <- Ma(rep(1L, 10) / 10)      # Moving average filter
#' yn <- filter(h, xn)
#' atfm <- freqz(h, fs = fs)
#' etfm <- pwelch(cbind(xn, yn), fs = fs)
#' op <- par(mfrow = c(2, 1))
#' xl <- "Frequency (Hz)"; yl <- "Magnitude"
#' plot(atfm$w, abs(atfm$h), type = "l", main = "Actual", xlab = xl, ylab = yl)
#' plot(etfm$freq, abs(etfm$trans), type = "l", main = "Estimated",
#'      xlab = xl, ylab = yl)
#' par(op)
#'
#' @note Unlike the 'Octave' function 'pwelch', the current implementation
#'   does not compute confidence intervals because they can be inaccurate in
#'   case of overlapping segments.
#'
#' @author Peter V. Lanspeary \email{pvl@@mecheng.adelaide.edu.au}.\cr
#' Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#'
#' @references [1] Welch, P.D. (1967). The use of Fast Fourier Transform for
#'   the estimation of power spectra: A method based on time averaging over
#'   short, modified periodograms. IEEE Transactions on Audio and
#'   Electroacoustics, AU-15 (2): 70–73.\cr
#' @references [2] \url{https://en.wikipedia.org/wiki/Welch\%27s_method}
#'
#' @rdname pwelch
#' @export

pwelch <- function(x, window = nextpow2(sqrt(NROW(x))), overlap = 0.5,
                   nfft = if (isScalar(window)) window else length(window),
                   fs = 1,
                   detrend = c("long-mean", "short-mean",
                               "long-linear", "short-linear", "none"),
                   range = if (is.numeric(x)) "half" else "whole") {

  # check parameters
  if (!(is.vector(x) || is.matrix(x)) &&
      !(is.numeric(x) || is.complex(x))) {
    stop("x must be a numeric or complex vector or matrix")
  }

  series <- deparse(substitute(x))
  if (is.vector(x)) {
    vec <- TRUE
    x <- as.matrix(x, ncol = 1)
    snames <- ""
  } else {
    vec <- FALSE
    snames <- colnames(x)
  }
  x_len <- nrow(x)
  ns <- ncol(x)
  if (is.null(snames)) {
    snames <- colnames(x) <- as.character(seq_len(ns))
  }

  if (!is.vector(window) || !is.numeric(window)) {
    stop("window must be a numeric vector or scalar")
  } else {
    if (isPosscal(window)) {
      if (window <= 3) {
        stop("window must be a scalar > 3 or a vector with length > 3")
      } else {
        window <- hamming(window)
      }
    } else if (length(window) <= 3) {
      stop("window must be a scalar > 3 or a vector with length > 3")
    }
  }

  if (!isScalar(overlap) || !(overlap >= 0 && overlap < 1)) {
    stop("overlap must be a numeric value >= 0 and < 1")
  }

  if (!isPosscal(nfft) || !isWhole(nfft)) {
    stop("nfft must be a positive integer")
  }

  if (!isPosscal(fs) || fs <= 0) {
    stop("fs must be a numeric value > 0")
  }

  detrend <- match.arg(detrend)
  range <- match.arg(range, c("half", "whole"))

  # initialize variables
  seg_len <- length(window)
  overlap <- trunc(seg_len * overlap)
  # GvB 20220405 removed nextpow2
  # nfft <- nextpow2(max(nfft, seg_len))
  nfft <- max(nfft, seg_len)
  win_meansq <- as.vector(window %*% window / seg_len)
  if (overlap >= seg_len) {
    stop("overlap must be smaller than windowlength")
  }
  # Pad data with zeros if shorter than segment. This should not happen.
  if (x_len < seg_len) {
    x <- c(x, rep(0, seg_len - x_len))
    x_len <- seg_len
  }

  # MAIN CALCULATIONS
  # remove overall mean or linear trend
  if (detrend == "long-mean" || detrend == "long-linear") {
    n_ffts <- max(0, trunc((x_len - seg_len) / (seg_len - overlap))) + 1
    x_len  <- min(x_len, (seg_len - overlap) * (n_ffts - 1) + seg_len)
    if (detrend == "long-mean") {
      x <- detrend(x, p = 0)
    } else if (detrend == "long-linear") {
      x <- detrend(x, p = 1)
    }
  }

  # Calculate and accumulate periodograms
  xx <- Pxx <- matrix(0, nrow = nfft, ncol = ns)
  if (ns > 1) {
    Pxy <- phase <- matrix(0, nrow = nfft, ncol = ns * (ns - 1) / 2)
  }
  n_ffts <- 0
  for (start_seg in seq(1, x_len - seg_len + 1, seg_len - overlap)) {
    end_seg <- start_seg + seg_len - 1
    # Don't truncate/remove the zero padding in xx
    if (detrend == "short-mean") {
      xx[1:seg_len, ] <- window * detrend(x[start_seg:end_seg, ], p = 0)
    } else if (detrend == "short-linear") {
      xx[1:seg_len, ] <- window * detrend(x[start_seg:end_seg, ], p = 1)
    } else {
      xx[1:seg_len, ] <- window * x[start_seg:end_seg, ]
    }
    fft_x <- stats::mvfft(xx)
    # accumulate periodogram
    Pxx <- Pxx + Re(fft_x * Conj(fft_x))
    # acculumulate crossspectrum for all signals
    if (ns > 1) {
      ptr <- 1
      for (i in seq_len(ns - 1)) {
        for (j in seq(i + 1, ns)) {
          # corrected bug 20211021
          # ptr <- i + (j - 1) * (j - 2) / 2
          Pxy[, ptr] <- Pxy[, ptr] + fft_x[, i] * Conj(fft_x[, j])
          ptr <- ptr + 1
        }
      }
    }
    n_ffts <- n_ffts + 1
  }

  # Convert two-sided spectra to one-sided spectra if range == 'half'
  # (normally if the input is numeric). For one-sided spectra, contributions
  # from negative frequencies are added to the positive side of the spectrum
  # -- but not at zero or Nyquist (half sampling) frequencies.
  # This keeps power equal in time and spectral domains, as required by
  # Parseval theorem.
  if (range == "half") {
    if (nfft %% 2 == 0) {    # one-sided, nfft is even
      psd_len <- nfft / 2 + 1
      Pxx <- apply(Pxx, 2,
                   function(x)
                     x[1:psd_len] + c(0, x[seq(nfft, psd_len + 1, -1)], 0))
      if (ns > 1) {
        Pxy <- apply(Pxy, 2,
                     function(x)
                       x[1:psd_len] +
                       Conj(c(0, x[seq(nfft, psd_len + 1, -1)], 0)))
      }
    } else {                    # one-sided, nfft is odd
      psd_len <- (nfft + 1) / 2
      Pxx <- apply(Pxx, 2,
                   function(x)
                     x[1:psd_len] + c(0, x[seq(nfft, psd_len + 1, -1)]))
      if (ns > 1) {
        Pxy <- apply(Pxy, 2,
                    function(x)
                      x[1:psd_len] + Conj(c(0, x[seq(nfft, psd_len + 1, -1)])))
      }
    }
  } else {  # range equals 'whole'
    psd_len <- nfft
  }
  # end MAIN CALCULATIONS

  ## SCALING AND OUTPUT
  scale <- n_ffts * seg_len * fs * win_meansq
  spec <- Pxx / scale
  colnames(spec) <- snames
  if (ns > 1) {
    cross <- Mod(Pxy) / scale
    coh <- phase <- trans <- matrix(0, nrow = psd_len, ncol = ns * (ns - 1) / 2)
    cn <- NULL
    ptr <- 1
    for (i in seq_len(ns - 1)) {
      for (j in seq(i + 1, ns)) {
        # corrected bug 20211021
        # ptr <- i + (j - 1) * (j - 2) / 2
        coh[, ptr] <- Mod(Pxy[, ptr])^2 / Pxx[, i] / Pxx[, j]
        phase[, ptr] <- Arg(Pxy[, ptr])
        trans[, ptr] <- Pxy[, ptr] / Pxx[, i]
        cn <- c(cn, paste(snames[i], snames[j], sep = "-"))
        ptr <- ptr + 1
      }
    }
    colnames(phase) <- colnames(cross)  <-
      colnames(trans) <- colnames(coh) <- cn
  } else {
    cross <- coh <- phase <- trans <- NULL
  }
  freq <- seq.int(0, psd_len - 1) * (fs / nfft)

  if (vec) {
    spec <- as.vector(spec)
  }
  all <- list(freq = freq, spec = spec, cross = cross, phase = phase,
              coh = coh, trans = trans, x_len = x_len, seg_len = seg_len,
              psd_len = psd_len, nseries = ns, series = series,
              snames = snames, window = window, fs = fs, detrend = detrend)
  class(all) <- "pwelch"
  all
}

#' @rdname pwelch
#' @export

plot.pwelch <- function(
    x, xlab = NULL, ylab = NULL, main = NULL,
    plot.type = c("spectrum", "cross-spectrum",
                  "phase", "coherence", "transfer"),
    yscale = c("linear", "log", "dB"), ...) {

  if (!inherits(x, "pwelch")) {
    stop("invalid object type")
  }
  plot.type <- match.arg(plot.type)
  yscale <- match.arg(yscale)

  if (is.null(xlab)) {
    if (x$fs == 1) {
      xlab <- expression(paste("Normalized frequency (\u00D7 ", pi,
                               " rad/sample)"))
    } else if (x$fs == pi) {
      xlab <- "Frequency (rad/sample)"
    } else {
      xlab <- "Frequency (Hz)"
    }
  }
  sub <- paste("Resolution:", format(x$fs / x$psd_len, digits = 6, nsmall = 6))

  if (plot.type == "spectrum" || plot.type == "cross-spectrum") {
    if (is.null(ylab)) {
      ylab <- switch(yscale,
                     "linear" = "Power/Frequency",
                     "log" = expression(paste("log"[10], "(Power/Frequency)")),
                     "dB" = "Power/Frequency (dB)")
    }
    if (plot.type == "spectrum") {
      if (is.null(main)) {
        main <- paste("Welch Power Spectral Density Estimate\nSeries:",
                      x$series)
      }
      plt <- switch(yscale,
                    "linear" = x$spec,
                    "log" = log10(x$spec),
                    "dB" = 10 * log10(x$spec))
    }
    if (plot.type == "cross-spectrum") {
      if (is.null(main)) {
        main <- paste("Welch Cross Power Spectral Density Estimate\nSeries:",
                      x$series)
      }
      plt <- switch(yscale,
                    "linear" = x$cross,
                    "log" = log10(x$cross),
                    "dB" = 10 * log10(x$cross))
    }
  }
  if (plot.type == "phase") {
    if (is.null(main)) {
      main <- paste("Welch Cross Spectrum Phase\nSeries:", x$series)
    }
    if (is.null(ylab)) {
      ylab <- expression(paste(theta, " / Frequency"))
    }
    plt <- x$phase
  }
  if (plot.type == "coherence") {
    if (is.null(main)) {
      main <- paste("Squared coherence\nSeries:", x$series)
    }
    if (is.null(ylab)) {
      ylab <- "Magnitude Squared Coherence"
    }
    plt <- x$coh
  }
  if (plot.type == "transfer") {
    if (is.null(main)) {
      main <- paste("Welch Transfer Function Estimate\nSeries:", x$series)
    }
    if (is.null(ylab)) {
      ylab <- switch(yscale,
                     "linear" = "Magnitude",
                     "log" = expression(paste("log"[10], "(Magnitude)")),
                     "dB" = "Magnitude (dB)")
    }
    plt <- switch(yscale,
                  "linear" = Mod(x$trans),
                  "log" = log10(Mod(x$trans)),
                  "dB" = 10 * log10(Mod(x$trans)))

  }
  graphics::matplot(x$freq, plt, type = "l", xlab = xlab, ylab = "", ...)
  graphics::title(main = main, sub = sub)
  graphics::title(ylab = ylab, line = 2)

  if (plot.type == "phase") {
    graphics::abline(h = c(-pi, -pi / 2, 0, pi / 2, pi), col = "red", lty = 2)
    graphics::text(-1, -3, expression(-pi), col = "red")
    graphics::text(-1, -1.4, expression(-pi / 2), col = "red")
    graphics::text(-1, 1.4, expression(pi / 2), col = "red")
    graphics::text(-1, 3, expression(pi), col = "red")
  }
}

#' @rdname pwelch
#' @export

print.pwelch <-
  function(x, plot.type = c("spectrum", "cross-spectrum",
                            "phase", "coherence", "transfer"),
           yscale = c("linear", "log", "dB"),
           xlab = NULL, ylab = NULL, main = NULL, ...) {
    plot.pwelch(x, plot.type, yscale, xlab, ylab, main, ...)
  }
