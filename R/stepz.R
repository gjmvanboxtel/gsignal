# stepz.R
# Copyright (C) 2026 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Original Octave function:
# Copyright (C) 2026 Tang Chonghao <chadholton@qq.com>
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
# 20261012  GvB       setup for gsignal v0.4.0
#------------------------------------------------------------------------------

#' Step response of digital filter
#'
#' \code{stepz} filters a length n step sequence using \code{fftfilt(b, a,
#' rep(1, n))} in the FIR case, or else using \code{filter(b,a,rep(1, n))}. A
#' plot is generated when the function is called without return value. If
#' \code{n} is not specified, \code{stepz} computes it using
#' \code{impzlength()}.
#' 
#' @param filt for the default case, the moving-average coefficients of an ARMA
#'   model or filter. Generically, \code{filt} specifies an arbitrary model or
#'   filter operation.
#' @param a the autoregressive (recursive) coefficients of an ARMA filter.
#'   Default: 1
#' @param n	number of points over which to evaluate the response. If
#'   \code{n} is a vector with a length greater than 1, then evaluate the
#'   frequency response at these points. For fastest computation, \code{n}
#'   should factor into a small number of small primes. Default: computed
#'   using \code{impzlength()}
#' @param fs sampling frequency in Hz. Default: 1
#' @param x	object to be plotted.
#' @param t times at which step response is plotted
#' @param h step response
#' @param ...	 for methods of \code{stepz}, arguments are passed to the default
#'   method. For \code{stepz_plot}, additional arguments are passed through to
#'   plot.
#'
#' @return For \code{stepz}, a list of class \code{'stepz'} with items:
#' \describe{
#'   \item{h}{step response at times \code{t}.}
#'   \item{t}{array of sample times.}
#' }
#'
#' @examples
#' 
#' # FIR
#' fir <- fir1(40, 0.3)
#' freqz(fir)
#' stepz(fir)
#' 
#' # custom number of points
#' stepz (1, c(1, -0.9), 30)
#'
#' # Step response at specific sample indices
#' stepz (1, c(1, -0.9), c(0, 5, 10, 15, 20))
#' 
#' # elliptical filter
#' el <- ellip(4, 0.5, 20, 0.4, output = "Sos")
#' stepz(el)
#' 
#' # 40th order low-pass Chebyshev type II filter
#' ch <- cheby2(40, 50, 0.4, output = "Zpg")
#' stepz(ch, n = 64)
#' 
#' @author Copyright (C) 2026 Tang Chonghao \email{chadholton@@qq.com}.\cr
#'  Port to R by Geert van Boxtel, \email{gjmvanboxtel@@gmail.com}
#'  
#' @seealso [freqz()], [zplane()], [impz()], [impzlength()]
#'
#' @rdname stepz
#' @export

stepz <- function(filt, ...) UseMethod("stepz")

#' @rdname stepz
#' @export
stepz.default <- function(filt, a = 1, n = gsignal::impzlength(filt, a),
                          fs = 1, ...)  {

  if (!(is.vector(filt) && is.vector(a))) {
    stop("filt and a must be vectors")
  }
  if (anyNA(filt) || anyNA(a)) {
    stop("filt and a must not contain missing values")
  } else {
    b <- filt
  }
  if (!(isScalar(fs) && fs > 0)) {
    stop("fs must be a scalar > 0")
  }

  n_is_vector <- FALSE
  if (missing(n)) {
    ## Auto-compute length using impzlength
    # should not occur because of default
    n <- impzlength(b, a)
  } else if (isScalar(n)) {
    if (!isPosscal(n)) {
      stop("n must be a non-negative integer when specified as a scalar")
    }
  } else {
    if (!all(sapply(n, isPosscal))) {
      stop("n must be all non-negative integers when specified as a vector")
    }
    n_is_vector <- TRUE
  }

  ## Compute step response
  if (n_is_vector) {
    max_n <- max(n)
    if (length(a) == 1) {
      x_full <- gsignal::fftfilt(b / a, rep(1, max_n + 1))
    } else {
    x_full = gsignal::filter(b, a, rep(1, max_n + 1))
    }
    ## Extract the step response at the specified indices
    x <- x_full[n + 1]
    t <- n / fs
  } else {
    if (length(a) == 1) {
      x <- gsignal::fftfilt(b / a, rep(1, n))
    } else {
      x <- gsignal::filter(b, a, rep(1, n))
    }
    t <- seq(0, n - 1, 1) / fs
  }

  res <- list(h = x, t = t)
  class(res) <- "stepz"
  res
}

#' @rdname stepz
#' @export

stepz.Arma <- function(filt, ...)
  stepz.default(filt$b, filt$a, ...)

#' @rdname stepz
#' @export

stepz.Ma <- function(filt, ...) # FIR
  stepz.default(unclass(filt), ...)

#' @rdname stepz
#' @export

stepz.Sos <- function(filt, ...) {
  stepz.Arma(as.Arma(filt), ...)
}

#' @rdname stepz
#' @export

stepz.Zpg <- function(filt, ...)
  stepz.Arma(as.Arma(filt), ...)

#' @rdname stepz
#' @export

print.stepz <- plot.stepz <- function(x, ...)
  stepz_plot(x$t, x$h, ...)

#' @rdname stepz
#' @export


#' @rdname stepz
#' @export

stepz_plot <- function(t, h, ...) {

  graphics::plot(t, h, type = "h", xlab = "Time or samples",
                 ylab = "Amplitude", ...)
  graphics::points(t, h)
  graphics::title("Step Response")
}

