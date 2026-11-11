# firpm.R - Octave Signal 1.4.8-compatible Parks-McClellan FIR design
# Copyright (C) 2026 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Octave signal package:
# Copyright (C) 2014/15/21 Rob Sykes <robs@users.sourceforge.net>#
# The native mmfir implementation retains its original GPL-3-or-later terms.
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
# 20261111  Geert van Boxtel      First version for v0.4.0
#-------------------------------------------------------------------------------

#' Parks-McClellan Optimal FIR Filter Design
#'
#' Design a linear-phase finite impulse response (FIR) digital filter
#' using the Parks-McClellan minimax approximation algorithm.
#'
#' The implementation uses the native `mmfir` algorithm adapted from
#' GNU Octave Signal 1.4.8 and supports both numeric frequency-response
#' specifications and user-defined R response functions.
#'
#' @param n
#'   Positive integer specifying the filter order. The resulting FIR
#'   filter has `n + 1` coefficients.
#'
#' @param f
#'   Numeric vector of normalized frequency-band edges.
#'
#'   Frequencies range from `0` (DC) to `1` (Nyquist).
#'   The vector must contain an even number of elements, representing
#'   successive lower and upper band edges, and must be strictly
#'   increasing.
#'
#'   For example, `c(0, 0.3, 0.4, 1)` specifies a passband and a
#'   stopband separated by a transition region.
#'
#' @param a
#'   Desired frequency response, specified as either:
#'  \itemize{
#'    \item A numeric vector containing one amplitude per frequency band.
#'    \item A numeric vector containing one amplitude per band edge.
#'    \item An R function defining a custom frequency response.
#'  }
#'
#'   A response function must accept the arguments
#'   `function(n, f, g, w, ...)`, where:
#'   \itemize{
#'    \item `n` is the filter order.
#'    \item `f` contains the original normalized band edges.
#'    \item `g` contains the normalized frequencies at which the response
#'     must be evaluated.
#'    \item `w` contains the weights supplied to `firpm()`.
#'    \item `...` contains optional arguments supplied through
#'     `response_args`.
#'   }
#'   The response function must return either:
#'   \itemize{
#'    \item A numeric vector of desired amplitudes, with length equal to
#'     `length(g)`.
#'    \item A list containing the desired amplitudes and optionally the
#'     corresponding positive weights.
#'   }
#'   For list results, the first element contains the amplitudes and
#'   the second element, when present, contains the weights.
#'   List element names such as `desired` and `weight` are recommended
#'   for readability.
#'
#'   All returned values must be finite. Weights must be strictly
#'   positive.
#'
#'   When the callback returns amplitudes without weights, unit
#'   weights are used.
#'
#' @param w
#'   Optional numeric vector of strictly positive weights.
#'
#'   For numeric response specifications, weights may be supplied
#'   per frequency band or per band edge.
#'
#'   The default is equal weighting.
#'
#'   For a custom response function, `w` is passed to the function.
#'   The callback must return the desired grid weights explicitly
#'   if non-unit weighting is required.
#'
#' @param class
#'   Character string specifying the filter class:
#'  \itemize{
#'    \item `"symmetric"`: symmetric impulse response.
#'    \item `"antisymmetric"`: antisymmetric impulse response.
#'    \item `"hilbert"`: Hilbert-transformer design.
#'    \item `"differentiator"`: differentiator design.
#'  }
#'   The default is `"symmetric"`.
#'
#' @param control
#'   Numeric vector containing one to three nonnegative algorithm
#'   control parameters:
#'  \itemize{
#'    \item `control[1]`: accuracy.
#'    \item `control[2]`: persistence.
#'    \item `control[3]`: robustness.
#'  }
#'   Missing elements are set to zero.
#'
#'   The default is `c(0, 0, 0)`.
#'
#'   Accuracy values greater than or equal to `16` are interpreted
#'   using the compatibility conversion implemented by this wrapper:
#'   `log(accuracy) / log(4) - 2`.
#'
#' @param output
#'   Character string specifying the output format:
#'  \itemize{
#'    \item `"coefficients"`: return the FIR coefficients.
#'    \item `"full"`: return the coefficients, optimization status,
#'     minimax error, and diagnostic information.
#'  }
#'   The default is `"coefficients"`.
#'
#' @param response_args
#'   List of additional arguments passed to a custom response
#'   function after its four standard arguments.
#'
#'   Named and unnamed list elements are supported.
#'
#'   This argument is only valid when `a` is a function.
#'
#' @details
#'   The Parks-McClellan algorithm constructs an FIR filter by
#'   minimizing the maximum weighted approximation error over the
#'   specified frequency bands.
#'
#'   The resulting filter has `n + 1` coefficients.
#'
#'   For numeric specifications, amplitudes and weights are
#'   interpolated over their corresponding frequency bands.
#'
#'   For custom response functions, the native algorithm evaluates
#'   the supplied function on its frequency grid. The callback
#'   determines the desired response and weighting at those points.
#'
#'   The optimization may fail to converge for some specifications,
#'   including degenerate designs whose optimal approximation error
#'   is exactly zero.
#'
#'   For example, a constant desired response over the complete
#'   frequency interval can be represented exactly by a centered
#'   impulse in a suitable symmetric FIR filter. Such a specification
#'   may cause the exchange algorithm to fail.
#'
#'   The native implementation is derived from the GPL-licensed
#'   GNU Octave Signal `mmfir` implementation.
#'
#' @return
#'   If `output = "coefficients"`, a numeric vector of length
#'   `n + 1` containing the FIR coefficients.
#'
#'   If `output = "full"`, a list containing:
#'  \itemize{
#'    \item `b`: FIR coefficients.
#'    \item `minimax`: minimax approximation error reported by the
#'     optimization algorithm.
#'    \item `result`: integer optimization status.
#'    \item `message`: character description of the optimization status.
#'    \item `iterations`: number of optimization iterations.
#'    \item `res`: list of diagnostic information containing:
#'      \itemize{
#'        \item `fgrid`: normalized frequency grid.
#'        \item `des`: desired response on the frequency grid.
#'        \item `wt`: weights on the frequency grid.
#'        \item `H`: complex frequency response on the grid.
#'        \item `error`: approximation error on the grid.
#'        \item `iextr`: one-based indices of extremal frequencies.
#'        \item `fextr`: normalized extremal frequencies.
#'      }
#'  }
#'   A successful optimization returns `result = 0`.
#'
#' @examples
#' # Low-pass FIR filter
#' b <- firpm(
#'   30,
#'   c(0, 0.3, 0.4, 1),
#'   c(1, 0)
#' )
#'
#' # Low-pass filter with unequal band weights
#' b <- firpm(
#'   30,
#'   c(0, 0.3, 0.4, 1),
#'   c(1, 0),
#'   w = c(1, 10)
#' )
#'
#' # Custom linear desired response
#' linear_response <- function(n, f, g, w) {
#'   list(
#'     desired = 1 - 0.5 * g,
#'     weight = rep(1, length(g))
#'   )
#' }
#'
#' b <- firpm(
#'   30,
#'   c(0, 1),
#'   linear_response
#' )
#'
#' # Custom nonlinear response and frequency-dependent weighting
#' nonlinear_response <- function(n, f, g, w) {
#'   list(
#'     desired = exp(-2 * g^2),
#'     weight = 1 + 5 * g
#'   )
#' }
#'
#' z <- firpm(
#'   30,
#'   c(0, 1),
#'   nonlinear_response,
#'   output = "full"
#' )
#'
#' # Inspect optimization results
#' z$minimax
#' z$iterations
#' z$res$fextr
#'
#' # Additional callback arguments
#' gaussian_response <- function(n, f, g, w, alpha) {
#'   list(
#'     desired = exp(-alpha * g^2),
#'     weight = rep(1, length(g))
#'   )
#' }
#'
#' b <- firpm(
#'   30,
#'   c(0, 1),
#'   gaussian_response,
#'   response_args = list(alpha = 2)
#' )
#'
#' @seealso
#'   \code{\link{fir1}}, \code{\link{freqz}}
#'
#' @export

firpm <- function(n, f, a, w = NULL,
                  class = c("symmetric", "antisymmetric", "hilbert",
                            "differentiator"),
                  control = c(0, 0, 0),
                  output = c("coefficients", "full"),
                  response_args = list()) {

  class <- match.arg(class)
  output <- match.arg(output)

  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) ||
      n < 1 || n != floor(n) || n >= .Machine$integer.max)
    stop("n must be a positive integer")

  if (!is.numeric(f) || length(f) < 2L || length(f) %% 2L ||
      any(!is.finite(f)) || any(diff(f) <= 0) || any(f < 0 | f > 1))
    stop("f must be strictly increasing band edges in [0,1]")

  nband <- length(f) / 2L
  custom <- is.function(a)
  if (!custom &&
      (!is.numeric(a) || !(length(a) %in% c(nband, 2L * nband)) ||
       any(!is.finite(a))))
    stop("a must contain one amplitude per band or band edge, or be a function")

  if (is.null(w)) w <- rep(1, nband)
  if (!is.numeric(w) || !(length(w) %in% c(nband, 2L * nband)) ||
      any(!is.finite(w)) || any(w <= 0))
    stop("w must contain positive weights per band or per band edge")

  if (!is.list(response_args))
    stop("response_args must be a list")
  if (!custom && length(response_args))
    stop("response_args is only valid with a response function")

  if (!is.numeric(control) || length(control) < 1L ||
      length(control) > 3L || any(!is.finite(control)) ||
      any(control < 0))
    stop("control must contain 1 to 3 finite, nonnegative values")
  control <- c(control, rep(0, 3L - length(control)))
  if (control[1L] >= 16)
    control[1L] <- log(control[1L]) / log(4) - 2

  filterClass <- switch(class, symmetric = 0L, antisymmetric = 1L,
                        hilbert = 1L, differentiator = 2L)
  ans <- firpm_mmfir(
    as.integer(n), as.double(f),
    if (custom) numeric() else as.double(a),
    as.double(w), as.integer(filterClass),
    as.double(control[1]), as.double(control[2]), as.double(control[3]),
    output == "full",
    if (custom) a else NULL,
    response_args
  )
  if (output == "coefficients") return(ans$b)
  ans$res <- .firpm_diagnostics(ans, n, filterClass)
  ans[c("fgrid", "des_internal", "wt_internal", "peaks")] <- NULL
  ans
}


.firpm_diagnostics <- function(ans, n, filterClass) {
  fgrid <- ans$fgrid
  omega <- pi * fgrid
  des_internal <- ans$des_internal
  wt_internal <- ans$wt_internal
  peaks <- ans$peaks
  
  des <- if (filterClass == 2L) -des_internal else des_internal
  wt <- if (filterClass == 2L) 0.5 * wt_internal else wt_internal
  
  k <- seq_along(ans$b) - 1L
  H <- vapply(omega, function(om) {
    sum(ans$b * exp(-1i * om * k))
  }, complex(1))
  
  offset <- if (filterClass == 0L) 0 else pi / 2
  phase <- offset + n * omega / 2
  ampResp <- Re(H) * cos(phase) - Im(H) * sin(phase)
  error <- des_internal - ampResp
  
  iextr <- integer(length(peaks))
  for (i in seq_along(peaks)) {
    j <- which(fgrid >= peaks[i])[1L]
    if (is.na(j)) j <- length(fgrid)
    if (j > 1L && fgrid[j] != peaks[i] &&
        abs(error[j]) < abs(error[j - 1L])) j <- j - 1L
    iextr[i] <- j
  }
  
  list(fgrid = fgrid, des = des, wt = wt, H = H,
       error = error, iextr = iextr, fextr = fgrid[iextr])
}
