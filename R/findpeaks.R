# findpeaks.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Octave signal package:
# Copyright (C) 2012 Juan Pablo Carbajal <carbajal@ifi.uzh.ch>
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
# 20200105  GvB       setup for gsignal v0.1.0
#------------------------------------------------------------------------------

#' Find local extrema
#'
#' Return peak values and their locations of the vector \code{data}.
#'
#' Peaks of a positive array of \code{data} are defined as local maxima. For
#' double-sided data, they are maxima of the positive part and minima of the
#' negative part. \code{data} is expected to be a one-dimensional vector.
#'
#' @param data the data, expected to be a vector or one-dimensional array.
#' @param MinPeakHeight Minimum peak height (non-negative scalar). Only peaks
#'   that exceed this value will be returned. For data taking positive and
#'   negative values use the option \code{DoubleSided}. Default:
#'   \code{.Machine$double.eps}.
#' @param MinPeakDistance Minimum separation between peaks (positive integer).
#'   Peaks separated by less than this distance are considered a single peak.
#'   This distance is also used to fit a second order polynomial to the peaks to
#'   estimate their width, therefore it acts as a smoothing parameter. The
#'   neighborhood size is equal to the value of \code{MinPeakDistance}. Default:
#'   1.
#' @param MinPeakWidth Minimum width of peaks (positive integer). The width of
#'   the peaks is estimated using a parabola fitted to the neighborhood of each
#'   peak. The width is calculated with the formula \eqn{a * (width - x0)^{2} =
#'   1}, where a is the the concavity of the parabola and x0 its vertex.
#'   Default: 1.
#' @param MaxPeakWidth Maximum width of peaks (positive integer). Default:
#'   \code{Inf}.
#' @param DoubleSided Tells the function that data takes positive and negative
#'   values. The baseline for the peaks is taken as the mean value of the
#'   function. This is equivalent as passing the absolute value of the data
#'   after removing the mean. Default: FALSE
#'
#' @return A list containing the following elements:
#' \describe{
#'   \item{pks}{The value of data at the peaks.}
#'   \item{loc}{The index indicating the position of the peaks.}
#'   \item{parabol}{A list containing the parabola fitted to each returned peak.
#'   The list has two fields, \code{x} and \code{pp}. The field \code{pp}
#'   contains the coefficients of the 2nd degree polynomial and \code{x} the
#'   extrema of the interval where it was fitted.}
#'   \item{height}{The estimated height of the returned peaks (in units of
#'   data).}
#'   \item{baseline}{The height at which the roots of the returned peaks were
#'   calculated (in units of data).}
#'   \item{roots}{The abscissa values (in index units) at which the parabola
#'   fitted to each of the returned peaks realizes its width as defined below.}
#' }
#'
#' @examples
#' ### demo 1
#' t <- 2 * pi * seq(0, 1,length = 1024)
#' y <- sin(3.14 * t) + 0.5 * cos(6.09 * t) +
#'      0.1 * sin(10.11 * t + 1 / 6) + 0.1 * sin(15.3 * t + 1 / 3)
#'
#' data1 <- abs(y) # Positive values
#' peaks1 <- findpeaks(data1)
#'
#' data2 <- y # Double-sided
#' peaks2 <- findpeaks(data2, DoubleSided = TRUE)
#' peaks3 <- findpeaks (data2, DoubleSided = TRUE, MinPeakHeight = 0.5)
#'
#' op <- par(mfrow=c(1,2))
#' plot(t, data1, type="l", xlab="", ylab="")
#' points(t[peaks1$loc], peaks1$pks, col = "red", pch = 1)
#' plot(t, data2, type = "l", xlab = "", ylab = "")
#' points(t[peaks2$loc], peaks2$pks, col = "red", pch = 1)
#' points(t[peaks3$loc], peaks3$pks, col = "red", pch = 4)
#' legend ("topleft", "0: >2*sd, x: >0.5", bty = "n",
#'         text.col = "red")
#' par (op)
#' title("Finding the peaks of smooth data is not a big deal")
#'
#' ## demo 2
#' t <- 2 * pi * seq(0, 1, length = 1024)
#' y <- sin(3.14 * t) + 0.5 * cos(6.09 * t) + 0.1 *
#'      sin(10.11 * t + 1 / 6) + 0.1 * sin(15.3 * t + 1 / 3)
#' data <- abs(y + 0.1*rnorm(length(y),1))   # Positive values + noise
#' peaks1 <- findpeaks(data, MinPeakHeight=1)
#' dt <- t[2]-t[1]
#' peaks2 <- findpeaks(data, MinPeakHeight=1, MinPeakDistance=round(0.5/dt))
#' op <- par(mfrow=c(1,2))
#' plot(t, data, type="l", xlab="", ylab="")
#' points (t[peaks1$loc],peaks1$pks,col="red", pch=1)
#' plot(t, data, type="l", xlab="", ylab="")
#' points (t[peaks2$loc],peaks2$pks,col="red", pch=1)
#' par (op)
#' title(paste("Noisy data may need tuning of the parameters.\n",
#'             "In the 2nd example, MinPeakDistance is used\n",
#'             "as a smoother of the peaks"))
#'
#' @author Juan Pablo Carbajal, \email{carbajal@@ifi.uzh.ch}.\cr
#' Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#
#' @export

findpeaks <- function(data,
                      MinPeakHeight = .Machine$double.eps,
                      MinPeakDistance = 1,
                      MinPeakWidth = 1,
                      MaxPeakWidth = Inf,
                      DoubleSided = FALSE) {

  # check function arguments
  ld <- length(data)
  if (!is.numeric(data) ||
      !(is.vector(data) || is.array(data) || "ts" %in% class(data))
      || ld < 3)
    stop("data must be a numeric vector of at least 3 elements")
  if (!isPosscal(MinPeakHeight))
    stop("MinPeakHeight must be a positive scalar")
  if (!isPosscal(MinPeakDistance))
    stop("MinPeakDistance must be a positive scalar")
  if (!isPosscal(MinPeakWidth))
    stop("MinPeakWidth must be a positive scalar")
  if (!is.logical(DoubleSided))
    stop("DoubleSided should a a logical value TRUE or FALSE")

  wdata <- abs(detrend(data, 0))
  if (DoubleSided) {
    tmp <- data
    data <- wdata
    wdata <- tmp
  } else {
    if (min(data, na.rm = TRUE) < 0) {
      stop("Data contains negative values. Use the 'DoubleSided' option?")
    }
  }

  # Rough estimates of first and second derivative
  df1 <- diff(data, differences = 1)[c(1, 1:(ld - 1))]
  df2 <- diff(data, differences = 2)[c(1, 1, 1:(ld - 2))]

  # check for changes of sign of 1st derivative and negativity of 2nd deriv.
  # <= in 1st derivative includes the case of oversampled signals.
  idx <- which(df1 * c(df1[2:length(df1)], 0) <= 0 &
                 c(df2[2:length(df2)], 0) < 0)

  # Get peaks that are beyond given height
  tf  <- which(data[idx] > MinPeakHeight)
  idx <- idx[tf]
  if (length(idx) <= 0) return(NULL)

  # sort according to magnitude
  tmp <- sort(data[idx], decreasing = TRUE, index = TRUE)
  idx_s <- idx[tmp$ix]

  ## Treat peaks separated less than MinPeakDistance as one
  D <- with(expand.grid(A = idx_s, B = t(idx_s)), abs(A - B))
  dim(D) <- c(length(idx_s), length(idx_s))
  diag(D) <- NA                     # eliminate diagonal comparison
  if (any(D) < MinPeakDistance) {
    i <- 1
    node2visit <- seq_along(idx_s)
    visited <- NULL
    idx_pruned <- idx_s
    while (length(node2visit) > 0) {
      d <- D[node2visit[1], ]
      visited <- c(visited, node2visit[1])
      node2visit <- node2visit[-1]
      neighs <- setdiff(which(d < MinPeakDistance), visited)
      if (length(neighs) > 0) {
        idx_pruned <- setdiff(idx_pruned, idx_s[neighs])
        visited <- c(visited, neighs)
        node2visit <- setdiff(node2visit, visited)
      }
    }
    idx <- idx_pruned
  }
  idx <- sort(idx)

  extra_x <- extra_pp <- extra_roots <-
    extra_height <- extra_baseline <- data.frame()

  # Estimate widths of peaks and filter for:
  # width smaller than given.
  # wrong concavity.
  # not high enough
  # data at peak is lower than parabola by 1%

  idx.pruned <- idx
  n  <- length(idx)
  for (i in 1:n) {
    ind <- round(max(idx[i] - MinPeakDistance / 2, 1)) :
      round(min(idx[i] + MinPeakDistance / 2, ld))
    pp <- rep(0L, 3)
    if (any(data[ind] > data[idx[i]])) {
      pp <- pracma::polyfit(ind, data[ind], 2)
      xm <- -pp[2]^2 / (2 * pp[1])       # position of extrema
      H <- pracma::polyval(pp, xm)       # value at extrema
    } else {                             # use it as vertex of parabola
      H <- data[idx[i]]
      xm <- idx[i]
      pp <- rep(1L, 3)
      pp[1] <- pracma::mldivide((ind - xm)^2, (data[ind] - H))
      pp[2] <- -2 * pp[1] * xm
      pp[3] <- H + pp[1] * xm^2
    }
    width <- sqrt(abs(1 / pp[1])) + xm

    if ((width > MaxPeakWidth || width < MinPeakWidth) ||
        pp[1] > 0 || H < MinPeakHeight ||
        data[idx[i]] < 0.99 * H || abs(idx[i] - xm) > MinPeakDistance / 2) {
      idx_pruned <- setdiff(idx_pruned, idx[i])
    } else {
      extra_x <- rbind(extra_x, ind[c(1, length(ind))])
      extra_pp <- rbind(extra_pp, pp)
      extra_roots <- rbind(extra_roots, xm + c(-width, width) / 2)
      extra_height <- rbind(extra_height, H)
      extra_baseline <- rbind(extra_baseline, mean(c(H, MinPeakHeight)))
    }
  }
  idx <- idx.pruned

  # check for double sided
  if (DoubleSided) {
    pks <- wdata[idx]
  } else {
    pks <- data[idx]
  }

  # return values
  colnames(extra_x) <- c("from", "to")
  colnames(extra_pp) <- c("b2", "b", "a")
  colnames(extra_roots) <- c("a0", "a1")
  list(pks = pks, loc = idx,
       parabol = list(x = as.list(extra_x), pp = as.list(extra_pp)),
       height = as.numeric(extra_height[, 1]),
       baseline = as.numeric(extra_baseline[, 1]),
       roots = as.list(extra_roots))
}
