# buffer.R
# Copyright (C) 2019 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Matlab/Octave signal package:
# Copyright (C) 2008 David Bateman <adb014@gmail.com>
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
# 20191120 Geert van Boxtel          First version for v0.1.0
#------------------------------------------------------------------------------

#' Buffer signal vector into matrix of data segments
#'
#' Partition a signal vector into nonoverlapping, overlapping, or underlapping
#' data segments.
#'
#' \code{y <- buffer(x, n)} partitions a signal vector \code{x} of length
#' \code{L} into nonoverlapping data segments of length \code{n}. Each data
#' segment occupies one column of matrix output \code{y}, which has \code{n}
#' rows and \code{ceil(L / n)} columns. If \code{L} is not evenly divisible by
#' \code{n}, the last column is zero-padded to length \code{n}.
#'
#' \code{y <- buffer(x, n, p)} overlaps or underlaps successive frames in the
#' output matrix by \code{p} samples.
#' \itemize{
#' \item {For \code{0 < p < n} (overlap), buffer repeats the final \code{p}
#' samples of each segment at the beginning of the following segment. See the
#' example where \code{x = 1:30}, \code{n = 7}, and an overlap of \code{p = 3}.
#' In this case, the first segment starts with \code{p} zeros (the default
#' initial condition), and the number of columns in \code{y} is \code{ceil(L /
#' (n - p))}.}
#' \item  {For \code{p < 0} (underlap), buffer skips \code{p} samples between
#' consecutive segments. See the example where \code{x = 1:30}, \code{n = 7},
#' and \code{p = -3}. The number of columns in \code{y} is \code{ceil(L / (n -
#' p))}.}
#' }
#'
#' In \code{y <- buffer(x, n, p, opt)}, \code{opt} specifies a vector of samples
#' to precede \code{x[1]} in an overlapping buffer, or the number of initial
#' samples to skip in an underlapping buffer.
#' \itemize{
#'   \item {For \code{0 < p < n} (overlap), \code{opt} specifies a vector of
#'   length \code{p} to insert before \code{x[1]} in the buffer. This vector can
#'   be considered an initial condition, which is needed when the current
#'   buffering operation is one in a sequence of consecutive buffering
#'   operations. To maintain the desired segment overlap from one buffer to the
#'   next, \code{opt} should contain the final \code{p} samples of the previous
#'   buffer in the sequence. Set \code{opt} to \code{"nodelay"} to skip the
#'   initial condition and begin filling the buffer immediately with
#'   \code{x[1]}. In this case, \code{L} must be \code{length(p)} or longer. See
#'   the example where \code{x = 1:30}, \code{n = 7}, \code{p = 3}, and
#'   \code{opt = "nodelay"}.}
#'   \item {For \code{p < 0} (underlap), \code{opt} is an integer value in the
#'     range \code{0 : -p} specifying the number of initial input samples,
#'     \code{x[1:opt]}, to skip before adding samples to the buffer. The first
#'     value in the buffer is therefore \code{x[opt + 1]}.}
#' }
#' The \code{opt} option is especially useful when the current buffering
#' operation is one in a sequence of consecutive buffering operations. To
#' maintain the desired frame underlap from one buffer to the next, \code{opt}
#' should equal the difference between the total number of points to skip
#' between frames (\code{p}) and the number of points that were available to be
#' skipped in the previous input to buffer. If the previous input had fewer than
#' p points that could be skipped after filling the final frame of that buffer,
#' the remaining opt points need to be removed from the first frame of the
#' current buffer. See Continuous Buffering for an example of how this works in
#' practice.
#'
#' \code{buf <- buffer(..., zopt = TRUE)} returns the last \code{p} samples of a
#' overlapping buffer in output \code{buf$opt}. In an underlapping buffer,
#' \code{buf$opt} is the difference between the total number of points to skip
#' between frames (\code{-p}) and the number of points in \code{x} that were
#' available to be skipped after filling the last frame:
#' \itemize{
#'   \item {For \code{0 < p < n} (overlap), \code{buf$opt} contains the final
#'   \code{p} samples in the last frame of the buffer. This vector can be used
#'   as the initial condition for a subsequent buffering operation in a sequence
#'   of consecutive buffering operations. This allows the desired frame overlap
#'   to be maintained from one buffer to the next. See Continuous Buffering
#'   below.}
#'   \item {For \code{p < 0} (underlap), \code{buf$opt} is the difference
#'   between the total number of points to skip between frames \code{(-p)} and
#'   the number of points in \code{x} that were available to be skipped after
#'   filling the last frame: \code{buf$opt = m*(n-p) + opt - L} where \code{opt}
#'   on the right is the input argument to buffer, and \code{buf$opt} on the
#'   left is the output argument. Note that for an underlapping buffer output
#'   \code{buf$opt} is always zero when output \code{buf$z} contains data.\cr
#'   The opt output for an underlapping buffer is especially useful when the
#'   current buffering operation is one in a sequence of consecutive buffering
#'   operations. The \code{buf$opt} output from each buffering operation
#'   specifies the number of samples that need to be skipped at the start of the
#'   next buffering operation to maintain the desired frame underlap from one
#'   buffer to the next. If fewer than \code{p} points were available to be
#'   skipped after filling the final frame of the current buffer, the remaining
#'   opt points need to be removed from the first frame of the next buffer.}
#' }
#' In a sequence of buffering operations, the \code{buf$opt} output from each
#' operation should be used as the \code{opt} input to the subsequent buffering
#' operation. This ensures that the desired frame overlap or underlap is
#' maintained from buffer to buffer, as well as from frame to frame within the
#' same buffer. See Continuous Buffering below for an example of how this works
#' in practice.
#' \cr
#'
#' \strong{Continuous Buffering}\cr\cr
#' In a continuous buffering operation, the vector input to the buffer function
#' represents one frame in a sequence of frames that make up a discrete signal.
#' These signal frames can originate in a frame-based data acquisition process,
#' or within a frame-based algorithm like the FFT.\cr
#' As an example, you might acquire data from an A/D card in frames of 64
#' samples. In the simplest case, you could rebuffer the data into frames of 16
#' samples; \code{buffer} with \code{n = 16} creates a buffer of four frames
#' from each 64-element input frame. The result is that the signal of frame size
#' 64 has been converted to a signal of frame size 16; no samples were added or
#' removed.\cr
#' In the general case where the original signal frame size, \code{L}, is not
#' equally divisible by the new frame size, \code{n}, the overflow from the last
#' frame needs to be captured and recycled into the following buffer. You can do
#' this by iteratively calling buffer on input x with the \code{zopt} parameter
#' set to \code{TRUE}. This simply captures any buffer overflow in \code{buf$z},
#' and prepends the data to the subsequent input in the next call to buffer.\cr
#' Note that continuous buffering cannot be done without the \code{zopt}
#' parameter being set to \code{TRUE}, because the last frame of y (\code{buf$y}
#' in this case) is zero padded, which adds new samples to the signal.\cr
#' Continuous buffering in the presence of overlap and underlap is handled with
#' the \code{opt} parameter, which is used as both an input (\code{opt} and
#' output (\code{buf$opt}) to buffer. The two examples on this page demonstrate
#' how the \code{opt} parameter should be used.
#'
#' @param x The data to be buffered.
#' @param n The number of rows in the produced data buffer. This is an positive
#'   integer value and must be supplied.
#' @param p An integer less than \code{n} that specifies the under- or overlap
#'   between column in the data frame. Default 0.
#' @param opt In the case of an overlap, \code{opt} can be either a vector of
#'   length \code{p} or the string \code{'nodelay'}. If \code{opt} is a vector,
#'   then the first \code{p} entries in \code{y} will be filled with these
#'   values. If \code{opt} is the string \code{'nodelay'}, then the first value
#'   of \code{y} corresponds to the first value of \code{x}. In the case of an
#'   underlap, \code{opt} must be an integer between 0 and \code{-p}. The
#'   represents the initial underlap of the first \code{y}. The default value
#'   for \code{opt} the vector \code{matrix (0L, 1, p)} in the case of an
#'   overlap, or 0 otherwise.
#' @param zopt Logical. If TRUE, return values for \code{z} and \code{opt} in
#'   addition to \code{y}. Default is FALSE (return only \code{y}).
#'
#' @return If \code{zopt} equals FALSE (the default), this function returns a
#'   single numerical array containing the buffered data (\code{y}). If
#'   \code{zopt} equals TRUE, then a \code{list} containing 3 variables is
#'   returned: \code{y}: the buffered data, \code{z}: the over or underlap (if
#'   any), \code{opt}: the over- or underlap that might be used for a future
#'   call to \code{buffer} to allow continuous buffering.
#'
#' @examples
#' ## Examples without continuous buffering
#' y <- buffer(1:10, 5)
#' y <- buffer(1:10, 4)
#' y <- buffer(1:30, 7, 3)
#' y <- buffer(1:30, 7, -3)
#' y <- buffer(1:30, 7, 3, 'nodelay')
#'
#' ## Continuous buffering examples
#' # with overlap:
#' data <- buffer(1:1100, 11)
#' n <- 4
#' p <- 1
#' buf <- list(y = NULL, z = NULL, opt = -5)
#' for (i in 1:ncol(data)) {
#'   x <- data[,i]
#'   buf <- buffer(x = c(buf$z,x), n, p, opt=buf$opt, zopt = TRUE)
#' }
#' # with underlap:
#' data <- buffer(1:1100, 11)
#' n <- 4
#' p <- -2
#' buf <- list(y = NULL, z = NULL, opt = 1)
#' for (i in 1:ncol(data)) {
#'   x <- data[,i]
#'   buf <- buffer(x = c(buf$z,x), n, p, opt=buf$opt, zopt = TRUE)
#' }
#'
#' @author David Bateman, \email{adb014@@gmail.com}.\cr
#'  Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}
#'
#' @export

buffer <- function(x, n, p = 0, opt, zopt = FALSE) {

  # parameter checking etc.
  if (is.data.frame(x)) x <- as.vector(x[, 1])
  else x <- as.vector(x)
  if (!isScalar(n) || !isWhole(n)) stop("n must be an integer")
  if (!isScalar(p) || !isWhole(p) || p >= n)
    stop("p must be an integer less than n")
  if (missing(opt)) {
    if (p < 0) {
      opt <- 0
    } else {
      opt <- matrix(0L, 1, p)
    }
  }
  if (p < 0) {
    if (isScalar(opt) && isWhole(opt) && opt >= 0 && opt <= -p) {
      lopt <- opt
    } else {
      stop("expecting opt to be an integer between 0 and -p")
    }
  } else {
    lopt <- 0
  }
  if (!is.logical(zopt)) stop("zopt must be a logical")

  l <- length(x)
  m <- ceiling((l - lopt) / (n - p))
  y <- matrix(0L, n - p, m)
  y [1:(l - lopt)] <- x[(lopt + 1):l]

  if (p < 0) {
    y <- y[- ((nrow(y) + p + 1):nrow(y)), ]
  } else if (p > 0) {
    if (is.character(opt)) {
      if (opt == "nodelay") {
        y <- rbind(y, matrix(0L, p, m))
        if (p > n / 2) {
          iis <- n - p + 1
          iin <- n - p
          iie <- iis + iin - 1
          off <- 1
          while (iin > 0) {
            y[iis:iie, 1:(ncol(y) - off)] <- y[1:iin, (1 + off):ncol(y)]
            off <- off + 1
            iis <- iie + 1
            iie <- iie + iin
            if (iie > n) {
              iie <- n
            }
            iin <- iie - iis + 1
          }
          i <- ((l - 1) %% (n - p)) + 1
          j <- floor((l - 1) / (n - p)) + 1
          if (all(c(i, j) == c(n - p, m))) {
            off <- off - 1
          }
          y <- y[, -c((ncol(y) - off + 2):ncol(y))]
        } else {
          y[(nrow(y) - p + 1):nrow(y), 1:(ncol(y) - 1)] <- y[(1:p), 2:ncol(y)]
          if ((m - 1) * (n - p) + p >= l) {
            y <- y[, -ncol(y)]
          }
        }
      } else {
        stop(paste("Unexpected string argument to 'opt':", opt))
      }
    } else if (is.numeric(opt)) {
      if (length(opt) == p) {
        lopt <- p
        y <- rbind(matrix(0L, p, m), y)
        iin <- p
        off <- 1
        out <- 1
        while (iin > 0) {
          y[1:iin, off] <- opt[out:length(opt)]
          off <- off + 1
          iin <- iin - (n - p)
          out <- out + (n - p)
        }
        if (p > n / 2) {
          iin <- n - p
          iie <- p
          iis <- p - iin + 1
          off <- 1
          while (iie > 0) {
            y[iis:iie, (1 + off):ncol(y)] <-
              y[(nrow(y) - iin + 1):nrow(y), 1:(ncol(y) - off)]
            off <- off + 1
            iie <- iis - 1
            iis <- iis - iin
            if (iis < 1) {
              iis <- 1
            }
            iin <- iie - iis + 1
          }
        } else {
          y[1:p, 2:ncol(y)] <- y[(nrow(y) - p + 1):nrow(y), 1:(ncol(y) - 1)]
        }
      } else {
        stop("'opt' vector should be of length 'p'")
      }
    } else {
      stop("Unrecognized 'opt' argument")
    }
  }

  if (zopt) {
    if (p >= 0) {
      i <- (((l + lopt + p * (ncol(y) - 1)) - 1) %% nrow(y)) + 1
      j <- floor(((l + lopt + p * (ncol(y) - 1)) - 1) / nrow(y)) + 1
      if (any(c(i, j) != c(nrow(y), ncol(y)))) {
        z <- y[(1 + p):i, ncol(y)]
        y <- y[, -ncol(y)]
      } else {
        z <- NULL
      }
    } else {
      i <- ((l - lopt - 1) %% (nrow(y) - p)) + 1
      if (i < nrow(y)) {
        z <- y[1:i, ncol(y)]
        y <- y[, -ncol(y)]
      } else {
        z <- NULL
      }
    }
    if (p < 0) {
      opt <- max(0, ncol(y) * (n - p) + opt - l)
    } else if (p > 0) {
      opt <- y[(nrow(y) - p + 1):nrow(y), ncol(y)]
    } else {
      opt <- NA
    }
    return(list(y = y, z = z, opt = opt))
  } else {
    return(y)
  }
}
