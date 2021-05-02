# mpoles.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Original Octave function:
# Copyright (C) 2007-2017 Ben Abbott
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
# 20200606  GvB       setup for gsignal v0.1.0
#------------------------------------------------------------------------------

#' Multiplicity of poles
#'
#' Identify unique poles and their associated multiplicity.
#'
#' @param p vector of poles.
#' @param tol tolerance. If the relative difference of two poles is less than
#'   \code{tol} then they are considered to be multiples. The default value for
#'   \code{tol} is 0.001.
#' @param reorder logical. If \code{TRUE}, (default), the output is ordered from
#'   largest pole to smallest pole.
#' @param index.return logical indicating if index vector should be returned as
#'   well. See examples. Default: \code{FALSE}.
#'
#' @return If \code{index.return = TRUE}, a list consisting of two vectors:
#' \describe{
#'   \item{m}{vector specifying the multiplicity of the poles}
#'   \item{n}{index}
#' }
#' If \code{index.return = FALSE}, only \code{m} is returned (as a vector).
#'
#' @examples
#' p <- c(2, 3, 1, 1, 2)
#' ret <- mpoles(p, index = TRUE)
#'
#' @seealso \code{\link{poly}}, \code{\link{residue}}
#'
#' @author Ben Abbott, \email{bpabbott@@mac.com}.\cr
#'   Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}
#
#' @export

mpoles <- function(p, tol = 0.001, reorder = TRUE, index.return = FALSE) {

  ## Force the poles to be a vector.
  p <- as.vector(p)
  np <- length(p)
  # tol must be a positive scalar
  tol <- tol[1]
  if (!isPosscal(tol)) {
    stop("'tol' must be a positive scalar")
  }
  tol <- abs(tol)
  # reorder and index.return shuuld be logical
  if (!is.logical(reorder)) {
    stop("'reorder' should be TRUE or FALSE")
  }
  if (!is.logical(index.return)) {
    stop("'index.return' should be TRUE or FALSE")
  }

  ## Sort the poles according to their magnitidues, largest first.
  if (reorder) {
    ## Sort with smallest magnitude first.
    s <- sort(p, index.return = TRUE)
    ## Reverse order, largest maginitude first.
    n <- seq(np, 1, -1)
    p <- s$x[n]
    ordr <- s$ix[n]
  } else {
    ordr <- seq_len(np)
  }

  ## Find pole multiplicty by comparing the relative differnce in the
  ## poles.

  multp <- array(0L, np)
  indx <- NULL
  n <- which(multp == 0)[1]
  while (!is.na(n) & n > 0) {
    dp <- abs(p - p[n])
    if (p[n] == 0.0) {
      if (any(abs(p) > 0 & is.finite(p))) {
        p0 <- mean(abs(p[abs(p) > 0 & is.finite(p)]))
      } else {
        p0 <- 1
      }
    } else {
      p0 <- abs(p[n])
    }
    k <- which(dp < tol * p0)
    ## Poles can only be members of one multiplicity group.
    if (length(indx)) {
      k <- k[!(k %in% indx)]
    }
    m <- seq_len(length(k))
    multp[k] <- m
    indx <- c(indx, k)
    n <- which(multp == 0)[1]
  }
  multp <- as.numeric(multp[indx])
  indx <- as.numeric(ordr[indx])

  if (index.return) {
    ret <- list(m = multp, n = indx)
  } else {
    ret <- multp
  }
  ret
}
