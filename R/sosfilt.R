# sosfilt.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
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
# 20200413  GvB       setup for gsignal v0.1.0
# 20210329  GvB       different setup for v0.3.0 including initial conditions
# 20210712  GvB       copy attributes of input x to output y
# 20210724  GvB       allow filtering complex signals
#------------------------------------------------------------------------------

#' Second-order sections filtering
#'
#' One-dimensional second-order (biquadratic) sections IIR digital filtering.
#'
#' The filter function is implemented as a series of second-order filters
#' with direct-form II transposed structure. It is designed to minimize
#' numerical precision errors for high-order filters [1].
#'
#' @param sos Second-order section representation, specified as an nrow-by-6
#'   matrix, whose rows contain the numerator and denominator coefficients of
#'   the second-order sections:\cr \code{sos <- rbind(cbind(B1, A1), cbind(...),
#'   cbind(Bn, An))}, where \code{B1 <- c(b0, b1, b2)}, and \code{A1 <- c(a0,
#'   a1, a2)} for section 1, etc. The b0 entry must be nonzero for each section.
#' @param x the input signal to be filtered, specified as a numeric or complex
#'   vector or matrix. If \code{x} is a matrix, each column is filtered.
#' @param zi If \code{zi} is provided, it is taken as the initial state of the
#'   system and the final state is returned as zf. If \code{x} is a vector,
#'   \code{zi} must be a matrix with \code{nrow(sos)} rows and 2 columns. If
#'   \code{x} is a matrix, then \code{zi} must be a 3-dimensional array of size
#'   \code{(nrow(sos), 2, ncol(x))}. Alternatively, \code{zi} may be the
#'   character string \code{"zf"}, which specifies to return the final state
#'   vector even though the initial state vector is set to all zeros. Default:
#'   NULL.
#'
#' @return The filtered signal, of the same dimensions as the input signal. In
#'   case the \code{zi} input argument was specified, a list with two elements
#'   is returned containing the variables \code{y}, which represents the output
#'   signal, and \code{zf}, which contains the final state vector or matrix.
#'
#' @examples
#' fs <- 1000
#' t <- seq(0, 1, 1/fs)
#' s <- sin(2* pi * t * 6)
#' x <- s + rnorm(length(t))
#' plot(t, x, type = "l", col="light gray")
#' lines(t, s, col="black")
#' sosg <- butter(3, 0.02, output = "Sos")
#' sos <- sosg$sos
#' sos[1, 1:3] <- sos[1, 1:3] * sosg$g
#' y <- sosfilt(matrix(sos, ncol=6), x)
#' lines(t, y, col="red")
#'
#' ## using 'filter' will handle the gain for you
#' y2 <- filter(sosg, x)
#' all.equal(y, y2)
#'
#' ## The following example is from Python scipy.signal.sosfilt
#' ## It shows the instability that results from trying to do a
#' ## 13th-order filter in a single stage (the numerical error
#' ## pushes some poles outside of the unit circle)
#' arma <- ellip(13, 0.009, 80, 0.05, output='Arma')
#' sos <- ellip(13, 0.009, 80, 0.05, output='Sos')
#' x <- rep(0, 700); x[1] <- 1
#' y_arma <- filter(arma, x)
#' y_sos <- filter(sos, x)
#' plot(y_arma, type ="l")
#' lines (y_sos, col = 2)
#' legend("topleft", legend = c("Arma", "Sos"), lty = 1, col = 1:2)
#'
#' @seealso \code{\link{filter}}, \code{\link{filtfilt}}, \code{\link{Sos}}
#'
#' @references Smith III, J.O. (2012). Introduction to digital filters, with
#'   audio applications (3rd Ed.). W3K Publishing.
#'
#' @author Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#'
#' @export

sosfilt <- function(sos, x, zi = NULL) {

  # Check sos
  if (is.vector(sos)) {
    if (length(sos) == 6) {
      sos <- matrix(sos, ncol = 6)
    } else {
      stop("sos must a matrix with 6 columns")
    }
  } else if (is.matrix(sos)) {
    if (ncol(sos) != 6) {
      stop("sos must a matrix with 6 columns")
    }
  } else {
    stop("sos must a matrix with 6 columns")
  }
  a0 <- sos[, 4]
  if (any(a0 == 0)) {
    stop("invalid sos structure (sos[, 4] must not be zero)")
  }
  if (any(a0 != 1)) {
    sos <- sos / a0
  }

  # check x, coerce it to matrix
  if (is.null(x)) {
    return(NULL)
  }
  if (!(is.numeric(x) || is.complex(x))) {
    stop("x must be a numeric or complex vector or matrix")
  }
  
  #save attributes of x
  atx <- attributes(x)
  if (is.vector(x)) {
    x <- as.matrix(x, ncol = 1)
    vec <- TRUE
  } else {
    vec <- FALSE
  }
  nrx <- NROW(x)
  ncx <- NCOL(x)
  if (is.null(nrx) || nrx <= 0) {
    return(x)
  }

  # check zi, coerce to 3d array
  rzf <- (is.character(zi) && zi == "zf")
  if (!is.null(zi) && !is.numeric(zi) && !rzf) {
    stop("zi must be NULL, a numeric vector or matrix, or the string 'zf'")
  }
  if (is.null(zi) || rzf) {
    if (is.numeric(x)) {
      zi <- array(0, dim = c(nrow(sos), 2, ncx))
    } else if(is.complex(x)) {
      zi <- array(0 + 0i, dim = c(nrow(sos), 2, ncx))
    }
    if (is.null(zi)) {
      rzf <- FALSE
    }
  } else if (!is.null(zi)) {
    rzf <- TRUE
  }
  if (is.vector(zi)) {
    stop("zi must be NULL, a matrix or a 3-dimensional array")
  } else {
    dims_zi <- dim(zi)
    if (length(dims_zi) == 2) {
      dim(zi) <- c(dims_zi, 1)
      dims_zi <- c(dims_zi, 1)
    }
  }
  if (dims_zi[1] != nrow(sos)) {
    stop("zi must equal the number of sections in sos")
  }
  if (dims_zi[2] != 2) {
    stop("number of columns of zi must be 2")
  }
  if (dims_zi[3] != ncx) {
    stop("third dimension of zi must be equal the number of columns in x")
  }

  if (is.numeric(x)) {
    y <- matrix(0, nrx, ncx)
    zf <- array(0, dim = dims_zi)
    for (icol in seq_len(ncx)) {
      l <- .Call("_gsignal_rsosfilt", PACKAGE = "gsignal",
                 sos, x[, icol], matrix(zi[, , icol], ncol = 2))
      if (length(l) <= 0) {
        stop("Error filtering data")
      }
      y[, icol] <- l[["y"]]
      zf[, , icol] <- l[["zf"]]
    }
  } else if (is.complex(x)) {
    y <- matrix(0 + 0i, nrx, ncx)
    zf <- array(0 + 0i, dim = dims_zi)
    for (icol in seq_len(ncx)) {
      re <- .Call("_gsignal_rsosfilt", PACKAGE = "gsignal",
                 sos, Re(x[, icol]), Re(matrix(zi[, , icol], ncol = 2)))
      if (length(re) <= 0) {
        stop("Error filtering data")
      }
      im <- .Call("_gsignal_rsosfilt", PACKAGE = "gsignal",
                  sos, Im(x[, icol]), Im(matrix(zi[, , icol], ncol = 2)))
      if (length(im) <= 0) {
        stop("Error filtering data")
      }
      y[, icol] <- re[["y"]] + 1i * im[["y"]]
      zf[, , icol] <- re[["zf"]] + 1i * im[["zf"]]
    }
  }
  if (vec) {
    y <- as.vector(y)
    dim(zf) <- dims_zi[1:2]
  }
  if (is.complex(x)) {
    y <- zapIm(y)
    zf <- zapIm(zf)
  }
  # set attributes of y nd return
  attributes(y) <- atx
  
  if (rzf) {
    rv <- list(y = y, zf = zf)
  } else {
    rv <- y
  }
  rv
}
