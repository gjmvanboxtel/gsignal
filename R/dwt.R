# dwt.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Original Octave code:
# Copyright (C) 2013 Lukas F. Reichlin
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
# 20201020  GvB       setup for gsignal v0.1.0
#------------------------------------------------------------------------------

#' 1-D Discrete Wavelet Transform
#'
#' Compute the single-level discrete wavelet transform of a signal
#'
#' This function is only included because of compatibility with the 'Octave'
#' 'signal' package. Specialized packages exist in R to perform the discrete
#' wavelet transform, e.g., the \code{wavelets} package [1]. this function
#' recognizes only a few wavelet names, namely those for which scale
#' coefficients are available (Daubechies [2] and Coiflet [3]).
#'
#' The wavelet and scaling coefficients are returned by the function
#' \code{wfilters}, which returns the coefficients for reconstruction filters
#' associated with the wavelet \code{wname}. Decomposition filters are the time
#' reverse of the reconstruction filters (see examples).
#'
#' @param x input data, specified as a numeric vector.
#' @param wname analyzing wavelet, specified as a character string consisting of
#'   a class name followed by the wavelet length Only two classes of wavelets
#'   are supported; Daubechies (denoted by the prefix \code{'d'} of even
#'   lengths 2 - 20, and Coiflet (denoted by the prefix '\code{'c'} of
#'   lengths 6, 12, 18, 24, and 30. The wavelet name \code{'haar'} is
#'   the equivalent of \code{'d2'}. Default: d8.
#' @param lo scaling (low-pass) filter, specified as an even-length numeric
#'   vector. \code{lo} must be the same length as \code{hi}. Ignored when
#'   \code{wname != NULL}.
#' @param hi wavelet (high-pass) filter, specified as an even-length numeric
#'   vector. \code{hi} must be the same length as \code{lo}, Ignored when
#'   \code{wname != NULL}.
#'
#' @note The notations \code{g} and \code{h} are often used to denote low-pass
#'   (scaling) and high-pass (wavelet) coefficients, respectively, but
#'   inconsistently. Ref [4] uses it, as does the R \code{wavelets} package.
#'   'Octave' uses the reverse notation. To avoid confusion, more neutral terms
#'   are used here.
#'
#' @note There are two naming schemes for wavelet names in use. For instance for
#'   Daubechies wavelets (d), dN using the length or number of taps, and dbA
#'   referring to the number of vanishing moments. So d4 and db2 are the same
#'   wavelet transform. This function uses the formed (dN) notation; 'Matlab'
#'   uses the latter (dbA).
#'
#' @return A list containing two numeric vectors:
#' \describe{
#'   \item{a}{approximation (average) coefficients, obtained from convolving
#'     \code{x} with the scaling (low-pass) filter \code{lo}, and then
#'     downsampled (keep the even-indexed elements).}
#'   \item{d}{detail (difference) coefficients, obtained from convolving
#'     \code{x} with the wavelet (high-pass) filter \code{hi}, and then
#'     downsampled (keep the even-indexed elements).}
#' }
#'
#' @examples
#' # get Coiflet 30 coefficients
#' wv <- wfilters('c30')
#' lo <- rev(wv$lo)
#' hi <- rev(wv$hi)
#'
#' # general time-varying signal
#' time <- 1
#' fs <- 1000
#' x <- seq(0,time, length.out=time*fs)
#' y <- c(cos(2*pi*100*x)[1:300], cos(2*pi*50*x)[1:300],
#'        cos(2*pi*25*x)[1:200], cos(2*pi*10*x)[1:200])
#' op <- par(mfrow = c(3,1))
#' plot(x, y, type = "l", xlab = "Time", ylab = "Amplitude",
#'      main = "Original signal")
#' wt <- dwt(y, wname = NULL, lo, hi)
#'
#' x2 <- seq(1, length(x) - length(hi) + 1, 2)
#' plot(x2, wt$a, type = "h", xlab = "Time", ylab = "",
#'     main = "Approximation coefficients")
#' plot(x2, wt$d, type = "h", xlab = "Time", ylab = "",
#'      main = "Detail coefficients")
#' par (op)
#'
#' @author Lukas F. Reichlin.\cr
#' Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#'
#' @references [1] \url{https://CRAN.R-project.org/package=wavelets}
#' @references [2] \url{https://en.wikipedia.org/wiki/Daubechies_wavelet}
#' @references [3] \url{https://en.wikipedia.org/wiki/Coiflet}
#' @references [4]
#'   \url{https://en.wikipedia.org/wiki/Discrete_wavelet_transform}
#'
#' @rdname dwt
#' @export

 dwt <- function(x, wname = "d8", lo = NULL, hi = NULL) {

  # check parameters
  if (!is.vector(x) || !is.numeric(x)) {
    stop("x must be a numeric vector")
  }

  if (is.null(wname)) {
    if (is.null(lo) || is.null(hi)) {
      stop("both lo and hi needed when wname not specified")
    } else {
      if (length(lo) != length(hi)) {
        stop("lo and hi must be of equal length")
      }
    }
  } else {
    cf <- wfilters(wname)
    lo <- cf$lo
    hi <- cf$hi
  }

  tmp <- wconv("1d", x, lo, "valid")
  u <- tmp[seq(1, length(tmp), 2)]

  tmp <- wconv("1d", x, hi, "valid")
  v <- tmp[seq(1, length(tmp), 2)]

  list(a = u, d = v)
 }


 #' @rdname dwt
 #' @export

 wfilters <- function(wname) {

  lo <- switch(wname,
              "haar"   = c(1, 1),
              "d2"     = c(1, 1),
              "d4"     = c(0.6830127, 1.1830127, 0.3169873, -0.1830127),
              "d6"     = c(0.47046721, 1.14111692, 0.650365,
                           -0.19093442, -0.12083221, 0.0498175),
              "d8"     = c(0.32580343, 1.01094572, 0.89220014, -0.03957503,
                           -0.26450717, 0.0436163, 0.0465036, -0.01498699),
              "d10"    = c(0.22641898, 0.85394354, 1.02432694, 0.19576696,
                           -0.34265671, -0.04560113, 0.10970265, -0.00882680,
                           -0.01779187, 4.71742793e-3),
              "d12"    = c(0.15774243, 0.69950381, 1.06226376,
                           0.44583132, -0.31998660, -0.18351806,
                           0.13788809, 0.03892321, -0.04466375,
                           7.83251152e-4, 6.75606236e-3, -1.52353381e-3),
              "d14"    = c(0.11009943, 0.56079128, 1.03114849, 0.66437248,
                           -0.20351382, -0.31683501, 0.1008467, 0.11400345,
                           -0.05378245, -0.02343994, 0.01774979, 6.07514995e-4,
                           -2.54790472e-3, 5.00226853e-4),
              "d16"    = c(0.07695562, 0.44246725, 0.95548615,
                           0.82781653, -0.02238574, -0.40165863,
                           6.68194092e-4, 0.18207636, -0.02456390,
                           -0.06235021, 0.01977216, 0.01236884,
                           -6.88771926e-3, -5.54004549e-4, 9.55229711e-4,
                           -1.66137261e-4),
              "d18"    = c(0.05385035, 0.34483430, 0.85534906, 0.92954571,
                           0.18836955, -0.41475176, -0.13695355, 0.21006834,
                           0.043452675, -0.09564726, 3.54892813e-4, 0.03162417,
                           -6.67962023e-3, -6.05496058e-3, 2.61296728e-3,
                           3.25814671e-4, -3.56329759e-4, 5.5645514e-5),
              "d20"    = c(0.03771716, 0.26612218, 0.74557507, 0.97362811,
                           0.39763774, -0.35333620, -0.27710988, 0.18012745,
                           0.13160299, -0.10096657, -0.04165925, 0.04696981,
                           5.10043697e-3, -0.01517900, 1.97332536e-3,
                           2.81768659e-3, -9.69947840e-4, -1.64709006e-4,
                           1.32354367e-4, -1.875841e-5),
              "c6"     = c(-0.1028594569415370, 0.4778594569415370,
                           1.2057189138830700, 0.5442810861169260,
                           -0.1028594569415370, -0.0221405430584631),
              "c12"    = c(0.0231751934774337, -0.0586402759669371,
                           -0.0952791806220162, 0.5460420930695330,
                           1.1493647877137300, 0.5897343873912380,
                           -0.1081712141834230, -0.0840529609215432,
                           0.0334888203265590, 0.0079357672259240,
                           -0.0025784067122813, -0.0010190107982153),
              "c18"    = c(-0.0053648373418441, 0.0110062534156628,
                           0.0331671209583407, -0.0930155289574539,
                           -0.0864415271204239, 0.5730066705472950,
                           1.1225705137406600, 0.6059671435456480,
                           -0.1015402815097780, -0.1163925015231710,
                           0.0488681886423339, 0.0224584819240757,
                           -0.0127392020220977, -0.0036409178311325,
                           0.0015804102019152, 0.0006593303475864,
                           -0.0001003855491065, -0.0000489314685106),
              "c24"    = c(0.0012619224228619, -0.0023044502875399,
                           -0.0103890503269406, 0.0227249229665297,
                           0.0377344771391261, -0.1149284838038540,
                           -0.0793053059248983, 0.5873348100322010,
                           1.1062529100791000, 0.6143146193357710,
                           -0.0942254750477914, -0.1360762293560410,
                           0.0556272739169390, 0.0354716628454062,
                           -0.0215126323101745, -0.0080020216899011,
                           0.0053053298270610, 0.0017911878553906,
                           -0.0008330003901883, -0.0003676592334273,
                           0.0000881604532320, 0.0000441656938246,
                           -0.0000046098383254, -0.0000025243583600),
              "c30"    = c(-0.0002999290456692, 0.0005071055047161,
                           0.0030805734519904, -0.0058821563280714,
                           -0.0143282246988201, 0.0331043666129858,
                           0.0398380343959686, -0.1299967565094460,
                           -0.0736051069489375, 0.5961918029174380,
                           1.0950165427080700, 0.6194005181568410,
                           -0.0877346296564723, -0.1492888402656790,
                           0.0583893855505615, 0.0462091445541337,
                           -0.0279425853727641, -0.0129534995030117,
                           0.0095622335982613, 0.0034387669687710,
                           -0.0023498958688271, -0.0009016444801393,
                           0.0004268915950172, 0.0001984938227975,
                           -0.0000582936877724, -0.0000300806359640,
                           0.0000052336193200, 0.0000029150058427,
                           -0.0000002296399300, -0.0000001358212135)
              )
  lo <- lo / sqrt(2)
  ll <- length(lo)
  hi <- lo[ll:1] * (-1) ^ ((1:ll) - 1)

  list(lo = lo, hi = hi)

 }
