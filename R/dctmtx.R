# dctmtx.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Original Octave code:
# Copyright (C) 2001 Paul Kienzle <pkienzle@users.sf.net>
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
# See also: http://www.gnu.org/licenses/gpl-2.0.txt
#
# Version history
# 20201016  GvB       setup for gsignal v0.1.0
#---------------------------------------------------------------------------------------------------------------------

#' Discrete Cosine Transform Matrix
#' 
#' Compute the discrete cosine transform matrix
#' 
#' A DCT transformation matrix is useful for doing things like JPEG image
#' compression, in which an 8x8 DCT matrix is applied to non-overlapping blocks
#' throughout an image and only a subblock on the top left of each block is
#' kept.  During restoration, the remainder of the block is filled with zeros
#' and the inverse transform is applied to the block.
#' 
#' The two-dimensional DCT of A can be computed as \code{D * A * t(D)}. This
#' computation is sometimes faster than using \code{dct2}, especially if you are
#' computing a large number of small DCTs, because D needs to be determined only
#' once. For example, in JPEG compression, the DCT of each 8-by-8 block is
#' computed. To perform this computation, use \code{dctmtx} to determine D of
#' imput image A, and then calculate each DCT using \eqn{D * A * t(D)} (where A
#' is each 8-by-8 block). This is faster than calling \code{dct2} for each
#' individual block. See the examples.
#'
#' @param n Size of DCT matrix, specified as a positive integer. Default: none.
#'  
#' @return Discrete cosine transform, returned as a vector or matrix.
#'  
#' @examples
#' D <- dctmtx(8)
#' 
#' \dontrun{
#' library(imager)
#' op <- par(mfrow = c(1, 2))
#' A <- grayscale(boats)[1:256, 65:320, 1, 1]  #256*256 just to show example
#' plot(as.cimg(im), axes = FALSE, main = "Original")
#' D <- dctmtx(nrow(A))
#' dct <- D * A * t(D)
#' plot(as.cimg(dct), axes = FALSE, main = "DCT2 of Original")
#' par(op)
#' }
#'         
#' @author Paul Kienzle \email{pkienzle@@users.sf.net}; port to R by Geert van
#'   Boxtel \email{G.J.M.vanBoxtel@@gmail.com}.
#' 
#' @seealso \code{\link{dct}}, \code{\link{dct2}}, \code{\link{idct}}, \code{\link{idct2}}
#' 
#' @export

 dctmtx <- function (n) {

  # check parameters
  if (!isPosscal(n) || !isWhole(n)) {
    stop('n must be a positive integer')
  }
  
  if (n > 1) {
    T <- rbind(sqrt(1 / n) * rep(1, n),
               sqrt(2 / n) * cos((pi / 2 / n) * seq(1, (n - 1)) %o% seq(1, 2 * n, 2)))
  } else if (n == 1) {
    T <- 1
  } else {
    stop("n must be >= 1")
  }
  T
}