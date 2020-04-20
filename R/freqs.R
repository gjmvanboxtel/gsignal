# freqs.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Original Octave function:
# Copyright (C) 2003 Julius O. Smith III <jos@ccrma.stanford.edu>
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
# 20200413  GvB       setup for gsignal v0.1.0
#---------------------------------------------------------------------------------------------------------------------

#' Frequency response of analog filters
#' 
#' Compute the s-plane frequency response of an IIR filter
#' 
#' The s-plane frequency response of the IIR filter \code{B(s)/A(s)} is computed
#' as \code{H = polyval(B, 1i * W) / polyval(A, 1i * W)}. If called with no
#' output argument, a plot of magnitude and phase are displayed.
#' 
#' @param b moving average (MA) polynomial coefficients, specified as a numeric
#'   vector or matrix. In case of a matrix, then each row corresponds to an
#'   output of the system. The number of columns of \code{b} must be less than
#'   or equal to the length of \code{a}.
#' @param a autoregressive (AR) polynomial coefficients, specified as a vector.
#' @param w angular frequencies, specified as a positive real vector expressed
#'   in rad/second.
#' @param plot logical. If \code{TRUE} (default), plots the magnitude and phase
#'   responses as a function of angular frequency, otherwise a vector of the
#'   frequency response is returned.
#' 
#' @return Frequency response, returned as a complex vector.
#' 
#' @examples
#' b <- c(1, 2); a <- c(1, 1)
#' w <- seq(0, 4, length.out = 128)
#' freqs (b, a, w)
#' 
#' @author Julius O. Smith III \email{jos@@ccrma.stanford.edu}, port to R by
#'   Geert van Boxtel \email{gjmvanboxtel@@gmail.com}
#' 
#' @export

freqs <- function(b, a, w, plot = TRUE) {
  
  h <- pracma::polyval(b, 1i * w) / pracma::polyval(a, 1i * w)
  
  if (plot) {
    freqs_plot(w, h)
    invisible(h)
  } else {
    h
  }
}