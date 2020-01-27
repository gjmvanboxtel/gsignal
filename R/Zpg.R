# Zpg.R
# Copyright (C) 2006 EPRI Solutions, Inc.
# by Tom Short, tshort@eprisolutions.com
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
# 20200127  GvB       setup for gsignal v0.1.0
#---------------------------------------------------------------------------------------------------------------------

#' Zero pole gain model
#' 
#' Create an zero pole gain model of an ARMA filter
#' 
#' \code{as.Zpg} converts from other forms, including \code{Arma} and \code{Ma}.
#' 
#' @param zero complex vector of the zeros of the model.
#' @param pole complex vector of the poles of the model.
#' @param gain gain of the model.
#' @param x model to be converted.
#' @param ...	additional arguments (ignored).
#' 
#' @return A list of class Zpg with the following list elements:
#' \describe{
#'   \item{zero}{complex vector of the zeros of the model}
#'   \item{pole}{complex vector of the poles of the model}
#'   \item{gain}{gain of the model}
#' }
#' 
#' @seealso See also \code{\link{Arma}}
#' 
#' @examples
#' filt <- Zpg(c(-1, -1), -1, 1/3)
#' #zplane(filt)
#' 
#' @author Tom Short \email{tshort@@eprisolutions.com}
#' @rdname Zpg
#' @export

Zpg <- function(zero, pole, gain) {
  res <- list(zero = zero, pole = pole, gain = gain)
  class(res) <- "Zpg"
  res
}

#' @rdname Zpg
#' @export
as.Zpg <- function(x, ...) UseMethod("as.Zpg")

#' @rdname Zpg
#' @usage
#' ## S3 method for class 'Arma'
#' as.Zpg(x, ...)
#' @export
as.Zpg.Arma <- function(x, ...) {
  Zpg(pole = roots(x$a), zero = roots(x$b), gain = x$b[1] / x$a[1])
}

#' @rdname Zpg
#' @usage
#' ## S3 method for class 'Ma'
#' as.Zpg(x, ...)
#' @export
as.Zpg.Ma <- function(x, ...) {
  as.Zpg(as.Arma(x))
}

#' @rdname Zpg
#' @usage
#' ## S3 method for class 'Zpg'
#' as.Zpg(x, ...)
#' @export
as.Zpg.Zpg <- function(x, ...) x


