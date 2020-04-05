# Sos.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
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
# 20200402  GvB       setup for gsignal v0.1.0
#---------------------------------------------------------------------------------------------------------------------

#' Second-order sections
#' 
#' Create or convert filter models to second-order sections form
#' 
#' \code{as.Sos} converts from other forms, including \code{Arma}, \code{Ma},
#' and \code{Zpg}.
#' 
#' @param sos second-order sections representation of the model
#' @param g overall gain factor
#' @param x model to be converted.
#' @param ...	additional arguments (ignored).
#' 
#' @return A list of class \code{'Sos'} with the following list elements:
#' \describe{
#'   \item{sos}{second-order section representation of the model, returned as an
#'     \code{L x 6} matrix, one row for each section \code{1:L}. Each row
#'     constists of an \code{[B, A]}, pair, where \code{B = c(b0, b1, b2)}, and
#'     \code{A = c(1, a1, a2)}, the filer coeffients for each section. Each
#'     \code{b0} entry must be nonzero for each section.}
#'   \item{g}{overall gain factor that scales any one of the \eqn{B_i} vectors}
#' }
#' 
#' @seealso \code{\link{Arma}}, \code{\link{Ma}}, \code{\link{Zpg}} 
#' 
#' @author Geert van Boxtel \email{G.J.M.vanBoxtel@@gmail.com}.
#' 
#' @rdname Sos
#' @export

Sos <- function(sos, g) {
  res <- list(sos = sos, g = g)
  class(res) <- "Sos"
  res
}

#' @rdname Sos
#' @export
as.Sos <- function(x, ...) UseMethod("as.Sos")

#' @rdname Sos
#' @usage
#' ## S3 method for class 'Arma'
#' as.Sos(x, ...)
#' @export
as.Sos.Arma <- function(x, ...) {
  sosg <- tf2sos(x$b, x$a)
  Sos(sosg$sos, sosg$g)
}

#' @rdname Sos
#' @usage
#' ## S3 method for class 'Ma'
#' as.Sos(x, ...)
#' @export
as.Sos.Ma <- function(x, ...) {
  as.Sos.Arma(as.Arma.Ma(x))
}

#' @rdname Sos
#' @usage
#' ## S3 method for class 'Sos'
#' as.Sos(x, ...)
#' @export
as.Sos.Sos <- function(x, ...) x

#' @rdname Sos
#' @usage
#' ## S3 method for class 'Zpg'
#' as.Sos(x, ...)
#' @export
as.Sos.Zpg <- function(x, ...) {
  ret <- zp2sos(x$z, x$p, x$k)
  class(ret) <- "Sos"
  ret
}