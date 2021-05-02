# Ma.R
# Copyright (C) 2006 EPRI Solutions, Inc.
# by Tom Short, tshort@eprisolutions.com
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
# 20200127  GvB       setup for gsignal v0.1.0
#------------------------------------------------------------------------------

#' Moving average (MA) model
#'
#' Create an MA model representing a filter or system model
#'
#' @param b moving average (MA) polynomial coefficients.
#'
#' @return A list of class \code{Ma} with the polynomial coefficients
#'
#' @seealso See also \code{\link{Arma}}
#'
#' @examples
#' f <- Ma(b = c(1, 2, 1) / 3)
#' freqz(f)
#' zplane(f)
#'
#' @author Tom Short, \email{tshort@@eprisolutions.com}
#' @export

Ma <- function(b) {
  class(b) <- "Ma"
  b
}
