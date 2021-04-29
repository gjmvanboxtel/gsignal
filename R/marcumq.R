# marcumq.R
# Copyright (C) 2020 William Asquith <william.asquith at ttu.edu>
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
# 20201123  GvB       setup for gsignal v0.1.0
#------------------------------------------------------------------------------

#' Marcum Q function
#'
#' Compute the generalized Marcum Q function
#'
#' The code for this function was taken from the help file of the \code{cdfkmu}
#' function in the \code{lmomco} package, based on a suggestion of Daniel
#' Wollschlaeger.
#'
#' @param a,b input arguments, specified as non-negative real numbers.
#' @param m order, specified as a positive integer
#'
#' @return Marcum Q function.
#'
#' @examples
#' mq <- marcumq(12.4, 12.5)
#'
#' @author William Asquith, \email{william.asquith@@ttu.edu}.
#'
#' @references \url{https://cran.r-project.org/package=lmomco}
#
#' @export

marcumq <- function(a, b, m = 1) {
  stats::pchisq(b^2, df = 2 * m, ncp = a^2, lower.tail = FALSE)
}
