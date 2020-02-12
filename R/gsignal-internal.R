# gsignal-internal.R - internal or barely commented functions not exported from the namespace
#
# Copyright (C) 2019  Geert van Boxtel,
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
# 20191029  GvB       Initial setup
# 20200112  GvB       Added ssq() msq() rmsq()
#
#---------------------------------------------------------------------------------------------------------------------

#' Internal functions not exported to the namespace
#' 
#' @keywords internal
#' @noRd

# test if x is a scalar
isScalar <- function(x) ifelse(is.character(x), nchar(x) == 1L, (is.atomic(x) && length(x) == 1L))

# test if x is a positive scalar    
isPosscal <- function(x) isScalar(x) && is.numeric(x) && x >= 0

# test if x is a whole number
isWhole <- function(x, tol = .Machine$double.eps^0.5)  !(is.null(x) || is.character(x)) && abs(x - round(x)) < tol

# convert factor to numeric
unfactor <- function(f) if (is.factor(f)) as.numeric(levels(f)[as.integer(f)]) else NULL

# sinc function
sinc <- function(x) ifelse(x == 0, 1, sin(x) / x)

# sum of squares (assume input is a vector)
ssq <- function(x) ifelse(is.complex(x), sum(Re(x * Conj(x))), sum(x * x))

# mean sum of squares (assume input is a vector)
msq <- function(x) ssq(x) / length(x)

# root mean square (assume input is a vector)
rmsq <- function(x) sqrt(msq(x))

# # compute next power of 2 
nextpow <- function (x) 2^ceiling(log2(x))
