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
#
#---------------------------------------------------------------------------------------------------------------------

#' Internal functions not exported to the namespace
#' 
#' @keywords internal
#' @noRd

# test if x is a scalar
isScalar <- function(x) is.atomic(x) && length(x) == 1L

# test if x is a positive scalar    
isPosscal <- function (x) isScalar(x) && is.numeric(x) && x >= 0

# test if x is a whole number
isWhole <- function(x, tol = .Machine$double.eps^0.5)  !(is.null(x) || is.character(x)) && abs(x - round(x)) < tol

# convert factor to numeric
unfactor <- function (f) if (is.factor(f)) as.numeric(levels(f)[as.integer(f)]) else NULL

# # compute next power of 2 
# nextpow <- function (x) 2^ceiling(log2(x))