# signals.R
# Copyright (C) 2021  Geert van Boxtel, <G.J.M.vanBoxtel@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
# Version history
# 20210414  GvB       Initial setup (v0.3-1)
#------------------------------------------------------------------------------

#' signals
#'
#' Sample EEG and ECG data.
#'
#' @author Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}
#'
#' @examples
#' data(signals)
#' time <- seq(0, 10, length.out = nrow(signals))
#' op <- par(mfcol = c(2, 1))
#' plot(time, signals[, 1], type = "l", xlab = "Time", ylab = "EEG (uV)")
#' plot(time, signals[, 2], type = "l", xlab = "Time", ylab = "ECG (uV)")
#' par(op)
#'
#' @format A \code{\link{data.frame}} containing 10 seconds of data
#'   electrophysiological data, sampled at 256 Hz with a 24 bit A/D converter,
#'   measured in microVolts. The data frame consists of 2 columns (channels):
#' \describe{
#'   \item{eeg}{electroencephalogram (EEG) data measured from electrode Pz
#'   according to the 10-20 system, referred to algebraically linked mastoids
#'   (the brain's alpha rhythm is clearly visible).}
#'   \item{ecg}{electrocardiogram (ECG) data, recorded bipolarly with a V6
#'   versus V1 chest lead (this lead maximizes the R wave of the ECG with
#'   respect to the P, Q, S, T and U waves of the cardiac cycle).}
#' }
#'
"signals"
