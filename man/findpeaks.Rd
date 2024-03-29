% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/findpeaks.R
\name{findpeaks}
\alias{findpeaks}
\title{Find local extrema}
\usage{
findpeaks(
  data,
  MinPeakHeight = .Machine$double.eps,
  MinPeakDistance = 1,
  MinPeakWidth = 1,
  MaxPeakWidth = Inf,
  DoubleSided = FALSE
)
}
\arguments{
\item{data}{the data, expected to be a vector or one-dimensional array.}

\item{MinPeakHeight}{Minimum peak height (non-negative scalar). Only peaks
that exceed this value will be returned. For data taking positive and
negative values use the option \code{DoubleSided}. Default:
\code{.Machine$double.eps}.}

\item{MinPeakDistance}{Minimum separation between peaks (positive integer).
Peaks separated by less than this distance are considered a single peak.
This distance is also used to fit a second order polynomial to the peaks to
estimate their width, therefore it acts as a smoothing parameter. The
neighborhood size is equal to the value of \code{MinPeakDistance}. Default:
1.}

\item{MinPeakWidth}{Minimum width of peaks (positive integer). The width of
the peaks is estimated using a parabola fitted to the neighborhood of each
peak. The width is calculated with the formula \eqn{a * (width - x0)^{2} =
1}, where a is the the concavity of the parabola and x0 its vertex.
Default: 1.}

\item{MaxPeakWidth}{Maximum width of peaks (positive integer). Default:
\code{Inf}.}

\item{DoubleSided}{Tells the function that data takes positive and negative
values. The baseline for the peaks is taken as the mean value of the
function. This is equivalent as passing the absolute value of the data
after removing the mean. Default: FALSE}
}
\value{
A list containing the following elements:
\describe{
  \item{pks}{The value of data at the peaks.}
  \item{loc}{The index indicating the position of the peaks.}
  \item{parabol}{A list containing the parabola fitted to each returned peak.
  The list has two fields, \code{x} and \code{pp}. The field \code{pp}
  contains the coefficients of the 2nd degree polynomial and \code{x} the
  extrema of the interval where it was fitted.}
  \item{height}{The estimated height of the returned peaks (in units of
  data).}
  \item{baseline}{The height at which the roots of the returned peaks were
  calculated (in units of data).}
  \item{roots}{The abscissa values (in index units) at which the parabola
  fitted to each of the returned peaks realizes its width as defined below.}
}
}
\description{
Return peak values and their locations of the vector \code{data}.
}
\details{
Peaks of a positive array of \code{data} are defined as local maxima. For
double-sided data, they are maxima of the positive part and minima of the
negative part. \code{data} is expected to be a one-dimensional vector.
}
\examples{
### demo 1
t <- 2 * pi * seq(0, 1,length = 1024)
y <- sin(3.14 * t) + 0.5 * cos(6.09 * t) +
     0.1 * sin(10.11 * t + 1 / 6) + 0.1 * sin(15.3 * t + 1 / 3)

data1 <- abs(y) # Positive values
peaks1 <- findpeaks(data1)

data2 <- y # Double-sided
peaks2 <- findpeaks(data2, DoubleSided = TRUE)
peaks3 <- findpeaks (data2, DoubleSided = TRUE, MinPeakHeight = 0.5)

op <- par(mfrow=c(1,2))
plot(t, data1, type="l", xlab="", ylab="")
points(t[peaks1$loc], peaks1$pks, col = "red", pch = 1)
plot(t, data2, type = "l", xlab = "", ylab = "")
points(t[peaks2$loc], peaks2$pks, col = "red", pch = 1)
points(t[peaks3$loc], peaks3$pks, col = "red", pch = 4)
legend ("topleft", "0: >2*sd, x: >0.5", bty = "n",
        text.col = "red")
par (op)
title("Finding the peaks of smooth data is not a big deal")

## demo 2
t <- 2 * pi * seq(0, 1, length = 1024)
y <- sin(3.14 * t) + 0.5 * cos(6.09 * t) + 0.1 *
     sin(10.11 * t + 1 / 6) + 0.1 * sin(15.3 * t + 1 / 3)
data <- abs(y + 0.1*rnorm(length(y),1))   # Positive values + noise
peaks1 <- findpeaks(data, MinPeakHeight=1)
dt <- t[2]-t[1]
peaks2 <- findpeaks(data, MinPeakHeight=1, MinPeakDistance=round(0.5/dt))
op <- par(mfrow=c(1,2))
plot(t, data, type="l", xlab="", ylab="")
points (t[peaks1$loc],peaks1$pks,col="red", pch=1)
plot(t, data, type="l", xlab="", ylab="")
points (t[peaks2$loc],peaks2$pks,col="red", pch=1)
par (op)
title(paste("Noisy data may need tuning of the parameters.\n",
            "In the 2nd example, MinPeakDistance is used\n",
            "as a smoother of the peaks"))

}
\author{
Juan Pablo Carbajal, \email{carbajal@ifi.uzh.ch}.\cr
Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@gmail.com}.
}
