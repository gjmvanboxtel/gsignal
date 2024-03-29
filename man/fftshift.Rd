% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fftshift.R
\name{fftshift}
\alias{fftshift}
\title{Zero-frequency shift}
\usage{
fftshift(x, MARGIN = 2)
}
\arguments{
\item{x}{input data, specified as a vector or matrix.}

\item{MARGIN}{dimension to operate along, 1 = row, 2 = columns (default).
Specifying \code{MARGIN = c(1, 2)} centers along both rows and columns.
Ignored when \code{x} is a vector.}
}
\value{
vector or matrix with centered frequency.
}
\description{
Perform a shift in order to move the frequency 0 to the center of the input.
}
\details{
If \code{x} is a vector of \code{N} elements corresponding to \code{N} time
samples spaced by \code{dt}, then \code{fftshift(x)} corresponds to
frequencies \code{f = c(-seq(ceiling((N-1)/2), 1, -1), 0, (1:floor((N-1)/2)))
* df}, where \code{df = 1 / (N * dt)}. In other words, the left and right
halves of \code{x} are swapped.

If \code{x} is a matrix, then \code{fftshift} operates on the rows or columns
of \code{x}, according to the \code{MARGIN} argument, i.e. it swaps the the
upper and lower halves of the matrix \code{(MARGIN = 1)}, or the left and
right halves of the matrix \code{(MARGIN = 2)}. Specifying \code{MARGIN =
c(1, 2)} swaps along both dimensions, i.e., swaps the first quadrant with the
fourth, and the second with the third.
}
\examples{
Xeven <- 1:6
ev <- fftshift(Xeven)   # returns 4 5 6 1 2 3

Xodd <- 1:7
odd <- fftshift(Xodd)   # returns 5 6 7 1 2 3 4

fs <- 100                      # sampling frequency
t <- seq(0, 10 - 1/fs, 1/fs)   # time vector
S <- cos(2 * pi * 15 * t)
n <- length(S)
X <- fft(S)
f <- (0:(n - 1)) * (fs / n);   # frequency range
power <- abs(X)^2 / n          # power
plot(f, power, type="l")
Y <- fftshift(X)
fsh <- ((-n/2):(n/2-1)) * (fs / n)  # zero-centered frequency range
powersh <- abs(Y)^2 / n             # zero-centered power
plot(fsh, powersh, type = "l")

}
\seealso{
\code{ifftshift}
}
\author{
Vincent Cautaerts, \email{vincent@comf5.comm.eng.osaka-u.ac.jp},\cr
  adapted by John W. Eaton.\cr
  Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@gmail.com}.
}
