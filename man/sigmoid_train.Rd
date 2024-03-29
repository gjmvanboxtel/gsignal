% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sigmoid_train.R
\name{sigmoid_train}
\alias{sigmoid_train}
\title{Sigmoid Train}
\usage{
sigmoid_train(t, ranges, rc)
}
\arguments{
\item{t}{Vector (or coerced to a vector) of time values at which the sigmoids
are calculated.}

\item{ranges}{Matrix or array with 2 columns containing the time values
within \code{t} at which each sigmoid is evaluated. The number of sigmoids
is determined by the number of rows in \code{ranges}.}

\item{rc}{Time constant. Either a scalar or a matrix or array with 2 columns
containing the rising and falling time constants of each sigmoid. If a
matrix or array is passed in \code{rc}, its size must equal the size of
\code{ranges}. If a single scalar is passed in \code{rc}, then all sigmoids
have the same time constant and are symmetrical.}
}
\value{
A list consisting two variables; \code{y} the combined sigmoid train
  (length identical to \code{t}), and \code{s}, the individual sigmoids
  (number of rows equal to number of rows in \code{ranges} and \code{rc}.
}
\description{
Evaluate a train of sigmoid functions at \code{t}.
}
\details{
The number and duration of each sigmoid is determined from ranges. Each row
of \code{ranges} represents a real interval, e.g. if sigmoid \code{i} starts
at \code{t = 0.1} and ends at \code{t = 0.5}, then \code{ranges[i, ] = c(0.1,
0.5)}. The input \code{rc} is an array that defines the rising and falling
time constants of each sigmoid. Its size must equal the size of ranges.

The individual sigmoids are returned in \code{s}. The combined sigmoid train
is returned in the vector \code{y} of length equal to \code{t}, and such that
\code{y = max(s)}.
}
\examples{

t <- seq(0, 2, length.out = 500)
ranges <- rbind(c(0.1, 0.4), c(0.6, 0.8), c(1, 2))
rc <- rbind(c(1e-2, 1e-3), c(1e-3, 2e-2), c(2e-2, 1e-2))
st <- sigmoid_train (t, ranges, rc)
plot(t, st$y[1,], type="n", xlab = "Time(s)", ylab = "S(t)",
     main = "Vectorized use of sigmoid train")
for (i in 1:3) rect(ranges[i, 1], 0, ranges[i, 2], 1,
                    border = NA, col="pink")
for (i in 1:3) lines(t, st$y[i,])
# The colored regions show the limits defined in range.

t <- seq(0, 2, length.out = 500)
ranges <- rbind(c(0.1, 0.4), c(0.6, 0.8), c(1, 2))
rc <- rbind(c(1e-2, 1e-3), c(1e-3, 2e-2), c(2e-2, 1e-2))
amp <- c(4, 2, 3)
st <- sigmoid_train (t, ranges, rc)
y <- amp \%*\% st$y
plot(t, y[1,], type="l", xlab = 'time', ylab = 'signal',
     main = 'Varying amplitude sigmoid train', col="blue")
lines(t, st$s, col = "orange")
legend("topright", legend = c("Sigmoid train", "Components"),
       lty = 1, col = c("blue", "orange"))

}
\author{
Juan Pablo Carbajal, \email{carbajal@ifi.uzh.ch}.\cr
Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@gmail.com}.
}
