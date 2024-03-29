% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/square.R
\name{square}
\alias{square}
\title{Square wave}
\usage{
square(t, duty = 50)
}
\arguments{
\item{t}{Time array, specified as a vector.}

\item{duty}{Duty cycle, specified as a real scalar from 0 to 100. Default:
50.}
}
\value{
Square wave, returned as a vector.
}
\description{
Generate a square wave of period \eqn{2\pi} with limits +1 and -1.
}
\details{
\code{y <- square(t)} generates a square wave with period \eqn{2\pi} for the
elements of the time array \code{t}.
\code{square} is similar to the sine function but creates a square wave with
values of –1 and 1.

\code{y <- square(t, duty)} generates a square wave with specified duty cycle
\code{duty}. The duty cycle is the percent of the signal period in which the
square wave is positive.
\if{latex}{
  \deqn{duty cycle = \frac{ontime * 100}{ontime + offtime}}
}
\if{html}{\preformatted{
                     ontime * 100
     duty cycle =  ----------------
                   ontime + offtime
}}
}
\examples{

## Create a vector of 100 equally spaced numbers from 0 to 3pi.
## Generate a square wave with a period of 2pi.
t <- seq(0, 3*pi, length.out = 100)
y <- square(t)
plot(t/pi, y, type="l", xlab = expression(t/pi), ylab = "")
lines (t/pi, sin(t), col = "red")

## Generate a 30 Hz square wave sampled at 1 kHz for 70 ms.
## Specify a duty cycle of 37\%.
## Add white Gaussian noise with a variance of 1/100.
t <- seq(0, 0.07, 1/1e3)
y <- square(2 * pi * 30 * t, 37) + rnorm(length(t)) / 10
plot(t, y, type="l", xlab = "", ylab = "")

}
\author{
Paul Kienzle.\cr
Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@gmail.com}.
}
