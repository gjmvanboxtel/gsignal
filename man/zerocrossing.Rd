% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zerocrossing.R
\name{zerocrossing}
\alias{zerocrossing}
\title{Zero Crossing}
\usage{
zerocrossing(x, y)
}
\arguments{
\item{x}{the x-coordinates of points in the function.}

\item{y}{the y-coordinates of points in the function.}
}
\value{
Zero-crossing points
}
\description{
Estimate zero crossing points of waveform.
}
\examples{
x <- seq(0, 1, length.out = 100)
y <- runif(100) - 0.5
x0 <- zerocrossing(x, y)
plot(x, y, type ="l", xlab = "", ylab = "")
points(x0, rep(0, length(x0)), col = "red")

}
\author{
Carlo de Falco, \email{carlo.defalco@gmail.com}.\cr
Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@gmail.com}.
}
