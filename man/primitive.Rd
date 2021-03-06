% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/primitive.R
\name{primitive}
\alias{primitive}
\title{Primitive}
\usage{
primitive(FUN, t, C = 0)
}
\arguments{
\item{FUN}{the function to calculate the primitive of.}

\item{t}{points at which the function \code{FUN} is evaluated, specified as a
vector of ascending values}

\item{C}{constant of integration. Default: 0}
}
\value{
Vector of integrated function values.
}
\description{
Calculate the indefinitive integral of a function.
}
\details{
This function is a fancy way of calculating the cumulative sum.
}
\examples{
f <- function(t) sin(2 * pi * 3 * t)
t <- c(0, sort(runif(100)))
F <- primitive (f, t, 0)
t_true <- seq(0, 1, length.out = 1e3)
F_true <- (1 - cos(2 * pi * 3 * t_true)) / (2 * pi * 3)
plot (t, F, xlab = "", ylab = "")
lines (t_true, F_true, col = "red")
legend("topright", legend = c("Numerical primitive", "True primitive"),
  lty = c(0, 1), pch = c(1, NA), col = 1:2)

}
\seealso{
\code{\link{cumsum}}
}
\author{
Juan Pablo Carbajal.\cr
Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@gmail.com}
}
