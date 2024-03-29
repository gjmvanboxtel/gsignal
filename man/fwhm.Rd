% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fwhm.R
\name{fwhm}
\alias{fwhm}
\title{Full width at half maximum}
\usage{
fwhm(
  x = seq_len(length(y)),
  y,
  ref = c("max", "zero", "middle", "min", "absolute"),
  level = 0.5
)
}
\arguments{
\item{x}{samples at which \code{y} is measured, specified as a vector. I.e.,
\code{y} is sampled as \code{y[x]}. Default: \code{seq_len(length(y))}.}

\item{y}{signal to find the width of. If \code{y} is a matrix, widths of all
columns are computed.}

\item{ref}{reference. Compute the width with reference to:
\describe{
  \item{\code{"max" | "zero"}}{\code{max(y)}}
  \item{\code{"middle" | "min"}}{\code{min(y) + max(y)}}
  \item{\code{"absolute"}}{an absolute level of \code{y}}
}}

\item{level}{the level at which to compute the width. Default: 0.5.}
}
\value{
Full width at half maximum, returned as a vector with a length equal
  to the number of columns in \code{y}, or 1 in case of a vector.
}
\description{
Compute peak full-width at half maximum or at another level of peak maximum
for a vector or matrix.
}
\examples{
x <- seq(-pi, pi, 0.001)
y <- cos(x)
w <- fwhm(x, y)
m <- x[which.max(y)]
f <- m - w/2
t <- m + w/2
plot(x, y, type="l",
     panel.first = {
       usr <- par('usr')
       rect(f, usr[3], t, usr[4], col = rgb(0, 1, 0, 0.4), border = NA)
     })
abline(h = max(y) / 2, lty = 2, col = "gray")

}
\author{
Petr Mikulik.\cr
 Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@gmail.com}.
}
