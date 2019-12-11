% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/triang.R
\name{triang}
\alias{triang}
\title{Triangular window}
\usage{
triang(n)
}
\arguments{
\item{n}{Window length, specified as a positive integer.}
}
\value{
triangular window, returned as a vector.
}
\description{
Return the filter coefficients of a triangular window of length \code{n}
}
\details{
Unlike the Bartlett window, \code{triang} does not go to zero at the edges of the window.
For odd \code{n}, \code{triang(n)} is equal to \code{bartlett(m + 2)} except for the zeros
at the edges of the window.
 
If you specify a one-point window \code{(n = 1)}, the value 1 is returned.
}
\examples{

t <- triang(64)
plot (t, type = "l", xlab = "Samples", ylab =" Amplitude")

}
\seealso{
\code{\link{bartlett}}
}
\author{
Original Octave code Copyright (C) 1995-2017 Andreas Weingessel \email{Andreas.Weingessel@ci.tuwien.ac.at}.
Port to R by Geert van Boxtel \email{G.J.M.vanBoxtel@gmail.com}.
}