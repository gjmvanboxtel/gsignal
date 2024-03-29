% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dct2.R
\name{dct2}
\alias{dct2}
\title{2-D Discrete Cosine Transform}
\usage{
dct2(x, m = NROW(x), n = NCOL(x))
}
\arguments{
\item{x}{2-D numeric matrix}

\item{m}{Number of rows, specified as a positive integer. \code{dct2} pads or
truncates \code{x} so that is has \code{m} rows. Default: \code{NROW(x)}.}

\item{n}{Number of columns, specified as a positive integer. \code{dct2} pads
or truncates \code{x} so that is has \code{n} columns. Default:
\code{NCOL(x)}.}
}
\value{
\code{m}-by-\code{n} numeric discrete cosine transformed matrix.
}
\description{
Compute the two-dimensional discrete cosine transform of a matrix.
}
\details{
The discrete cosine transform (DCT) is closely related to the discrete
Fourier transform. It is a separable linear transformation; that is, the
two-dimensional transform is equivalent to a one-dimensional DCT performed
along a single dimension followed by a one-dimensional DCT in the other
dimension.
}
\examples{
A <- matrix(runif(100), 10, 10)
B <- dct2(A)

}
\seealso{
\code{\link{idct2}}
}
\author{
Paul Kienzle, \email{pkienzle@users.sf.net}.\cr Conversion to R by
  Geert van Boxtel, \email{G.J.M.vanBoxtel@gmail.com}.
}
