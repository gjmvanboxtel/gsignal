% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sos2zp.R
\name{sos2zp}
\alias{sos2zp}
\title{Sos to zero-pole-gain}
\usage{
sos2zp(sos, g = 1)
}
\arguments{
\item{sos}{Second-order section representation, specified as an nrow-by-6
matrix, whose rows contain the numerator and denominator coefficients of
the second-order sections:\cr \code{sos <- rbind(cbind(B1, A1), cbind(...),
cbind(Bn, An))}, where \code{B1 <- c(b0, b1, b2)}, and \code{A1 <- c(a0,
a1, a2)} for section 1, etc. The b0 entry must be nonzero for each section.}

\item{g}{Overall gain factor that effectively scales the output \code{b}
vector (or any one of the input \code{B_i} vectors). Default: 1.}
}
\value{
A list of class "Zpg" with the following list elements:
\describe{
  \item{z}{complex vector of the zeros of the model (roots of \code{B(z)})}
  \item{p}{complex vector of the poles of the model (roots of \code{A(z)})}
  \item{k}{overall gain (\code{B(Inf)})}
}
}
\description{
Convert digital filter second-order section data to zero-pole-gain form.
}
\examples{
sos <- rbind(c(1, 0, 1, 1, 0, -0.81), c(1, 0, 0, 1, 0, 0.49))
zpk <- sos2zp(sos)

}
\seealso{
\code{\link{filter}}
}
\author{
Julius O. Smith III \email{jos@ccrma.stanford.edu}.\cr
Conversion to R by, Geert van Boxtel \email{G.J.M.vanBoxtel@gmail.com}
}
