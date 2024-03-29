% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Zpg.R
\name{Zpg}
\alias{Zpg}
\alias{as.Zpg}
\alias{as.Zpg.Arma}
\alias{as.Zpg.Ma}
\alias{as.Zpg.Sos}
\alias{as.Zpg.Zpg}
\title{Zero pole gain model}
\usage{
Zpg(z, p, g)

as.Zpg(x, ...)

## S3 method for class 'Arma'
as.Zpg(x, ...)

## S3 method for class 'Ma'
as.Zpg(x, ...)

## S3 method for class 'Sos'
as.Zpg(x, ...)

## S3 method for class 'Zpg'
as.Zpg(x, ...)
}
\arguments{
\item{z}{complex vector of the zeros of the model.}

\item{p}{complex vector of the poles of the model.}

\item{g}{overall gain of the model.}

\item{x}{model to be converted.}

\item{...}{additional arguments (ignored).}
}
\value{
A list of class Zpg with the following list elements:
\describe{
  \item{z}{complex vector of the zeros of the model}
  \item{p}{complex vector of the poles of the model}
  \item{g}{gain of the model}
}
}
\description{
Create an zero pole gain model of an ARMA filter, or convert other forms to a
Zpg model.
}
\details{
\code{as.Zpg} converts from other forms, including \code{Arma} and \code{Ma}.
}
\examples{
## design notch filter at pi/4 radians = 0.5/4 = 0.125 * fs
w = pi/4
# 2 poles, 2 zeros
# zeroes at r = 1
r <- 1
z1 <- r * exp(1i * w)
z2 <- r * exp(1i * -w)
# poles at r = 0.9
r = 0.9
p1 <- r * exp(1i * w)
p2 <- r * exp(1i * -w)

zpg <- Zpg(c(z1, z2), c(p1, p2), 1)
zplane(zpg)
freqz(zpg)

## Sharper edges: increase distance between zeros and poles
r = 0.8
p1 <- r * exp(1i * w)
p2 <- r * exp(1i * -w)
zpg <- Zpg(c(z1, z2), c(p1, p2), 1)
zplane(zpg)
freqz(zpg)

}
\seealso{
See also \code{\link{Arma}}
}
\author{
Tom Short, \email{tshort@eprisolutions.com},\cr
 adapted by Geert van Boxtel, \email{gjmvanboxtel@gmail.com}.
}
