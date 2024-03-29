% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mexihat.R
\name{mexihat}
\alias{mexihat}
\title{Mexicat Hat}
\usage{
mexihat(lb = -5, ub = 5, n = 1000)
}
\arguments{
\item{lb, ub}{Lower and upper bounds of the interval to evaluate the wavelet
on. Default: -5 to 5.}

\item{n}{Number of points on the grid between \code{lb} and \code{ub} (length
of the wavelet). Default: 1000.}
}
\value{
A list containing 2 variables; \code{x}, the grid on which the
  complex Mexican Hat wavelet was evaluated, and \code{psi} (\eqn{\Psi}), the
  evaluated wavelet on the grid \code{x}.
}
\description{
Generate a Mexican Hat (Ricker) wavelet sampled on a regular grid.
}
\details{
The Mexican Hat or Ricker wavelet is the negative normalized second
derivative of a Gaussian function, i.e., up to scale and normalization, the
second Hermite function. It is a special case of the family of continuous
wavelets (wavelets used in a continuous wavelet transform) known as Hermitian
wavelets. The Ricker wavelet is frequently employed to model seismic data,
and as a broad spectrum source term in computational electrodynamics. It is
usually only referred to as the Mexican hat wavelet in the Americas, due to
taking the shape of a sombrero when used as a 2D image processing kernel. It
is also known as the Marr wavelet (source: Wikipedia)
}
\examples{

mh <- mexihat(-5, 5, 1000)
plot(mh$x, mh$psi, type="l", main = "Mexican Hat Wavelet",
     xlab = "", ylab = "")

}
\author{
Sylvain Pelissier, \email{sylvain.pelissier@gmail.com}.\cr
Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@gmail.com}.
}
