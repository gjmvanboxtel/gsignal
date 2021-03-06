% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ellip.R
\name{ellip}
\alias{ellip}
\alias{ellip.FilterSpecs}
\alias{ellip.default}
\title{Elliptic filter design}
\usage{
ellip(n, ...)

\method{ellip}{FilterSpecs}(n, Rp = n$Rp, Rs = n$Rs, w = n$Wc, type = n$type, plane = n$plane, ...)

\method{ellip}{default}(
  n,
  Rp,
  Rs,
  w,
  type = c("low", "high", "stop", "pass"),
  plane = c("z", "s"),
  ...
)
}
\arguments{
\item{n}{filter order.}

\item{...}{additional arguments passed to ellip, overriding those given by n
of class \code{FilterSpecs}.}

\item{Rp}{dB of passband ripple.}

\item{Rs}{dB of stopband ripple.}

\item{w}{critical frequencies of the filter. \code{w} must be a scalar for
low-pass and high-pass filters, and \code{w} must be a two-element vector
\code{c(low, high)} specifying the lower and upper bands in radians/second.
For digital filters, w must be between 0 and 1 where 1 is the Nyquist
frequency.}

\item{type}{filter type, one of \code{"low"}, \code{"high"}, \code{"stop"},
or \code{"pass"}.}

\item{plane}{"z" for a digital filter or "s" for an analog filter.}
}
\value{
list of class \code{'\link{Arma}'} with list elements:
\describe{
  \item{b}{moving average (MA) polynomial coefficients}
  \item{a}{autoregressive (AR) polynomial coefficients}
}
}
\description{
Compute the transfer function coefficients of an elliptic filter.
}
\details{
An elliptic filter is a filter with equalized ripple (equiripple) behavior in
both the passband and the stopband. The amount of ripple in each band is
independently adjustable, and no other filter of equal order can have a
faster transition in gain between the passband and the stopband, for the
given values of ripple.

As the ripple in the stopband approaches zero, the filter becomes a type I
Chebyshev filter. As the ripple in the passband approaches zero, the filter
becomes a type II Chebyshev filter and finally, as both ripple values
approach zero, the filter becomes a Butterworth filter.

Because \code{ellip} is generic, it can be extended to accept other inputs, using
\code{ellipord} to generate filter criteria for example.
}
\examples{
# compare the frequency responses of 5th-order Butterworth and elliptic filters.
bf <- butter(5, 0.1)
ef <- ellip(5, 3, 40, 0.1)
bfr <- freqz(bf)
efr <- freqz(ef)
plot(bfr$w, 20 * log10(abs(bfr$h)), type = "l", ylim = c(-80, 0),
     xlab = "Frequency (Rad)", ylab = c("dB"))
lines(efr$w, 20 * log10(abs(efr$h)), col = "red")

}
\references{
\url{https://en.wikipedia.org/wiki/Elliptic_filter}
}
\seealso{
\code{\link{Arma}}, \code{\link{filter}}, \code{\link{butter}}, \code{\link{cheby1}}, \code{\link{ellipord}}
}
\author{
Paulo Neis, \email{p_neis@yahoo.com.br},\cr
  adapted by Doug Stewart, \email{dastew@sympatico.ca}.\cr
  Conversion to R Tom Short,\cr
  adapted by Geert van Boxtel, \email{G.J.M.vanBoxtel@gmail.com}.
}
