% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/butter.R
\name{butter}
\alias{butter}
\alias{butter.FilterSpecs}
\alias{butter.default}
\title{Butterworth filter design}
\usage{
butter(n, ...)

\method{butter}{FilterSpecs}(n, ...)

\method{butter}{default}(
  n,
  w,
  type = c("low", "high", "stop", "pass"),
  plane = c("z", "s"),
  output = c("Arma", "Zpg", "Sos"),
  ...
)
}
\arguments{
\item{n}{filter order.}

\item{...}{additional arguments passed to butter, overriding those given by
\code{n} of class \code{\link{FilterSpecs}}.}

\item{w}{critical frequencies of the filter. \code{w} must be a scalar for
low-pass and high-pass filters, and \code{w} must be a two-element vector
c(low, high) specifying the lower and upper bands in radians/second. For
digital filters, w must be between 0 and 1 where 1 is the Nyquist
frequency.}

\item{type}{filter type, one of \code{"low"}, (default) \code{"high"},
\code{"stop"}, or \code{"pass"}.}

\item{plane}{"z" for a digital filter or "s" for an analog filter.}

\item{output}{Type of output, one of:
\describe{
  \item{"Arma"}{Autoregressive-Moving average (aka numerator/denominator, aka
  b/a)}
  \item{"Zpg"}{Zero-pole-gain format}
  \item{"Sos"}{Second-order sections}
}
Default is \code{"Arma"} for compatibility with the 'signal' package and the
'Matlab' and 'Octave' equivalents, but \code{"Sos"} should be preferred for
general-purpose filtering because of numeric stability.}
}
\value{
Depending on the value of the \code{output} parameter, a list of
  class \code{\link{Arma}}, \code{\link{Zpg}}, or \code{\link{Sos}}
  containing the filter coefficients
}
\description{
Compute the transfer function coefficients of a Butterworth filter.
}
\details{
Butterworth filters have a magnitude response that is maximally flat in the
passband and monotonic overall. This smoothness comes at the price of
decreased rolloff steepness. Elliptic and Chebyshev filters generally provide
steeper rolloff for a given filter order.

Because butter is generic, it can be extended to accept other inputs, using
\code{buttord} to generate filter criteria for example.
}
\examples{
## 50 Hz notch filter
fs <- 256
bf <- butter(4, c(48, 52) / (fs / 2), "stop")
freqz(bf, fs = fs)

## EEG alpha rhythm (8 - 12 Hz) bandpass filter
fs <- 128
fpass <- c(8, 12)
wpass <- fpass / (fs / 2)
but <- butter(5, wpass, "pass")
freqz(but, fs = fs)

## filter to remove vocals from songs, 25 dB attenuation in stop band
## (not optimal with a Butterworth filter)
fs <- 44100
specs <- buttord(230/(fs/2), 450/(fs/2), 1, 25)
bf <- butter(specs)
freqz(bf, fs = fs)
zplane(bf)

}
\references{
\url{https://en.wikipedia.org/wiki/Butterworth_filter}
}
\seealso{
\code{\link{Arma}}, \code{\link{Zpg}}, \code{\link{Sos}},
  \code{\link{filter}}, \code{\link{cheby1}}, \code{\link{ellip}},
  \code{\link{buttord}}.
}
\author{
Paul Kienzle, \email{pkienzle@users.sf.net},\cr
  Doug Stewart, \email{dastew@sympatico.ca},\cr
  Alexander Klein, \email{alexander.klein@math.uni-giessen.de},\cr
  John W. Eaton.\cr
  Conversion to R by Tom Short,\cr
  adapted by Geert van Boxtel \email{G.J.M.vanBoxtel@gmail.com}.
}
