% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fir1.R
\name{fir1}
\alias{fir1}
\title{Window-based FIR filter design}
\usage{
fir1(
  n,
  w,
  type = c("low", "high", "stop", "pass", "DC-0", "DC-1"),
  window = hamming(n + 1),
  scale = TRUE
)
}
\arguments{
\item{n}{filter order (1 less than the length of the filter).}

\item{w}{band edges, strictly increasing vector in the range c(0, 1), where 1
is the Nyquist frequency. A scalar for highpass or lowpass filters, a
vector pair for bandpass or bandstop, or a vector for an alternating
pass/stop filter.}

\item{type}{character specifying filter type, one of \code{"low"} for a
low-pass filter, \code{"high"} for a high-pass filter, \code{"stop"} for a
stop-band (band-reject) filter, \code{"pass"} for a pass-band filter,
\code{"DC-0"} for a bandpass as the first band of a multiband filter, or
\code{"DC-1"} for a bandstop as the first band of a multiband filter.
Default: \code{"low"}.}

\item{window}{smoothing window. The returned filter is the same shape as the
smoothing window. Default: \code{hamming(n + 1)}.}

\item{scale}{whether to normalize or not. Use \code{TRUE} (default) or
\code{"scale"} to set the magnitude of the center of the first passband to
1, and \code{FALSE} or \code{"noscale"} to not normalize.}
}
\value{
The FIR filter coefficients, a vector of length \code{n + 1}, of
  class \code{Ma}.
}
\description{
FIR filter coefficients for a filter with the given order and frequency
cutoff.
}
\examples{
freqz(fir1(40, 0.3))
freqz(fir1(10, c(0.3, 0.5), "stop"))
freqz(fir1(10, c(0.3, 0.5), "pass"))

}
\references{
\url{https://en.wikipedia.org/wiki/Fir_filter}
}
\seealso{
\code{\link{Ma}}, \code{\link{filter}}, \code{\link{fftfilt}},
  \code{\link{fir2}}
}
\author{
Paul Kienzle, \email{pkienzle@users.sf.net},
 Conversion to R Tom Short,\cr
 adapted by Geert van Boxtel, \email{G.J.M.vanBoxtel@gmail.com}.
}
