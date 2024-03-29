% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cceps.R
\name{cceps}
\alias{cceps}
\title{Complex cepstral analysis}
\usage{
cceps(x)
}
\arguments{
\item{x}{input data, specified as a real vector.}
}
\value{
Complex cepstrum, returned as a vector.
}
\description{
Return the complex cepstrum of the input vector.
}
\details{
Cepstral analysis is a nonlinear signal processing technique that is applied
most commonly in speech and image processing, or as a tool to investigate
periodic structures within frequency spectra, for instance resulting from
echos/reflections in the signal or to the occurrence of harmonic frequencies
(partials, overtones).

The cepstrum is used in many variants. Most important are the power cepstrum,
the complex cepstrum, and real cepstrum. The function \code{cceps} implements
the complex cepstrum by computing the inverse of the log-transformed FFT,
i.e.,

\deqn{cceps(x) <- ifft(log(fft(x)))}

However, because taking the logarithm of a complex number can lead to
unexpected results, the phase of \code{fft(x)} needs to be unwrapped before
taking the log.
}
\note{
This function returns slightly different results in comparison with the
  'Matlab' and 'Octave' equivalents. The 'Octave' version does not apply phase
  unwrapping, but has an optional correction procedure in case of zero phase
  at \eqn{\pi} radians. The present implementation does apply phase
  unwrapping so that the correction procedure is unnecessary. The 'Matlab'
  implementation also applies phase unwrapping, and a circular shift if
  necessary to avoid zero phase at \eqn{\pi} radians. The circular shift is
  not done here. In addition, the 'Octave' version shifts the zero frequency to
  the center of the series, which neither the 'Matlab' nor the present
  implementation do.
}
\examples{
## Generate a sine of frequency 45 Hz, sampled at 100 Hz.
fs <- 100
t <- seq(0, 1.27, 1/fs)
s1 <- sin(2 * pi * 45 * t)
## Add an echo with half the amplitude and 0.2 s later.
s2 <- s1 + 0.5 * c(rep(0L, 20), s1[1:108])
## Compute the complex cepstrum of the signal. Notice the echo at 0.2 s.
cep <- cceps(s2)
plot(t, cep, type="l")

}
\references{
\url{https://en.wikipedia.org/wiki/Cepstrum}
}
\seealso{
\code{\link{rceps}}
}
\author{
Geert van Boxtel, \email{G.J.M.vanBoxtel@gmail.com}.
}
