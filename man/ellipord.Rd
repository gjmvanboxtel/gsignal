% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ellipord.R
\name{ellipord}
\alias{ellipord}
\title{Elliptic Filter Order}
\usage{
ellipord(Wp, Ws, Rp, Rs, plane = c("z", "s"))
}
\arguments{
\item{Wp, Ws}{pass-band and stop-band edges. For a low-pass or high-pass
filter, \code{Wp} and \code{Ws} are scalars. For a band-pass or
band-rejection filter, both are vectors of length 2. For a low-pass filter,
\code{Wp < Ws}. For a high-pass filter, \code{Ws > Wp}. For a band-pass
\code{(Ws[1] < Wp[1] < Wp[2] < Ws[2])} or band-reject \code{(Wp[1] < Ws[1]
< Ws[2] < Wp[2])} filter design, \code{Wp} gives the edges of the pass
band, and \code{Ws} gives the edges of the stop band. For digital filters,
frequencies are normalized to [0, 1], corresponding to the range [0, fs/2].
In case of an analog filter, all frequencies are specified in radians per
second.}

\item{Rp}{allowable decibels of ripple in the pass band.}

\item{Rs}{minimum attenuation in the stop band in dB.}

\item{plane}{"z" for a digital filter or "s" for an analog filter.}
}
\value{
A list of class \code{\link{FilterSpecs}} with the following list
  elements:
\describe{
  \item{n}{filter order}
  \item{Wc}{cutoff frequency}
  \item{type}{filter type, one of \code{"low"}, \code{"high"}, \code{"stop"},
  or \code{"pass"}.}
  \item{Rp}{dB of passband ripple.}
  \item{Rs}{dB of stopband ripple.}
}
}
\description{
Compute elliptic filter order and cutoff for the desired
response characteristics.
}
\examples{
fs <- 10000
spec <- ellipord(1000/(fs/2), 1200/(fs/2), 0.5, 29)
ef <- ellip(spec)
hf <- freqz(ef, fs = fs)
plot(c(0, 1000, 1000, 0, 0), c(0, 0, -0.5, -0.5, 0),
     type = "l", xlab = "Frequency (Hz)", ylab = "Attenuation (dB)",
     col = "red", ylim = c(-35,0), xlim = c(0,2000))
lines(c(5000, 1200, 1200, 5000, 5000), c(-1000, -1000, -29, -29, -1000),
      col = "red")
lines(hf$w, 20*log10(abs(hf$h)))

}
\seealso{
\code{\link{buttord}}, \code{\link{cheb1ord}},
  \code{\link{cheb2ord}}, \code{\link{ellip}}
}
\author{
Paulo Neis, \email{p_neis@yahoo.com.br},\cr
  adapted by Charles Praplan.\cr
  Conversion to R by Tom Short,\cr
  adapted by Geert van Boxtel, \email{G.J.M.vanBoxtel@gmail.com}.
}
