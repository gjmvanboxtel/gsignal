% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/buttord.R
\name{buttord}
\alias{buttord}
\title{Butterworth filter order and cutoff frequency}
\usage{
buttord(Wp, Ws, Rp, Rs, plane = c("z", "s"))
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
  \item{type}{filter type, normally one of "low", "high", "stop", or "pass".}
}
}
\description{
Compute the minimum filter order of a Butterworth filter with the desired
response characteristics.
}
\details{
Deriving the order and cutoff is based on:
\if{latex}{
  \deqn{|H(W)|^{2} = 1/[1 + (W / Wc)^{(2N)}] = 10^{(-R / 10)}}
}
\if{html}{\preformatted{
       2                  (2N)       (-R / 10)
 |H(W)|  = 1/[1 + (W / Wc)    ] = 10

}}

With some algebra, you can solve simultaneously for \code{Wc} and \code{N}
given \code{Ws}, \code{Rs} and Wp,Rp. Rounding N to the next greater integer,
one can recalculate the allowable range for \code{Wc} (filter characteristic
touching the pass band edge or the stop band edge).

For other types of filter, before making the above calculation, the
requirements must be transformed to lowpass requirements. After the
calculation, \code{Wc} must be transformed back to the original filter type.
}
\examples{
## low-pass 30 Hz filter
fs <- 128
butspec <- buttord(30/(fs/2), 40/(fs/2), 0.5, 40)
but <- butter(butspec)
freqz(but, fs = fs)

}
\seealso{
\code{\link{butter}}, \code{\link{FilterSpecs}}
}
\author{
Paul Kienzle,\cr
 adapted by Charles Praplan.\cr
 Conversion to R by Tom Short,\cr
 adapted by Geert van Boxtel,  \email{G.J.M.vanBoxtel@gmail.com}.
}
