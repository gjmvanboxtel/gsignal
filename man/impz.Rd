% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/impz.R
\name{impz}
\alias{impz}
\alias{print.impz}
\alias{impz.Arma}
\alias{impz.Ma}
\alias{impz.Sos}
\alias{impz.Zpg}
\alias{impz.default}
\title{Impulse response of digital filter}
\usage{
impz(filt, ...)

\method{print}{impz}(x, ...)

\method{impz}{Arma}(filt, ...)

\method{impz}{Ma}(filt, ...)

\method{impz}{Sos}(filt, ...)

\method{impz}{Zpg}(filt, ...)

\method{impz}{default}(filt, a = 1, n = NULL, fs = 1, ...)
}
\arguments{
\item{filt}{for the default case, the moving-average coefficients of an ARMA
model or filter. Generically, \code{filt} specifies an arbitrary model or
filter operation.}

\item{...}{for methods of \code{freqz}, arguments are passed to the default
method. For \code{plot.impz}, additional arguments are passed through to
plot.}

\item{x}{object to be printed or plotted.}

\item{a}{the autoregressive (recursive) coefficients of an ARMA filter.}

\item{n}{number of points at which to evaluate the frequency response. If
\code{n} is a vector with a length greater than 1, then evaluate the
frequency response at these points. For fastest computation, \code{n}
should factor into a small number of small primes. Default: 512.}

\item{fs}{sampling frequency in Hz. If not specified (default = 2 * pi), the
frequencies are in radians.}
}
\value{
For \code{impz}, a list of class \code{"impz"} with items:
\describe{
  \item{x}{impulse response signal.}
  \item{t}{time.}
}
}
\description{
Compute the z-plane impulse response of an ARMA model or rational IIR
filter. A plot of the impulse and step responses is generated.
}
\note{
When results of \code{impz} are printed, \code{plot} will be called to
  display a plot of the impulse response against frequency. As with lattice
  plots, automatic printing does not work inside loops and function calls, so
  explicit calls to print or plot are needed there.
}
\examples{
## elliptic low-pass filter
elp <- ellip(4, 0.5, 20, 0.4)
impz(elp)

xt <- impz(elp)

}
\author{
Paul Kienzle, \email{pkienzle@users.sf.net}.\cr
Conversion to R by Tom Short;\cr
 adapted by Geert van Boxtel, \email{gjmvanboxtel@gmail.com}
}
