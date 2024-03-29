% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cl2bp.R
\name{cl2bp}
\alias{cl2bp}
\title{Constrained L2 bandpass FIR filter design}
\usage{
cl2bp(m = 30, w1, w2, up, lo, L = 2048)
}
\arguments{
\item{m}{degree of cosine polynomial, resulting in a filter of length
\code{2 * m + 1}. Must be an even number. Default: 30.}

\item{w1, w2}{bandpass filter cutoffs in the range \code{0 <= w1 < w2 <= pi},
where pi is the Nyquist frequency.}

\item{up}{vector of 3 upper bounds for c(stopband1, passband, stopband2).}

\item{lo}{vector of 3 lower bounds for c(stopband1, passband, stopband2).}

\item{L}{search grid size; larger values may improve accuracy, but greatly
increase calculation time. Default: 2048, maximum: 1e6.}
}
\value{
The FIR filter coefficients, a vector of length \code{2 * m + 1}, of
  class \code{Ma}.
}
\description{
Constrained least square band-pass FIR filter design without specified
transition bands.
}
\details{
This is a fast implementation of the algorithm cited below. Compared to
\code{remez}, it offers implicit specification of transition bands, a higher
likelihood of convergence, and an error criterion combining features of both
L2 and Chebyshev approaches
}
\examples{
w1 <- 0.3 * pi
w2 <- 0.6 * pi
up <- c(0.02, 1.02, 0.02)
lo <- c(-0.02, 0.98, -0.02)
h  <- cl2bp(30, w1, w2, up, lo, 2^11)
freqz(h)

}
\references{
Selesnick, I.W., Lang, M., and Burrus, C.S. (1998) A modified
  algorithm for constrained least square design of multiband FIR filters
  without specified transition bands. IEEE Trans. on Signal Processing,
  46(2), 497-501. \cr
  \url{https://www.ece.rice.edu/dsp/software/cl2.shtml}
}
\seealso{
\code{\link{Ma}}, \code{\link{filter}}, \code{\link{remez}}
}
\author{
Ivan Selesnick, Rice University, 1995,
  downloaded from \url{https://www.ece.rice.edu/dsp/software/cl2.shtml}.\cr
  Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@gmail.com}.
}
