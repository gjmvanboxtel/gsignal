% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tukeywin.R
\name{tukeywin}
\alias{tukeywin}
\title{Tukey (tapered cosine) window}
\usage{
tukeywin(n, r = 1/2)
}
\arguments{
\item{n}{Window length, specified as a positive integer.}

\item{r}{Cosine fraction, specified as a real scalar. The Tukey window is a
rectangular window with the first and last \code{r / 2} percent of the
samples equal to parts of a cosine. For example, setting \code{r = 0.5}
(default) produces a Tukey window where 1/2 of the entire window length
consists of segments of a phase-shifted cosine with period 2r = 1. If you
specify r <= 0, an n-point rectangular window is returned. If you specify r
>= 1, an n-point von Hann window is returned.}
}
\value{
Tukey window, returned as a vector.
}
\description{
Return the filter coefficients of a Tukey window (also known as the
cosine-tapered window) of length \code{n}.
}
\details{
The Tukey window, also known as the tapered cosine window, can be regarded as
a cosine lobe that is convolved with a rectangular window.  \code{r} defines
the ratio between the constant section and and the cosine section. It has to
be between 0 and 1. The function returns a Hann window for \code{r} equal to
1 and a rectangular window for \code{r} equal to 0.
}
\examples{

n <- 128
t0 <- tukeywin(n, 0)        # Equivalent to a rectangular window
t25 <- tukeywin(n, 0.25)
t5 <- tukeywin(n)           # default r = 0.5
t75 <- tukeywin(n, 0.75)
t1 <- tukeywin(n, 1)         # Equivalent to a Hann window
plot(t0, type = "l", xlab = "Samples", ylab =" Amplitude", ylim=c(0,1.2))
lines(t25, col = 2)
lines(t5, col = 3)
lines(t75, col = 4)
lines(t1, col = 5)

}
\author{
Laurent Mazet, \email{mazet@crm.mot.com}.\cr
Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@gmail.com}.
}
