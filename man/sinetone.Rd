% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sinetone.R
\name{sinetone}
\alias{sinetone}
\title{Sine tone}
\usage{
sinetone(freq, rate = 8000, sec = 1, ampl = 64)
}
\arguments{
\item{freq}{frequency of the tone, specified as a vector of positive numeric
values. The length of \code{freq} should equal the length of the
\code{ampl} vector; the shorter of the two is recycled to the longer
vector.}

\item{rate}{sampling frequency, specified as a positive scalar. Default:
8000.}

\item{sec}{length of the generated tone in seconds. Default: 1}

\item{ampl}{amplitude of the tone, specified as a vector of positive numeric
values. The length of \code{ampl} should equal the length of the
\code{freq} vector; the shorter of the two is recycled to the longer
vector. Default: 64.}
}
\value{
Sine tone, returned as a vector of length \code{rate * sec}, or as a
  matrix with \code{rate * sec} columns and  \code{max(length(freq),
  length(ampl))} columns.
}
\description{
Generate discrete sine tone.
}
\examples{
fs <- 1000
sec <- 2
y <- sinetone(10, fs, sec, 1)
plot(seq(0, sec, length.out = sec * fs), y, type = "l", xlab = "", ylab = "")

y <- sinetone(c(10, 15), fs, sec, c(1, 2))
matplot(seq(0, sec, length.out = sec * fs), y, type = "l",
        xlab = "", ylab = "")

}
\author{
Friedrich Leisch.\cr
Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@gmail.com}.
}
