% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/uencode.R
\name{uencode}
\alias{uencode}
\title{Uniform encoder}
\usage{
uencode(u, n, v = 1, signed = FALSE)
}
\arguments{
\item{u}{Input, a multidimensional array of numbers, real or complex, single
or double precision.}

\item{n}{Number of levels used in \eqn{2^{n}}-level quantization. \code{n}
must be between 2 and 32}

\item{v}{Limit on the range of \code{u} to the range from \code{-v} to
\code{v} before saturating them. Default 1.}

\item{signed}{Logical indicating signed or unsigned output. See Details.
Default: FALSE.}
}
\value{
Multidimensional array of the same size as \code{u} containing signed
  or unsigned integers.
}
\description{
Quantize and encode floating-point inputs to integer outputs.
}
\details{
\code{y <- uencode(u, n)} quantizes the entries in a multidimensional array
of floating-point numbers \code{u} and encodes them as integers using
\eqn{2^{n}}-level quantization. \code{n} must be an integer between 2 and 32
(inclusive). Inputs can be real or complex, double- or single-precision. The
output \code{y} and the input \code{u} are arrays of the same size. The
elements of the output \code{y} are unsigned integers with magnitudes in the
range 0 to  \eqn{2^{n} - 1}. Elements of the input \code{u} outside of the
range -1 to 1 are treated as overflows and are saturated.
\itemize{
  \item For entries in the input u that are less than -1, the value of the
  output of uencode is 0.
  \item For entries in the input u that are greater than 1, the value of the
  output of uencode is \eqn{2^{n}-1}.
}

\code{y <- uencode(u, n, v)} allows the input \code{u} to have entries with
floating-point values in the range \code{-v} to \code{v} before saturating
them (the default value for \code{v} is 1). Elements of the input \code{u}
outside of the range \code{-v} to \code{v} are treated as overflows and are
saturated:
\itemize{
  \item For input entries less than \code{-v}, the value of the output of
  uencode is 0.
  \item For input entries greater than \code{v}, the value of the output of
  uencode is \eqn{2^{n} - 1}.
}

\code{y <- uencode(u, n, v, signed)} maps entries in a multidimensional array
of floating-point numbers \code{u} whose entries have values in the range
\code{-v} to \code{v} to an integer output \code{y}. Input entries outside
this range are saturated. The integer type of the output depends on the
number of quantization levels \eqn{2^{n}} and the value of \code{signed},
which can be one of the following:
\itemize{
  \item TRUE: Outputs are signed integers with magnitudes in the range
  \eqn{-2^{n} / 2} to \eqn{(2^{n} / 2) - 1}.
  \item FALSE (default): Outputs are unsigned integers with magnitudes in the
  range 0 to \eqn{2^{n} - 1}.
}
}
\examples{

u <- seq(-1, 1, 0.01)
y <- uencode(u, 3)
plot(u, y)

}
\author{
Georgios Ouzounis, \email{ouzounis_georgios@hotmail.com}.\cr
Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@gmail.com}.
}
