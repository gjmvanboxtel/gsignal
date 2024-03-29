% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/xcorr2.R
\name{xcorr2}
\alias{xcorr2}
\title{2-D cross-correlation}
\usage{
xcorr2(a, b = a, scale = c("none", "biased", "unbiased", "coeff"))
}
\arguments{
\item{a}{Input matrix, coerced to numeric. Must not be missing.}

\item{b}{Input matrix, coerced to numeric. Default: \code{a}.}

\item{scale}{Character string. Specifies the type of scaling applied to the
correlation matrix. matched to one of:
\describe{
  \item{"none"}{no scaling}
  \item{"biased"}{Scales the raw cross-correlation by the maximum number of
  elements of \code{a} and \code{b} involved in the generation of any
  element of the output matrix.}
  \item{"unbiased"}{Scales the raw correlation by dividing each element in
  the cross-correlation matrix by the number of products \code{a} and
  \code{b} used to generate that element. }
  \item{"coeff"}{Scales the normalized cross-correlation on the range of [0
  1] so that a value of 1 corresponds to a correlation coefficient of 1. }
}}
}
\value{
2-D cross-correlation or autocorrelation matrix, returned as a matrix
}
\description{
Compute the 2D cross-correlation of matrices \code{a} and \code{b}.
}
\details{
If \code{b} is not specified, computes autocorrelation of \code{a},
i.e., same as \code{xcorr2 (a, a)}.
}
\examples{
m1 <- matrix(c(17, 24,  1,  8, 15,
               23,  5,  7, 14, 16,
                4,  6, 13, 20, 22,
               10, 12, 19, 21,  3,
               11, 18, 25,  2,  9), 5, 5, byrow = TRUE)
m2 <- matrix(c(8, 1, 6,
               3, 5, 7,
               4, 9, 2), 3, 3, byrow = TRUE)
R <- xcorr2(m1, m2)

}
\seealso{
\code{\link{conv2}}, \code{\link{xcorr}}.
}
\author{
Dave Cogdell, \email{cogdelld@asme.org},\cr
Paul Kienzle, \email{pkienzle@users.sf.net},\cr
Carne Draug, \email{carandraug+dev@gmail.com}.\cr
Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@gmail.com}.
}
