% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clustersegment.R
\name{clustersegment}
\alias{clustersegment}
\title{Cluster Segments}
\usage{
clustersegment(x)
}
\arguments{
\item{x}{input data, specified as a numeric vector or matrix, coerced to
contain only 0's and 1's, i.e., every nonzero element in \code{x} will
be replaced by 1.}
}
\value{
A list of size \code{nr}, where \code{nr} is the number
  of rows in \code{x}. Each element of the list contains a matrix with two
  rows. The first row is the initial index of a sequence of 1’s and the
  second row is the end index of that sequence.
}
\description{
Calculate boundary indexes of clusters of 1’s.
}
\details{
The function calculates the initial index and end index of sequences of 1s
rising and falling phases of the signal in \code{x}. The clusters are sought
in the rows of the array \code{x}. The function works by finding the indexes
of jumps between consecutive values in the rows of \code{x}.
}
\examples{
(x <- c(0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1))
(ranges <- clustersegment(x))
# The first sequence of 1's in x lies in the interval
(r <- ranges[1,1]:ranges[2,1])

x <- matrix(as.numeric(runif(30) > 0.4), 3, 10)
ranges <- clustersegment(x)

x <- c(0, 1.2, 3, -8, 0)
ranges <- clustersegment(x)

}
\author{
Juan Pablo Carbajal, \email{carbajal@ifi.uzh.ch}.\cr
 Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@gmail.com}.
}
