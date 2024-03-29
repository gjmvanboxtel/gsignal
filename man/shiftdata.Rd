% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shiftdata.R
\name{shiftdata}
\alias{shiftdata}
\title{Shift data to operate on specified dimension}
\usage{
shiftdata(x, dimx)
}
\arguments{
\item{x}{The data to be shifted. Can be of any type.}

\item{dimx}{Dimension of \code{x} to be shifted to the first column. Named
"dimx" instead of "dim" to avoid confusion with R's dim() function.
Default: NULL (shift the first nonsingleton dimension)}
}
\value{
A list containing 3 variables; \code{x}, the shifted data,
  \code{perm}, the permutation vector, and \code{nshifts}, the number of
  shifts
}
\description{
Shift data in to permute the dimension \code{dimx} to the first column.
}
\details{
\code{shiftdata(x, dimx)} shifts data \code{x} to permute dimension
\code{dimx} to the first column using the same permutation as the built-in
\code{filter} function. The vector \code{perm} in the output list returns the
permutation vector that is used.

If \code{dimx} is missing or empty, then the first nonsingleton dimension is
shifted to the first column, and the number of shifts is returned in
\code{nshifts}.

\code{shiftdata} is meant to be used in tandem with \code{unshiftdata}, which
shifts the data back to its original shape. These functions are useful for
creating functions that work along a certain dimension, like
\code{\link{filter}}, \code{\link{sgolayfilt}}, and \code{\link{sosfilt}}.
}
\examples{

## create a 3x3 magic square
x <- pracma::magic(3)
## Shift the matrix x to work along the second dimension.
## The permutation vector, perm, and the number of shifts, nshifts,
## are returned along with the shifted matrix.
sd <- shiftdata(x, 2)

## Shift the matrix back to its original shape.
y <- unshiftdata(sd)

## Rearrange Array to Operate on First nonsingleton Dimension
x <- 1:5
sd <- shiftdata(x)
y <- unshiftdata(sd)

}
\seealso{
\code{\link{unshiftdata}}
}
\author{
Georgios Ouzounis, \email{ouzounis_georgios@hotmail.com}.\cr
Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@gmail.com}.
}
