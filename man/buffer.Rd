% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/buffer.R
\name{buffer}
\alias{buffer}
\title{Buffer signal vector into matrix of data segments}
\usage{
buffer(x, n, p = 0, opt, zopt = FALSE)
}
\arguments{
\item{x}{The data to be buffered.}

\item{n}{The number of rows in the produced data buffer. This is an positive
integer value and must be supplied.}

\item{p}{An integer less than \code{n} that specifies the under- or overlap
between column in the data frame. Default 0.}

\item{opt}{In the case of an overlap, \code{opt} can be either a vector of
length \code{p} or the string \code{'nodelay'}. If \code{opt} is a vector,
then the first \code{p} entries in \code{y} will be filled with these
values. If \code{opt} is the string \code{'nodelay'}, then the first value
of \code{y} corresponds to the first value of \code{x}. In the case of an
underlap, \code{opt} must be an integer between 0 and \code{-p}. The
represents the initial underlap of the first \code{y}. The default value
for \code{opt} the vector \code{matrix (0L, 1, p)} in the case of an
overlap, or 0 otherwise.}

\item{zopt}{Logical. If TRUE, return values for \code{z} and \code{opt} in
addition to \code{y}. Default is FALSE (return only \code{y}).}
}
\value{
If \code{zopt} equals FALSE (the default), this function returns a
  single numerical array containing the buffered data (\code{y}). If
  \code{zopt} equals TRUE, then a \code{list} containing 3 variables is
  returned: \code{y}: the buffered data, \code{z}: the over or underlap (if
  any), \code{opt}: the over- or underlap that might be used for a future
  call to \code{buffer} to allow continuous buffering.
}
\description{
Partition a signal vector into nonoverlapping, overlapping, or underlapping
data segments.
}
\details{
\code{y <- buffer(x, n)} partitions a signal vector \code{x} of length
\code{L} into nonoverlapping data segments of length \code{n}. Each data
segment occupies one column of matrix output \code{y}, which has \code{n}
rows and \code{ceil(L / n)} columns. If \code{L} is not evenly divisible by
\code{n}, the last column is zero-padded to length \code{n}.

\code{y <- buffer(x, n, p)} overlaps or underlaps successive frames in the
output matrix by \code{p} samples.
\itemize{
\item {For \code{0 < p < n} (overlap), buffer repeats the final \code{p}
samples of each segment at the beginning of the following segment. See the
example where \code{x = 1:30}, \code{n = 7}, and an overlap of \code{p = 3}.
In this case, the first segment starts with \code{p} zeros (the default
initial condition), and the number of columns in \code{y} is \code{ceil(L /
(n - p))}.}
\item  {For \code{p < 0} (underlap), buffer skips \code{p} samples between
consecutive segments. See the example where \code{x = 1:30}, \code{n = 7},
and \code{p = -3}. The number of columns in \code{y} is \code{ceil(L / (n -
p))}.}
}

In \code{y <- buffer(x, n, p, opt)}, \code{opt} specifies a vector of samples
to precede \code{x[1]} in an overlapping buffer, or the number of initial
samples to skip in an underlapping buffer.
\itemize{
  \item {For \code{0 < p < n} (overlap), \code{opt} specifies a vector of
  length \code{p} to insert before \code{x[1]} in the buffer. This vector can
  be considered an initial condition, which is needed when the current
  buffering operation is one in a sequence of consecutive buffering
  operations. To maintain the desired segment overlap from one buffer to the
  next, \code{opt} should contain the final \code{p} samples of the previous
  buffer in the sequence. Set \code{opt} to \code{"nodelay"} to skip the
  initial condition and begin filling the buffer immediately with
  \code{x[1]}. In this case, \code{L} must be \code{length(p)} or longer. See
  the example where \code{x = 1:30}, \code{n = 7}, \code{p = 3}, and
  \code{opt = "nodelay"}.}
  \item {For \code{p < 0} (underlap), \code{opt} is an integer value in the
    range \code{0 : -p} specifying the number of initial input samples,
    \code{x[1:opt]}, to skip before adding samples to the buffer. The first
    value in the buffer is therefore \code{x[opt + 1]}.}
}
The \code{opt} option is especially useful when the current buffering
operation is one in a sequence of consecutive buffering operations. To
maintain the desired frame underlap from one buffer to the next, \code{opt}
should equal the difference between the total number of points to skip
between frames (\code{p}) and the number of points that were available to be
skipped in the previous input to buffer. If the previous input had fewer than
p points that could be skipped after filling the final frame of that buffer,
the remaining opt points need to be removed from the first frame of the
current buffer. See Continuous Buffering for an example of how this works in
practice.

\code{buf <- buffer(..., zopt = TRUE)} returns the last \code{p} samples of a
overlapping buffer in output \code{buf$opt}. In an underlapping buffer,
\code{buf$opt} is the difference between the total number of points to skip
between frames (\code{-p}) and the number of points in \code{x} that were
available to be skipped after filling the last frame:
\itemize{
  \item {For \code{0 < p < n} (overlap), \code{buf$opt} contains the final
  \code{p} samples in the last frame of the buffer. This vector can be used
  as the initial condition for a subsequent buffering operation in a sequence
  of consecutive buffering operations. This allows the desired frame overlap
  to be maintained from one buffer to the next. See Continuous Buffering
  below.}
  \item {For \code{p < 0} (underlap), \code{buf$opt} is the difference
  between the total number of points to skip between frames \code{(-p)} and
  the number of points in \code{x} that were available to be skipped after
  filling the last frame: \code{buf$opt = m*(n-p) + opt - L} where \code{opt}
  on the right is the input argument to buffer, and \code{buf$opt} on the
  left is the output argument. Note that for an underlapping buffer output
  \code{buf$opt} is always zero when output \code{buf$z} contains data.\cr
  The opt output for an underlapping buffer is especially useful when the
  current buffering operation is one in a sequence of consecutive buffering
  operations. The \code{buf$opt} output from each buffering operation
  specifies the number of samples that need to be skipped at the start of the
  next buffering operation to maintain the desired frame underlap from one
  buffer to the next. If fewer than \code{p} points were available to be
  skipped after filling the final frame of the current buffer, the remaining
  opt points need to be removed from the first frame of the next buffer.}
}
In a sequence of buffering operations, the \code{buf$opt} output from each
operation should be used as the \code{opt} input to the subsequent buffering
operation. This ensures that the desired frame overlap or underlap is
maintained from buffer to buffer, as well as from frame to frame within the
same buffer. See Continuous Buffering below for an example of how this works
in practice.
\cr

\strong{Continuous Buffering}\cr\cr
In a continuous buffering operation, the vector input to the buffer function
represents one frame in a sequence of frames that make up a discrete signal.
These signal frames can originate in a frame-based data acquisition process,
or within a frame-based algorithm like the FFT.\cr
As an example, you might acquire data from an A/D card in frames of 64
samples. In the simplest case, you could rebuffer the data into frames of 16
samples; \code{buffer} with \code{n = 16} creates a buffer of four frames
from each 64-element input frame. The result is that the signal of frame size
64 has been converted to a signal of frame size 16; no samples were added or
removed.\cr
In the general case where the original signal frame size, \code{L}, is not
equally divisible by the new frame size, \code{n}, the overflow from the last
frame needs to be captured and recycled into the following buffer. You can do
this by iteratively calling buffer on input x with the \code{zopt} parameter
set to \code{TRUE}. This simply captures any buffer overflow in \code{buf$z},
and prepends the data to the subsequent input in the next call to buffer.\cr
Note that continuous buffering cannot be done without the \code{zopt}
parameter being set to \code{TRUE}, because the last frame of y (\code{buf$y}
in this case) is zero padded, which adds new samples to the signal.\cr
Continuous buffering in the presence of overlap and underlap is handled with
the \code{opt} parameter, which is used as both an input (\code{opt} and
output (\code{buf$opt}) to buffer. The two examples on this page demonstrate
how the \code{opt} parameter should be used.
}
\examples{
## Examples without continuous buffering
y <- buffer(1:10, 5)
y <- buffer(1:10, 4)
y <- buffer(1:30, 7, 3)
y <- buffer(1:30, 7, -3)
y <- buffer(1:30, 7, 3, 'nodelay')

## Continuous buffering examples
# with overlap:
data <- buffer(1:1100, 11)
n <- 4
p <- 1
buf <- list(y = NULL, z = NULL, opt = -5)
for (i in 1:ncol(data)) {
  x <- data[,i]
  buf <- buffer(x = c(buf$z,x), n, p, opt=buf$opt, zopt = TRUE)
}
# with underlap:
data <- buffer(1:1100, 11)
n <- 4
p <- -2
buf <- list(y = NULL, z = NULL, opt = 1)
for (i in 1:ncol(data)) {
  x <- data[,i]
  buf <- buffer(x = c(buf$z,x), n, p, opt=buf$opt, zopt = TRUE)
}

}
\author{
David Bateman, \email{adb014@gmail.com}.\cr
 Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@gmail.com}
}
