\name{fftshift}
\alias{fftshift}
\title{FFT Shift}
\description{
This function exchanges the left halves of a vector with the right halves.
}
\usage{
fftshift(x)
}
\arguments{
\item{x}{a vector}
}
\details{
This function exchanges the left halves of a vector with the right halves. This function is adapted from Matlab.
}
\value{ 
shifted vector
}
\seealso{
\code{\link{ifftshift}}.
}
\examples{
x <- 1:4
fftshift(fftshift(x))
ifftshift(fftshift(x))

x <- 1:5
fftshift(fftshift(x))
ifftshift(fftshift(x))
}
\keyword{nonparametric}
