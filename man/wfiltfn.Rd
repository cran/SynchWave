\name{wfiltfn}
\alias{wfiltfn}
\title{Wavelet Transform Function of the Wavelet Filter}
\description{
This function provides 
wavelet transform function of the wavelet filter in question,
fourier domain. 
This code is translated from MATLAB Synchrosqueezing Toolbox, version 1.1 developed by Eugene Brevdo (http://www.math.princeton.edu/~ebrevdo/).
}
\usage{
wfiltfn(type, opt)
}
\arguments{
\item{type}{wavelet type. `gauss', `mhat', `cmhat', `morlet', `shannon', `hshannon', `hhat', `hhhat', `bump'}
\item{opt}{list of options for wavelet type. 
For `gauss', \code{s} and \code{mu} with default 1/6 and 2. For `mhat', \code{s} with default 1. For `cmhat', \code{s} and \code{mu} with default 1 and 1.
For `morlet', \code{mu} with default \code{2*pi}. For `hhhat', \code{mu} with default 5. For `bump', \code{s} and \code{mu} with default 1 and 5.
The wavelet types `mhat', `shannon' and `hhat' are not used for synchrosqueezed transform.
The wavelet types `shannon' and `hshannon' are NOT recommended for analysis.}
}
\details{
This function provides 
wavelet transform function of the wavelet filter in question,
fourier domain.
}
\value{ 
wavelet transform function
}
\seealso{
\code{\link{wfilth}}, \code{\link{cwt_fw}}, \code{\link{cwt_iw}}.
}
\examples{
psihfn <- wfiltfn("bump", list(mu=1, s=.5))
plot(psihfn(seq(-5, 5, by=0.01)), type="l", xlab="", ylab="")
}
\keyword{nonparametric}
