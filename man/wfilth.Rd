\name{wfilth}
\alias{wfilth}
\title{FFT of Wavelet Transform Function}
\description{
This function outputs the FFT of the wavelet. 
This code is translated from MATLAB Synchrosqueezing Toolbox, version 1.1 developed by Eugene Brevdo (http://www.math.princeton.edu/~ebrevdo/).
}
\usage{
wfilth(type, N, a, opt)
}
\arguments{
\item{type}{wavelet type. See \code{wfiltfn}. `gauss', `mhat', `cmhat', `morlet', `shannon', `hshannon', `hhat', `hhhat', `bump'}
\item{N}{number of samples to calculate}
\item{a}{wavelet scale parameter (default = 1)}
\item{opt}{list of options for wavelet type. See \code{wfiltfn}.
\code{opt.dt}: sampling period, default = 1.}
}
\details{
This function outputs the FFT of the wavelet of family 'type' with parameters
in 'opt', of length N at scale a: (psi(-t/a)).

Note that the output is made so that the inverse fft of the
result is zero-centered in time.  This is important for
convolving with the derivative(dpsih).  To get the correct
output, perform an \code{ifftshift}. That is,

\code{psi = ifftshift(fft(psih, inverse=TRUE) / length(psih))}, 

\code{xfilt = ifftshift(fft(fft(x) * psih, inverse=TRUE) / length(fft(x) * psih))}
}
\value{ 
\item{psih}{wavelet sampling in frequency domain (for use in \code{fft})}
\item{dpsih}{derivative of same wavelet, sampled in frequency domain (for \code{fft})}
\item{xi}{associated fourier domain frequencies of the samples.}
}
\seealso{
\code{\link{wfiltfn}}, \code{\link{cwt_fw}}, \code{\link{cwt_iw}}.
}
\examples{
tmp <- wfilth("morlet", 1024, 4)
plot(fftshift(tmp$xi/(2*pi)), fftshift(abs(tmp$psih)), type="l", col="blue", xlab="", ylab="")
}
\keyword{nonparametric}
