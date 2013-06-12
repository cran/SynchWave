\name{cwt_iw}
\alias{cwt_iw}
\title{Inverse Wavelet Transform}
\description{
This function performs 
the inverse wavelet transform of signal \code{Wx}.

This code is translated from MATLAB Synchrosqueezing Toolbox, version 1.1 developed by Eugene Brevdo (http://www.math.princeton.edu/~ebrevdo/).
}
\usage{
cwt_iw(Wx, type, opt=NULL)
}
\arguments{
\item{Wx}{wavelet transform of a signal, see \code{cwt_fw}.}
\item{type}{wavelet used to take the wavelet transform,
see \code{cwt_fw} and \code{wfiltfn}.}
\item{opt}{options list used for forward wavelet transform.}
}
\details{
This function performs 
the inverse wavelet transform of signal \code{Wx}, and 
implements Eq. (4.67) of [1].
}
\value{ 
\item{x}{the signal, as reconstructed from \code{Wx}.}
}
\references{
[1] Mallat, S (2009) \emph{A Wavelet Tour of Signal Processing}, New York: Academic Press.
}
\seealso{
\code{\link{cwt_fw}}, \code{\link{wfiltfn}}, \code{\link{est_riskshrink_thresh}}.
}
\examples{
tt <- seq(0, 10, , 1024)
nv <- 32
f0 <- (1+0.6*cos(2*tt))*cos(4*pi*tt+1.2*tt^2)
sigma <- 0.5
f <- f0 + sigma*rnorm(length(tt))

# Continuous wavelet transform
opt <- list(type = "bump")
cwtfit <- cwt_fw(f, opt$type, nv, tt[2]-tt[1], opt)

# Hard thresholing
thresh <- est_riskshrink_thresh(cwtfit$Wx, nv)
cwtfit$Wx[which(abs(cwtfit$Wx) < thresh)] <- 0.0

# Reconstruction 
opt$gamma <- thresh
cwtrec <- cwt_iw(cwtfit$Wx, opt$type, opt)

par(mfrow=c(1,1))
plot(tt, f, type="p", lty=2, xlab="time", ylab="f", col="red", cex=0.1)
lines(tt, f0, col="blue")
lines(tt, cwtrec)
}
\keyword{nonparametric}
