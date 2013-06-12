\name{synsq_cwt_fw}
\alias{synsq_cwt_fw}
\title{Synchrosqueezing Transform}
\description{
This function calculates the synchrosqueezing transform. 

This code is translated from MATLAB Synchrosqueezing Toolbox, version 1.1 developed by Eugene Brevdo (http://www.math.princeton.edu/~ebrevdo/).
}
\usage{
synsq_cwt_fw(tt, x, nv=16, opt=NULL) 
}
\arguments{
\item{tt}{vector of times samples are taken (e.g. \code{seq(0, 1,length=2)})}
\item{x}{vector of signal samples (e.g. \code{x = cos(20*pi*tt)})}
\item{nv}{number of voices (default: 16, recommended: 32 or 64)}
\item{opt}{list of options. \code{opt.type}: type of wavelet (see \code{wfiltfn}), \code{opt$s}, \code{opt$mu}, etc (wavelet parameters: see \code{opt} from \code{wfiltfn}),
\code{opt.disp}: display debug information?,
\code{opt.gamma}: wavelet hard thresholding value (see \code{cwt_freq_direct}, default: sqrt(machine epsilon))
\code{opt.dtype}: direct or numerical differentiation (default: 1, uses direct).  if \code{dtype=0}, uses differentiation algorithm 
instead (see \code{cwt_freq}), which is faster and
uses less memory, but may be less accurate.}
}
\details{
This function
calculates the synchrosqueezing transform of vector \code{x}, with
samples taken at times given in vector \code{tt}.  Uses \code{nv} voices.  \code{opt}
is a struct of options for the way synchrosqueezing is done, the
type of wavelet used, and the wavelet parameters.  This
implements the algorithm described in Sec. III of [1].

Note that \code{Wx} itself is not thresholded by \code{opt$gamma}.
The instantaneoues frequency \code{w} is calculated using \code{Wx} by \code{cwt_freq_direct} and
\code{cwt_freq}. After the calculation, instantaneoues frequency \code{w} is treated as \code{NA} where \code{abs(Wx) < opt$gamma}.
}
\value{ 
\item{Tx}{synchrosqueezed output of \code{x} (columns associated with time \code{tt})}
\item{fs}{frequencies associated with rows of \code{Tx}}
\item{Wx}{wavelet transform of \code{x} (see \code{cwt_fw})}
\item{asc}{scales associated with rows of \code{Wx}}
\item{w}{instantaneous frequency estimates for each element of \code{Wx}
(see \code{cwt_freq_direct})}
}
\references{
[1] Thakur, G., Brevdo, E., Fuckar, N. S. and Wu, H-T. (2013) The Synchrosqueezing algorithm for time-varying spectral analysis: 
Robustness properties and new paleoclimate applications. 
\emph{Signal Processing}, \bold{93}, 1079--1094.
}
\seealso{
\code{\link{cwt_fw}}, \code{\link{est_riskshrink_thresh}}, \code{\link{wfiltfn}}.
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
#[1] 0.0593984
#opt$gamma <- 10^-5
cwtrec <- cwt_iw(cwtfit$Wx, opt$type, opt)

par(mfrow=c(1,1))
plot(tt, f, type="p", lty=2, xlab="time", ylab="f", col="red", cex=0.1)
lines(tt, f0, col="blue")
lines(tt, cwtrec)

# Synchrosqueezed wavelet transform
sstfit <- synsq_cwt_fw(tt, f, nv, opt)

#par(mfrow=c(2,2))
#plot(tt, f, type="p", lty=2, xlab="time", ylab="f", col="red", cex=0.1)
#lines(tt, f0, col="blue")

#image.plot(list(x=tt, y=sstfit$asc, z=t(abs(sstfit$Wx))), log="y", 
#    xlab="Time", ylab="Scale", main="Time-Scale Representation by CWT",  
#    col=designer.colors(64, c("azure", "cyan", "blue", "darkblue")), ylim=rev(range(sstfit$asc)))
#image.plot(list(x=tt, y=sstfit$fs, z=t(abs(sstfit$Tx))), log="y", 
#    xlab="Time", ylab="Frequency", main="Time-Frequency Representation by SST", 
#    col=designer.colors(64, c("azure", "cyan", "blue", "darkblue")), ylim=c(0.5, 25))
#image.plot(list(x=tt, y=sstfit$asc, z=t(sstfit$w)), log="y", 
#    xlab="Time", ylab="Scale", main="Instantaneous Frequency",  
#    col=designer.colors(64, c("azure", "cyan", "blue", "darkblue")), ylim=rev(range(sstfit$asc)))
}
\keyword{nonparametric}
