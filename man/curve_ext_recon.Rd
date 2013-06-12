\name{curve_ext_recon}
\alias{curve_ext_recon}
\title{Reconstruct Curves}
\description{
This function reconstructs the curves given by \code{curve_ext} or \code{curve_ext_multi}.

This code is translated from MATLAB Synchrosqueezing Toolbox, version 1.1 developed by Eugene Brevdo (http://www.math.princeton.edu/~ebrevdo/).
}
\usage{
curve_ext_recon(Tx, fs, Cs, opt=NULL, clwin=4)
}
\arguments{
\item{Tx}{same as input to \code{curve_ext_multi}. synchrosqueezed output of \code{x} (columns associated with time \code{t})}
\item{fs}{same as input to \code{curve_ext_multi}. frequencies associated with rows of \code{Tx}}
\item{opt}{same as input to \code{synsq_cwt_fw} and \code{synsq_cwt_iw}}
\item{Cs}{same as output to \code{curve_ext_multi}. \code{[N x nc]} matrix of curve indices (\code{N=ncol(Tx)})}
\item{clwin}{same as input to \code{curve_ext_multi}. frequency clearing window; after curve extraction, a window
    of frequencies \code{(Cs[,i]-clwin) : (Cs[,i]+clwin)} is removed from the
    original representation. (default = 2)}
}
\details{
This function reconstructs the curves given by \code{curve_ext} or \code{curve_ext_multi}.
}
\value{ 
\item{xrs}{\code{[n x Nc]} matrix - the \code{Nc} reconstructed signals (each length \code{n})}
}
\references{
[1] Thakur, G., Brevdo, E., Fuckar, N. S. and Wu, H-T. (2013) The Synchrosqueezing algorithm for time-varying spectral analysis: 
Robustness properties and new paleoclimate applications. 
\emph{Signal Processing}, \bold{93}, 1079--1094.
}
\seealso{
\code{\link{synsq_cwt_fw}}, \code{\link{curve_ext}}, \code{\link{curve_ext_multi}}, \code{\link{curve_ext_recon}}.
}
\examples{
set.seed(7)
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

# Denoised signal
opt$gamma <- thresh
fr <- cwt_iw(cwtfit$Wx, opt$type, opt)

# Synchrosqueezed wavelet transform using denoised signal
sstfit2 <- synsq_cwt_fw(tt, fr, nv, opt)

# Ridge extraction
lambda <- 1e+04
nw <- 16
imtfit <- curve_ext_multi(sstfit2$Tx, log2(sstfit2$fs), 1, lambda, nw)

# Reconstruction
curvefit <- curve_ext_recon(sstfit2$Tx, sstfit2$fs, imtfit$Cs, opt, nw)

par(mfrow=c(2,1))
image.plot(list(x=tt, y=sstfit2$fs, z=t(abs(sstfit2$Tx))), log="y", 
    xlab="Time", ylab="Frequency", main="Time-Frequency Representation by SST", 
    col=designer.colors(64, c("azure", "cyan", "blue", "darkblue")), ylim=c(0.5, 25))
lines(tt, sstfit2$fs[imtfit$Cs[,1]], col="red", lty=3, lwd=2)

plot(tt, f0, type="l")
lines(tt, curvefit, lty=2)
}
\keyword{nonparametric}
