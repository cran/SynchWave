\name{est_riskshrink_thresh}
\alias{est_riskshrink_thresh}
\title{Estimate the RiskShrink Hard Thresholding Level}
\description{
This function estimates the RiskShrink hard thresholding level.

This code is translated from MATLAB Synchrosqueezing Toolbox, version 1.1 developed by Eugene Brevdo (http://www.math.princeton.edu/~ebrevdo/).
}
\usage{
est_riskshrink_thresh(Wx, nv)
}
\arguments{
\item{Wx}{wavelet transform of a signal, see \code{cwt_fw}}
\item{nv}{number of voices}
}
\details{
This function implements Defn. 1 of Sec. 2.4 in [1], using the suggested
noise estimator from the discussion "Estimating the noise level"
in that same section.
}
\value{ 
the RiskShrink hard threshold estimate
}
\references{
[1] Donoho, D. L. and Johnstone, I. M. (1994) Ideal spatial adaptation by 
wavelet shrinkage. \emph{Biometrika}, \bold{81}, 425--455.
}
\seealso{
\code{\link{cwt_fw}}, \code{\link{cwt_iw}}.
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
}
\keyword{nonparametric}
