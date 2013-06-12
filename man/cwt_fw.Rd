\name{cwt_fw}
\alias{cwt_fw}
\title{Forward Continuous Wavelet Transform}
\description{
This function performs 
forward continuous wavelet transform, discretized, as described
in Sec. 4.3.3 of [1] and Sec. IIIA of [2]. 
This code is translated from MATLAB Synchrosqueezing Toolbox, version 1.1 developed by Eugene Brevdo (http://www.math.princeton.edu/~ebrevdo/).
}
\usage{
cwt_fw(x, type, nv, dt, opt)
}
\arguments{
\item{x}{input signal vector, length \code{n} (need not be dyadic length)}
\item{type}{wavelet type, string (see \code{wfiltfn})}
\item{nv}{number of voices (suggest \code{nv >= 32})}
\item{dt}{sampling period (default, \code{dt = 1})}
\item{opt}{list of options. \code{opt$padtype}: type of padding, `symmetric', `replicate', or `circular' (default = `symmetric').
\code{opt$rpadded}: return padded \code{Wx} and \code{dWx}?  (default = 0).
\code{opt$s}, \code{opt$mu}, etc: wavelet options (see \code{wfiltfn}).}
}
\details{
This function performs 
forward continuous wavelet transform, discretized, as described
in Sec. 4.3.3 of [1] and Sec. IIIA of [2].  This algorithm uses
the FFT and samples the wavelet atoms in the fourier domain.
Options such as padding of the original signal are allowed.
Returns the vector of scales and, if requested, the analytic
time-derivative of the wavelet transform (as described in
Sec. IIIB of [2].
}
\value{ 
\item{Wx}{\code{[na x n]} size matrix (rows = scales, cols = times) containing samples of the CWT of \code{x}.}
\item{asc}{\code{na} length vector containing the associated scales.}
\item{dWx}{\code{[na x n]} size matrix containing
samples of the time-derivatives of the CWT of \code{x}.} 
}
\references{
[1] Mallat, S (2009) \emph{A Wavelet Tour of Signal Processing}, New York: Academic Press.

[2] Thakur, G., Brevdo, E., Fuckar, N. S. and Wu, H-T. (2013) The Synchrosqueezing algorithm for time-varying spectral analysis: 
Robustness properties and new paleoclimate applications. 
\emph{Signal Processing}, \bold{93}, 1079--1094.
}
\seealso{
\code{\link{cwt_iw}}, \code{\link{wfiltfn}}, \code{\link{est_riskshrink_thresh}}.
}
\examples{
tt <- seq(0, 10, , 1024)
f0 <- (1+0.6*cos(2*tt))*cos(4*pi*tt+1.2*tt^2)
sigma <- 0.5
f <- f0 + sigma*rnorm(length(tt))

# Continuous wavelet transform
nv <- 32
opt <- list(type = "bump")
cwtfit <- cwt_fw(f, opt$type, nv, tt[2]-tt[1], opt)
}
\keyword{nonparametric}
