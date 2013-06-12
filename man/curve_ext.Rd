\name{curve_ext}
\alias{curve_ext}
\title{Extract a Maximum Energy / Minimum Curvature Curve}
\description{
This function extracts a maximum energy / minimum curvature curve from
Synchrosqueezed Representation.
This code is translated from MATLAB Synchrosqueezing Toolbox, version 1.1 developed by Eugene Brevdo (http://www.math.princeton.edu/~ebrevdo/).
}
\usage{
curve_ext(Tx, fs, lambda=0)
}
\arguments{
\item{Tx}{synchrosqueezed output of x (columns associated with time t)}
\item{fs}{frequencies associated with rows of \code{Tx}}
\item{lambda}{lambda should be greater than or equal to 0. Default: lambda=0}
}
\details{
This function extracts a maximum energy, minimum curvature, curve from
Synchrosqueezed Representation. Note, energy is given as:
\code{abs(Tx)^2}.

This implements the solution to Eq. (8) of [1].

Original Author: Jianfeng Lu
}
\value{  
\item{C}{the curve locations (indices)}
\item{E}{the (logarithmic) energy of the curve} 
}
\references{ 
[1] Thakur, G., Brevdo, E., Fuckar, N. S. and Wu, H-T. (2013) The Synchrosqueezing algorithm for time-varying spectral analysis: 
Robustness properties and new paleoclimate applications. 
\emph{Signal Processing}, \bold{93}, 1079--1094.
}
\seealso{
\code{\link{synsq_cwt_fw}}, \code{\link{curve_ext_multi}}, \code{\link{curve_ext_recon}}.
}
\keyword{nonparametric}
