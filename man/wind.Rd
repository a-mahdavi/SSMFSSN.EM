\name{wind}
\alias{wind}
\title{Wind speed data}
\usage{
data(wind)
}
\description{
Wind speed dataset
}
\examples{
 data(wind)
 y <- wind
 # EM output for MFSSN and MFSSTT distributions : 
 EM.SSMFSSN(y,xi=c(22,15,14), S=cov(y), la1=c(-1,1.2,1.3),la2=c(-.8,-.4,-.1))
 EM.SSMFSSN(y,xi=c(19,15,13), S=cov(y), la1=c(.3,1.2,1.3) , la2=c(3,-.4,-.1),
 nu1=5, nu2=2, family="MFSSTT" ,iter.max=500,tol=10^-9)
}
