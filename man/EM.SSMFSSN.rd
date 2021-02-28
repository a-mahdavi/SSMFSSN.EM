\name{EM.SSMFSSN}
\alias{EM.SSMFSSN}
\title{EM.SSMFSSN}
\usage{
EM.SSMFSSN(y, xi, S, la1, la2, nu1 , nu2, family="MFSSN", get.init = FALSE, iter.max=100,
 tol=10^-6, equal=FALSE)
}
\description{
 Fit the multivariate SSMFSSN distributions using EM-algorithm
 xi: the vector of location parameter.
 S: the cov-variance matrix.
 la1 and la2: the vector of shape parameters.
 nu1 and nu2: the flatness parameters.
 family: distribution family to be used in fitting ("MFSSN", "MFSSTN","MFSSLSN", "MFSSCN",
 "MFSST", "MFSSTT").
 get.init: if TRUE, the initial values are generated.
 iter.max: the maximum number of iterations of the EM algorithm. Default= 100.
 tol: the covergence maximum error. Default= 10^-6.
 equal: if TRUE the nu1 and nu2 assumed equal.
}
\examples{

 #  Example 1:
 # Simulating samples from MFSSTN distribution:
 y <- r.SSMFSSN(n=100, xi=c(0,5), S=matrix(c(1,.4,.4,4),2,2), la1=c(-2,3) , la2=c(.5,-.5),
 nu1=5, family="MFSSTN")
 # n: the number of random samples
 # EM output with specific initial values: 
 EM.SSMFSSN(y, xi=c(0,5), S=matrix(c(1,0.4,0.4,4),2,2), la1=c(-2,3) , la2=c(0.5,-0.5),
 family="MFSSTN", get.init=FALSE,  iter.max=100, tol=10^-6)
 # EM output without specific initial values: 
 EM.SSMFSSN(y, family="MFSSTN", get.init=TRUE)

 # Example 2:
 # Simulating samples from MFSSTT distribution:
 y <- r.SSMFSSN(n=100, xi=c(0,5), S=matrix(c(1,0.4,0.4,4),2,2), la1=c(-2,3) , la2=c(0.5,-0.5),
 nu1=5, nu2=10, family="MFSSTT")
 # EM output with specific initial values: 
 EM.SSMFSSN(y, xi=c(0,5), S=matrix(c(1,0.4,0.4,4),2,2), la1=c(-2,3) , la2=c(0.5,-0.5), nu1=5, nu2=10,
 family="MFSSTT", get.init=FALSE, equal=F)
 
 # Example 3:
 # Simulating samples from MFSSCNe distribution:
 y <- r.SSMFSSN(n=100, xi=c(0,5), S=matrix(c(1,0.4,0.4,4),2,2), la1=c(-2,3) , la2=c(0.5,-0.5),
 nu1=0.3, nu2=0.3, family="MFSSCN")
 # EM output assuming the equality for the flatness parameters : 
 EM.SSMFSSN(y, family="MFSSCN", get.init=TRUE, equal=T)

 # Example 4: wind speed data
 data(wind)
 y <- wind
 # EM output for MFSSN and MFSSTT distributions : 
 EM.SSMFSSN(y,xi=c(22,15,14), S=cov(y), la1=c(-1,1.2,1.3),la2=c(-.8,-.4,-.1))
 EM.SSMFSSN(y,xi=c(19,15,13), S=cov(y), la1=c(.3,1.2,1.3) , la2=c(3,-.4,-.1), nu1=5, nu2=2,
 family="MFSSTT" ,iter.max=500,tol=10^-9)

 # Example 5: Cost of living data (2016)
 data(cost)
 y <- cbind(cost$Cost.of.Living.Index,cost$Cappuccino.regular)
 EM.SSMFSSN(y, xi=colMeans(y), S=cov(y), la1=c(-1,2) , la2=c(-1,1), nu1=3, family="MFSST")
   }
