\name{r.SSMFSSN}
\alias{r.SSMFSSN}
\title{r.SSMFSSN function}
\usage{
r.SSMFSSN(n , xi, S, la1, la2, nu1=NULL, nu2=NULL, family="MFSSN")
}
\description{
 Generating random samples from multivariate SSMFSSN distributions 
 n: number of samples.
 xi: the vector of location parameter.
 S: the cov-variance matrix.
 la1 and la2: the vector of shape parameters.
 nu1 and nu2: the flatness parameters.
 family: distribution family to be used in fitting ("MFSSN", "MFSSTN","MFSSLSN",
 "MFSSCN", "MFSST", "MFSSTT").
}
\examples{
 # Example 1:
 # Simulating 100 samples from MFSSTN distribution:
 y <- r.SSMFSSN(n=100, xi=c(0,5), S=matrix(c(1,.4,.4,4),2,2), la1=c(-2,3),
 la2=c(.5,-.5), nu1=5, family="MFSSTN")

 # Example 2:
 # Simulating 100 samples from MFSSTT distribution:
 y <- r.SSMFSSN(n=100, xi=c(0,5), S=matrix(c(1,0.4,0.4,4),2,2), la1=c(-2,3),
 la2=c(0.5,-0.5), nu1=5, nu2=10, family="MFSSTT")

# Example 3:
 # Simulating 100 samples from MFSSCNe distribution:
 y <- r.SSMFSSN(n=100, xi=c(0,5), S=matrix(c(1,0.4,0.4,4),2,2), la1=c(-2,3),
 la2=c(0.5,-0.5), nu1=0.3, nu2=0.3, family="MFSSTT")
   }
