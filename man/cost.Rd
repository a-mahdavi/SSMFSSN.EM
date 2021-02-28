\name{cost}
\alias{cost}
\title{Cost of living data (2016)}
\usage{
data(cost)
}
\description{
Cost of living dataset
}
\examples{
 data(cost)
 y <- cbind(cost$Cost.of.Living.Index,cost$Cappuccino.regular)
 EM.SSMFSSN(y, xi=colMeans(y), S=cov(y), la1=c(-1,2) , la2=c(-1,1), nu1=3, family="MFSST")
}
