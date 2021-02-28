r.SSMFSSN <- function(n , xi, S, la1, la2, nu1=NULL, nu2=NULL, family="MFSSN"){
	if ((family != "MFSSN") && (family != "MFSSTN") && (family != 
        "MFSST") && (family != "MFSSSLN") && (family != "MFSSCN") && 
        (family != "MFSSTT"))
	  stop(paste("Family", family, "not recognized.\n", sep = " "))
		p <- length(xi)
	if( p!=ncol(S) || p!=length(la1) || p!=length(la2) )
	 stop("The size of the parametrs vectors are not compatibles.\n")
		r.SSN <- function(n , xi, S, la1, la2){
		require(MASS)
		sqrt.mt = function(S)	{
		p = ncol(S)
		if(p == 1) S.sqrt = as.matrix(sqrt(S))
		else{
		eig.S = eigen(S)
		S.sqrt = eig.S$ve %*% diag(sqrt(eig.S$va)) %*% t(eig.S$ve)}
		}
		p <- length(xi)
		y<-matrix(0,n,p) ; dec <- sqrt.mt(as.matrix(S))
		for(i in 1:n){
		x0 <- mvrnorm(1,mu=rep(0,p), Sigma=diag(1,p))
		x1 <- rnorm(1)
		if(x1<(t(la1)%*%x0+ t(la2)%*%x0^3)){
		y[i,] <- x0
		y[i,] <- xi+dec%*%y[i,]
			}
		else{
		y[i,] <- -x0
		y[i,] <- xi+dec%*%y[i,]
			}
				}
		return(y)
		}
		r.SSTN <- function(n , xi, S, la1, la2, nu){
		require(MASS)
		if (!is.numeric(nu))
		stop(paste("The nu1 parameter is not specified.\n"))
		sqrt.mt = function(S)	{
		p = ncol(S)
		if(p == 1) S.sqrt = as.matrix(sqrt(S))
		else{
		eig.S = eigen(S)
		S.sqrt = eig.S$ve %*% diag(sqrt(eig.S$va)) %*% t(eig.S$ve)}
		}
		p <- length(xi)
		y<-matrix(0,n,p) ; dec <-sqrt.mt(as.matrix(S))
		for(i in 1:n){
		x0 <- mvrnorm(1,mu=rep(0,p), Sigma=diag(1,p))
		x0 <- rgamma(1,nu/2,nu/2)^(-1/2)*x0
		x1 <- rnorm(1)
		if(x1<(t(la1)%*%x0+ t(la2)%*%x0^3)){
		y[i,] <- x0
		y[i,] <- xi+dec%*%y[i,]
			}
		else{
		y[i,] <- -x0
		y[i,] <- xi+dec%*%y[i,]
			}
				}
		return(y)
		}
		r.SSSLN <- function(n , xi, S, la1, la2, nu){
		require(MASS)
		if (!is.numeric(nu))
		stop(paste("The nu1 parameter is not specified.\n"))
		sqrt.mt = function(S)	{
		p = ncol(S)
		if(p == 1) S.sqrt = as.matrix(sqrt(S))
		else{
		eig.S = eigen(S)
		S.sqrt = eig.S$ve %*% diag(sqrt(eig.S$va)) %*% t(eig.S$ve)}
			}
		p <- length(xi)
		y<-matrix(0,n,p) ; dec <- sqrt.mt(as.matrix(S))
		for(i in 1:n){
		x0 <- mvrnorm(1,mu=rep(0,p), Sigma=diag(1,p))
		x0 <- rbeta(1,nu,1)^(-1/2)*x0
		x1 <- rnorm(1)
		if(x1<(t(la1)%*%x0+ t(la2)%*%x0^3)){
		y[i,] <- x0
		y[i,] <- xi+dec%*%y[i,]
			}
		else{
		y[i,] <- -x0
		y[i,] <- xi+dec%*%y[i,]
			}
				}
		return(y)
		}
		r.SSCN <- function(n , xi, S, la1, la2, nu1, nu2){
		require(MASS)
		if (!is.numeric(nu1)||!is.numeric(nu2))
		stop(paste("The nu1 and nu2 parameters are not specified.\n"))
		sqrt.mt = function(S)	{
		p = ncol(S)
		if(p == 1) S.sqrt = as.matrix(sqrt(S))
		else{
		eig.S = eigen(S)
		S.sqrt = eig.S$ve %*% diag(sqrt(eig.S$va)) %*% t(eig.S$ve)}
		}
		p <- length(xi)
		y<-matrix(0,n,p) ; dec <- sqrt.mt(as.matrix(S))
		for(i in 1:n){
		x0 <- mvrnorm(1,mu=rep(0,p), Sigma=diag(1,p))
		x0 <- sample(c(nu2,1),size=1,prob=c(nu1,(1-nu1)))^(-1/2)*x0
		x1 <- rnorm(1)
		if(x1<(t(la1)%*%x0+ t(la2)%*%x0^3)){
		y[i,] <- x0
		y[i,] <- xi+dec%*%y[i,]
			}
		else{
		y[i,] <- -x0
		y[i,] <- xi+dec%*%y[i,]
			}
				}
		return(y)
		}
		r.SST <- function(n , xi, S, la1, la2, nu){
		require(MASS)
		if (!is.numeric(nu))
		stop(paste("The nu1 parameter is not specified.\n"))
		sqrt.mt = function(S)	{
		p = ncol(S)
		if(p == 1) S.sqrt = as.matrix(sqrt(S))
		else{
		eig.S = eigen(S)
		S.sqrt = eig.S$ve %*% diag(sqrt(eig.S$va)) %*% t(eig.S$ve)}
		}
		p <- length(xi)
		y<-matrix(0,n,p) ; dec <- sqrt.mt(as.matrix(S))
		for(i in 1:n){
		w <- rgamma(1,nu/2,nu/2)^(-1/2)
		x0 <- mvrnorm(1,mu=rep(0,p), Sigma=diag(1,p))*w
		x1 <- rnorm(1)*w
		if(x1<(t(la1)%*%x0+ t(la2)%*%x0^3)){
		y[i,] <- x0
		y[i,] <- xi+dec%*%y[i,]
			}
		else{
		y[i,] <- -x0
		y[i,] <- xi+dec%*%y[i,]
			}
				}
		return(y)
		}
		r.SSTT <- function(n , xi, S, la1, la2, nu1, nu2){
		require(MASS)
		if (!is.numeric(nu1) || !is.numeric(nu2))
		stop(paste("The nu1 and nu2 parameters are not specified.\n"))
		sqrt.mt = function(S)	{
		p = ncol(S)
		if(p == 1) S.sqrt = as.matrix(sqrt(S))
		else{
		eig.S = eigen(S)
		S.sqrt = eig.S$ve %*% diag(sqrt(eig.S$va)) %*% t(eig.S$ve)}
		}
		p <- length(xi)
		y<-matrix(0,n,p) ; dec <- sqrt.mt(as.matrix(S))
		for(i in 1:n){
		x0 <- mvrnorm(1,mu=rep(0,p), Sigma=diag(1,p))
		x0 <- rgamma(1,nu1/2,nu1/2)^(-1/2)*x0
		x1 <- rnorm(1) ; x1 <- x1*rgamma(1,nu2/2,nu2/2)^(-1/2)
		if(x1<(t(la1)%*%x0+ t(la2)%*%x0^3)){
		y[i,] <- x0
		y[i,] <- xi+dec%*%y[i,]
			}
		else{
		y[i,] <- -x0
		y[i,] <- xi+dec%*%y[i,]
			}
				}
		return(y)
		}
	if (family=='MFSSN')
	y <- r.SSN(n=n,xi=xi, S=S, la1=la1 , la2=la2)
	if (family=='MFSSTN')
	y <- r.SSTN(n=n,xi=xi, S=S, la1=la1 , la2=la2, nu= nu1 )
	if (family=='MFSSSLN')
	y <- r.SSSLN(n=n,xi=xi, S=S, la1=la1 , la2=la2, nu= nu1)
	if (family=='MFSSCN')
	y <- r.SSCN(n=n,xi=xi, S=S, la1=la1 , la2=la2, nu1= nu1, nu2=nu2 )
	if (family=='MFSST')
	y <- r.SST(n=n,xi=xi, S=S, la1=la1 , la2=la2, nu= nu1)
	if (family=='MFSSTT')
	y <- r.SSTT(n=n,xi=xi, S=S, la1=la1 , la2=la2, nu1= nu1, nu2=nu2 )
	return(y)
	}









