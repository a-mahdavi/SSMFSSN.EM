EM.SSMFSSN <- function(y, xi=NULL, S=NULL, la1=NULL, la2=NULL, nu1=0.5 , nu2=0.5, family="MFSSN", get.init = FALSE, iter.max=100, tol=10^-6, equal=FALSE, CML=TRUE){
  	 	if ((family != "MFSSN") && (family != "MFSSTN") && (family != 
        "MFSST") && (family != "MFSSSLN") && (family != "MFSSCN") && 
        (family != "MFSSTT"))
	  stop(paste("Family", family, "not recognized.\n", sep = " "))
	y <- as.matrix(y) ; p <- ncol(y)
	if ( p <= 1 ) 
        stop("This function is developed for multivariate case.\n For univariate cas, see Mahdavi et al. (2021),
	Maximum likelihood estimation for scale-shape mixtures of flexible generalized skew normal 
	distributions via selection representation. Comput Stat (2021). https://doi.org/10.1007/s00180-021-01079-2")
	skewness <- function(y)
		sign(apply((y - matrix(rep(colMeans(y), nrow(y)), nrow = nrow(y), ncol = dim(y)[2], byrow = TRUE))^3, 2, sum))
	if (get.init) {
	xi <- colMeans(y) ; S=cov(y) ; la1=skewness(y) ; la2=rep(0,p)  }
 	if( p!=ncol(S) || p!=length(la1) || p!=length(la2) || p!=length(xi))
	 stop("The size of the initial values are not compatibles.\n")
	dmvt <- function (x, delta = rep(0, p), sigma = diag(p), df = 1, log = TRUE, 
  	  type = "shifted", checkSymmetry = TRUE) 
	{
   	 if (is.vector(x)) 
        x <- matrix(x, ncol = length(x))
   	 p <- ncol(x)
   	 if (df == 0 || is.null(df)) 
        return(dmvnorm(x, mean = delta, sigma = sigma, log = log))
    		if (!missing(delta)) {
        if (!is.null(dim(delta))) 
            dim(delta) <- NULL
        if (length(delta) != p) 
            stop("delta and sigma have non-conforming size")
   		 }
   	 if (!missing(sigma)) {
       	 if (p != ncol(sigma)) 
            stop("x and sigma have non-conforming size")
        	if (checkSymmetry && !isSymmetric(sigma, tol = sqrt(.Machine$double.eps), 
            check.attributes = FALSE)) 
            stop("sigma must be a symmetric matrix")
   		 }
    		type <- match.arg(type)
    		dec <- tryCatch(chol(sigma), error = function(e) e)
    		if (inherits(dec, "error")) {
      	  x.is.d <- colSums(t(x) != delta) == 0
      	  logretval <- rep.int(-Inf, nrow(x))
       	 logretval[x.is.d] <- Inf
   			 }
    			else {
      		  R.x_m <- backsolve(dec, t(x) - delta, transpose = TRUE)
     	 		  rss <- colSums(R.x_m^2)
       		 logretval <- lgamma((p + df)/2) - (lgamma(df/2) + sum(log(diag(dec))) + 
          		  p/2 * log(pi * df)) - 0.5 * (df + p) * log1p(rss/df)
   				 }
  			  names(logretval) <- rownames(x)
   			 if (log) 
      		  logretval
   			 else exp(logretval)
				}
	dmvnorm <-	function (x, mean = rep(0, p), sigma = diag(p), log = FALSE, 
 	   checkSymmetry = TRUE) 
		{
   	 if (is.vector(x)) 
   	     x <- matrix(x, ncol = length(x))
   		 p <- ncol(x)
   		 if (!missing(mean)) {
     	   if (!is.null(dim(mean))) 
            dim(mean) <- NULL
      	  if (length(mean) != p) 
            stop("x and mean have non-conforming size")
   		 }
   		 if (!missing(sigma)) {
   	     if (p != ncol(sigma)) 
            stop("x and sigma have non-conforming size")
       	 if (checkSymmetry && !isSymmetric(sigma, tol = sqrt(.Machine$double.eps), 
            check.attributes = FALSE)) 
            stop("sigma must be a symmetric matrix")
   		 }
   		 dec <- tryCatch(chol(sigma), error = function(e) e)
   		 if (inherits(dec, "error")) {
        x.is.mu <- colSums(t(x) != mean) == 0
        logretval <- rep.int(-Inf, nrow(x))
        logretval[x.is.mu] <- Inf
   		 }
  		  else {
        tmp <- backsolve(dec, t(x) - mean, transpose = TRUE)
        rss <- colSums(tmp^2)
        logretval <- -sum(log(diag(dec))) - 0.5 * p * log(2 * 
            pi) - 0.5 * rss
   		 }
   	 names(logretval) <- rownames(x)
   	 if (log) 
        logretval
  	  else exp(logretval)
	}
	EM.SSN <- function(y, xi, S, la1, la2, iter.max=100, tol=10^-6, CML=TRUE){ 
   	 begin <- proc.time()[3] 
	sqrt.mt = function(S,inverse=T)	{
	p = ncol(S)
	if(p == 1) S.sqrt = as.matrix(sqrt(S))
	else{
	eig.S = eigen(S)
	if(inverse)
	S.sqrt = eig.S$ve %*% diag(1/sqrt(eig.S$va)) %*% t(eig.S$ve)
	else 	S.sqrt = eig.S$ve %*% diag(sqrt(eig.S$va)) %*% t(eig.S$ve)
 			}		}
	  y <- as.matrix(y) ; n <- dim(y)[1] ; p <- dim(y)[2]
	  S <- as.matrix(S) ; xi.mat <- matrix(xi,n,p,byrow=T)
        dif <- 1
        count <- 0
	solve.S <- solve(S) ; S.s <- sqrt.mt(S)
	y.xi <- y-xi.mat ; one <- matrix(1,p,1)
	del <- matrix(t(la1)%*%S.s,p,1)  ; Del <- diag(abs(la2)^(1/3)*sign(la2))%*%S.s
	del.y <-  t(del)%*%t(y.xi)
	Del.y <-  matrix(apply(y.xi,1,function(x) Del%*%x),p,n)
	LL <- 1
	while ((dif > tol) && (count <= iter.max)) {
	Q <- apply(y.xi,1,function(x) t(x)%*%solve.S%*%(x))
	s1 <- s2 <- rep(1,n)
	s3 <- as.numeric(del.y + t(one)%*%Del.y^3 + dnorm(del.y + t(one)%*%Del.y^3)/pnorm(del.y + t(one)%*%Del.y^3))
	if (CML){
	xi <- optim(xi,function(x){
			xi.mat <- matrix(x,n,p,byrow=T)
			y.xi <- y-xi.mat
			del.y <-  t(del)%*%t(y.xi)
			Del.y <-  matrix(apply(y.xi,1,function(x) Del%*%x),p,n)
			-sum(log(2*dmvnorm(y,x,S)*pnorm(del.y + t(one)%*%Del.y^3))) 
			},method="BFGS")$par
	} else {
	xi <-  optim(xi,function(x){
	xi.mat <- matrix(x,n,p,byrow=T)
	y.xi <- y-xi.mat
	del.y <-  t(del)%*%t(y.xi)
	Del.y <-  matrix(apply(y.xi,1,function(x) Del%*%x),p,n)
	Sum <- 0 ; 	one.D <- as.numeric(t(one)%*%Del.y^3)
	for(i in 1:n)
	Sum <- Sum + s2[i]*(del.y[i]*del+3*as.numeric(one.D[i]+del.y[i])*t(Del)%*%Del.y[,i]^2+del*one.D[i])-
			s3[i]*(del+3*t(Del)%*%Del.y[,i]^2)
		norm(solve.S%*%colSums(s1*y.xi)+Sum)
		},method="BFGS")$par 
		}
	xi.mat <- matrix(xi,n,p,byrow=T)
	y.xi <- y-xi.mat ; 	one.D <- as.numeric(t(one)%*%Del.y^3)
	S <- 1/n*matrix(rowSums(matrix(s1,p^2,n,byrow=T)*apply(y.xi,1,function(x) x%*%t(x))),p,p)
	solve.S <- solve(S) ;  S.s <- sqrt.mt(S)
	A <- matrix(rowSums(matrix(s2,p^2,n,byrow=T)*apply(y.xi,1,function(x) x%*%t(x))),p,p)
	del <- matrix(solve(A)%*%colSums(s3*y.xi-s2*one.D*y.xi),p,1)
	del.y <-  t(del)%*%t(y.xi)
	if (CML){
	Del <- optim(as.vector(Del),function(x){
		Del <- matrix (x , p, p)
		Del.y <-  matrix(apply(y.xi,1,function(x) Del%*%x),p,n)
		-sum(log(2*dmvnorm(y,xi,S)*pnorm(del.y + t(one)%*%Del.y^3))) 
		},method="BFGS")$par
		} else {
	Del <-  optim(as.vector(Del),function(x){
	Del <- matrix (x , p, p)
	Del.y <-  matrix(apply(y.xi,1,function(x) Del%*%x),p,n)
	one.D <- as.numeric(t(one)%*%Del.y^3)
		eq <- 0
		for(i in 1:n)
		eq <- eq + ( s3[i]-s2[i]*(del.y[i]+one.D[i]) )*Del.y[,i]^2%*%t(y.xi[i,])
		norm( eq )
		},method="BFGS")$par
		}
	Del <- matrix(Del,p,p)
	Del.y <-  apply(y.xi,1,function(x) Del%*%x)
	LL.new <-  sum(log(2*dmvnorm(y,xi,S)*pnorm(del.y + t(one)%*%Del.y^3))) # log-likelihood function
	count <- count +1 
	dif <- abs(LL.new/LL-1)
	LL <- LL.new
	}
	la1 <- sqrt.mt(S,inverse=F)%*%del ; la2=(Del%*%(sqrt.mt(S,inverse=F)))^3%*%one
	aic <- -2 * LL.new + 2 * (3*p+p*(p+1)/2)
	bic <- -2 * LL.new + log(n) * (3*p+p*(p+1)/2)
	end <- proc.time()[3]
	time <- end-begin
	list(xi=xi, S=S, la1=la1 , la2=la2, loglik=LL.new, AIC=aic, BIC=bic, iter=count, elapsed=as.numeric(time))
	}
	EM.SSTN <- function(y, xi, S, la1, la2, nu, iter.max=100, tol=10^-6, CML=TRUE){  
 	begin <- proc.time()[3] 
	sqrt.mt = function(S,inverse=T)	{
	p = ncol(S)
	if(p == 1) S.sqrt = as.matrix(sqrt(S))
	else{
	eig.S = eigen(S)
	if(inverse)
	S.sqrt = eig.S$ve %*% diag(1/sqrt(eig.S$va)) %*% t(eig.S$ve)
	else 	S.sqrt = eig.S$ve %*% diag(sqrt(eig.S$va)) %*% t(eig.S$ve)
 			}		}
	  y <- as.matrix(y) ; n <- dim(y)[1] ; p <- dim(y)[2]
	  S <- as.matrix(S) ; xi.mat <- matrix(xi,n,p,byrow=T)
        dif <- 1
        count <- 0
	solve.S <- solve(S) ; S.s <- sqrt.mt(S)
	y.xi <- y-xi.mat ; one <- matrix(1,p,1)
	del <- matrix(t(la1)%*%S.s,p,1)  ; Del <- diag(abs(la2)^(1/3)*sign(la2))%*%S.s
	del.y <-  t(del)%*%t(y.xi)
	Del.y <-  matrix(apply(y.xi,1,function(x) Del%*%x),p,n)
	LL <- 1
	while ((dif > tol) && (count <= iter.max)) {
	Q <- apply(y.xi,1,function(x) t(x)%*%solve.S%*%(x))
	s1 <- (nu+p)/(nu+Q) ; s2 <- rep(1,n)
	s3 <- as.numeric(del.y + t(one)%*%Del.y^3 + dnorm(del.y + t(one)%*%Del.y^3)/pnorm(del.y + t(one)%*%Del.y^3))
	if (CML){
	xi <- optim(xi,function(x){
			xi.mat <- matrix(x,n,p,byrow=T)
			y.xi <- y-xi.mat
			del.y <-  t(del)%*%t(y.xi)
			Del.y <-  matrix(apply(y.xi,1,function(x) Del%*%x),p,n)
			-sum(log(2*dmvt(y,x,S,nu,log=F)*pnorm(del.y + t(one)%*%Del.y^3))) 
			},method="BFGS")$par
	} else {
	xi <-  optim(xi,function(x){
	xi.mat <- matrix(x,n,p,byrow=T)
	y.xi <- y-xi.mat
	del.y <-  t(del)%*%t(y.xi)
	Del.y <-  matrix(apply(y.xi,1,function(x) Del%*%x),p,n)
	Sum <- 0 ; 	one.D <- as.numeric(t(one)%*%Del.y^3)
	for(i in 1:n)
	Sum <- Sum + s2[i]*(del.y[i]*del+3*as.numeric(one.D[i]+del.y[i])*t(Del)%*%Del.y[,i]^2+del*one.D[i])-
			s3[i]*(del+3*t(Del)%*%Del.y[,i]^2)
		norm(solve.S%*%colSums(s1*y.xi)+Sum)
		},method="BFGS")$par 
		}
	xi.mat <- matrix(xi,n,p,byrow=T)
	y.xi <- y-xi.mat ; 	one.D <- as.numeric(t(one)%*%Del.y^3)
	S <- 1/n*matrix(rowSums(matrix(s1,p^2,n,byrow=T)*apply(y.xi,1,function(x) x%*%t(x))),p,p)
	solve.S <- solve(S) ;  S.s <- sqrt.mt(S) 
	A <- matrix(rowSums(matrix(s2,p^2,n,byrow=T)*apply(y.xi,1,function(x) x%*%t(x))),p,p)
	del <- matrix(solve(A)%*%colSums(s3*y.xi-s2*one.D*y.xi),p,1)
	del.y <-  t(del)%*%t(y.xi)
	if (CML){
	Del <- optim(as.vector(Del),function(x){
		Del <- matrix (x , p, p)
		Del.y <-  matrix(apply(y.xi,1,function(x) Del%*%x),p,n)
		-sum(log(2*dmvt(y,xi,S,nu,log=F)*pnorm(del.y + t(one)%*%Del.y^3))) 
		},method="BFGS")$par
		} else {
	Del <-  optim(as.vector(Del),function(x){
	Del <- matrix (x , p, p)
	Del.y <-  matrix(apply(y.xi,1,function(x) Del%*%x),p,n)
	one.D <- as.numeric(t(one)%*%Del.y^3)
		eq <- 0
		for(i in 1:n)
		eq <- eq + ( s3[i]-s2[i]*(del.y[i]+one.D[i]) )*Del.y[,i]^2%*%t(y.xi[i,])
		norm( eq )
		},method="BFGS")$par
		}
	Del <- matrix(Del,p,p)
	Del.y <-  apply(y.xi,1,function(x) Del%*%x)
	nu <- optim(nu,function(x){
		-sum(log(2*dmvt(y,xi,S,x,log=F)*pnorm(del.y + t(one)%*%Del.y^3))) 
		},method="L-BFGS-B",lower=.01,upper=100)$par
	LL.new <-  sum(log(2*dmvt(y,xi,S,nu,log=F)*pnorm(del.y + t(one)%*%Del.y^3))) # log-likelihood function
	count <- count +1 
	dif <- abs(LL.new/LL-1)
	LL <- LL.new
	}
	la1 <- sqrt.mt(S,inverse=F)%*%del ; la2=(Del%*%(sqrt.mt(S,inverse=F)))^3%*%one
	aic <- -2 * LL.new + 2 * (3*p+p*(p+1)/2+1)
	bic <- -2 * LL.new + log(n) * (3*p+p*(p+1)/2+1)
	end <- proc.time()[3]
	time <- end-begin
	list(xi=xi, S=S, la1=la1 , la2=la2, nu=nu, loglik=LL.new, AIC=aic, BIC=bic, iter=count,elapsed=as.numeric(time))
	}
	EM.SSSLN <- function(y, xi, S, la1, la2, nu, iter.max=100, tol=10^-6, CML=TRUE){  
	  begin <- proc.time()[3] 
	dmvs <- function(x,xi,S,nu){ # density of multivariate slash distribution
	p <- dim(x)[2] ; n <- dim(x)[1] ; x <- as.matrix(x); xi.mat <- matrix(xi,n,p,byrow=T)
	x.xi <- x-xi.mat ; S <- matrix(S,p,p) ; solve.S <- solve(S)
	Q <- apply(x.xi,1,function(x) t(x)%*%solve.S%*%(x))
	2^nu*nu*gamma(nu+p/2)/((pi)^(p/2)*det(S)^(1/2)*Q^(nu+p/2))*pgamma(Q/2,nu+p/2)
		}
	sqrt.mt = function(S,inverse=T)	{
	p = ncol(S)
	if(p == 1) S.sqrt = as.matrix(sqrt(S))
	else{
	eig.S = eigen(S)
	if(inverse)
	S.sqrt = eig.S$ve %*% diag(1/sqrt(eig.S$va)) %*% t(eig.S$ve)
	else 	S.sqrt = eig.S$ve %*% diag(sqrt(eig.S$va)) %*% t(eig.S$ve)
 			}		}	  
	y <- as.matrix(y) ; n <- dim(y)[1] ; p <- dim(y)[2]
	  S <- as.matrix(S) ; xi.mat <- matrix(xi,n,p,byrow=T)
        dif <- 1
        count <- 0
	solve.S <- solve(S) ; S.s <- sqrt.mt(S)
	y.xi <- y-xi.mat ; one <- matrix(1,p,1)
	del <- matrix(t(la1)%*%S.s,p,1)  ; Del <- diag(abs(la2)^(1/3)*sign(la2))%*%S.s
	del.y <-  t(del)%*%t(y.xi)
	Del.y <-  matrix(apply(y.xi,1,function(x) Del%*%x),p,n)
	LL <- 1
	while ((dif > tol) && (count <= iter.max)) {
	Q <- apply(y.xi,1,function(x) t(x)%*%solve.S%*%(x))
	s1 <- (2*nu+p)/Q*pgamma(Q/2,nu+p/2+1)/pgamma(Q/2,nu+p/2) ; s2 <- rep(1,n)
	s3 <- as.numeric(del.y + t(one)%*%Del.y^3 + dnorm(del.y + t(one)%*%Del.y^3)/pnorm(del.y + t(one)%*%Del.y^3))
	if (CML){
	xi <- optim(xi,function(x){
			xi.mat <- matrix(x,n,p,byrow=T)
			y.xi <- y-xi.mat
			del.y <-  t(del)%*%t(y.xi)
			Del.y <-  matrix(apply(y.xi,1,function(x) Del%*%x),p,n)
			-sum(log(2*dmvs(y,x,S,nu)*pnorm(del.y + t(one)%*%Del.y^3))) 
			},method="BFGS")$par
	} else {
	xi <-  optim(xi,function(x){
	xi.mat <- matrix(x,n,p,byrow=T)
	y.xi <- y-xi.mat
	del.y <-  t(del)%*%t(y.xi)
	Del.y <-  matrix(apply(y.xi,1,function(x) Del%*%x),p,n)
	Sum <- 0 ; 	one.D <- as.numeric(t(one)%*%Del.y^3)
	for(i in 1:n)
	Sum <- Sum + s2[i]*(del.y[i]*del+3*as.numeric(one.D[i]+del.y[i])*t(Del)%*%Del.y[,i]^2+del*one.D[i])-
			s3[i]*(del+3*t(Del)%*%Del.y[,i]^2)
		norm(solve.S%*%colSums(s1*y.xi)+Sum)
		},method="BFGS")$par 
		}
	xi.mat <- matrix(xi,n,p,byrow=T)
	y.xi <- y-xi.mat ; 	one.D <- as.numeric(t(one)%*%Del.y^3)
	S <- 1/n*matrix(rowSums(matrix(s1,p^2,n,byrow=T)*apply(y.xi,1,function(x) x%*%t(x))),p,p)
	solve.S <- solve(S) ;  S.s <- sqrt.mt(S)
	A <- matrix(rowSums(matrix(s2,p^2,n,byrow=T)*apply(y.xi,1,function(x) x%*%t(x))),p,p)
	del <- matrix(solve(A)%*%colSums(s3*y.xi-s2*one.D*y.xi),p,1)
	del.y <-  t(del)%*%t(y.xi)
	if (CML){
	Del <- optim(as.vector(Del),function(x){
		Del <- matrix (x , p, p)
		Del.y <-  matrix(apply(y.xi,1,function(x) Del%*%x),p,n)
		-sum(log(2*dmvs(y,xi,S,nu)*pnorm(del.y + t(one)%*%Del.y^3))) 
		},method="BFGS")$par
		} else {
	Del <-  optim(as.vector(Del),function(x){
	Del <- matrix (x , p, p)
	Del.y <-  matrix(apply(y.xi,1,function(x) Del%*%x),p,n)
	one.D <- as.numeric(t(one)%*%Del.y^3)
		eq <- 0
		for(i in 1:n)
		eq <- eq + ( s3[i]-s2[i]*(del.y[i]+one.D[i]) )*Del.y[,i]^2%*%t(y.xi[i,])
		norm( eq )
		},method="BFGS")$par
		}
	Del <- matrix(Del,p,p)
	Del.y <-  apply(y.xi,1,function(x) Del%*%x)
	nu <- optim(nu,function(x){
		-sum(log(2*dmvs(y,xi,S,x)*pnorm(del.y + t(one)%*%Del.y^3))) 
		},method="L-BFGS-B",lower=.01,upper=100)$par
	LL.new <-  sum(log(2*dmvs(y,xi,S,nu)*pnorm(del.y + t(one)%*%Del.y^3))) # log-likelihood function
	count <- count +1 
	dif <- abs(LL.new/LL-1)
	LL <- LL.new
	}
	la1 <- sqrt.mt(S,inverse=F)%*%del ; la2=(Del%*%(sqrt.mt(S,inverse=F)))^3%*%one
	aic <- -2 * LL.new + 2 * (3*p+p*(p+1)/2+1)
	bic <- -2 * LL.new + log(n) * (3*p+p*(p+1)/2+1)
	end <- proc.time()[3]
	time <- end-begin
	list(xi=xi, S=S, la1=la1 , la2=la2, nu=nu, loglik=LL.new, AIC=aic, BIC=bic, iter=count,elapsed=as.numeric(time))
	}
	EM.SSCN <- function(y, xi, S, la1, la2, nu1, nu2, equal=FALSE, iter.max=100, tol=10^-6, CML=TRUE){  
	 begin <- proc.time()[3] 
	sqrt.mt = function(S,inverse=T)	{
	p = ncol(S)
	if(p == 1) S.sqrt = as.matrix(sqrt(S))
	else{
	eig.S = eigen(S)
	if(inverse)
	S.sqrt = eig.S$ve %*% diag(1/sqrt(eig.S$va)) %*% t(eig.S$ve)
	else 	S.sqrt = eig.S$ve %*% diag(sqrt(eig.S$va)) %*% t(eig.S$ve)
 			}		}
	  y <- as.matrix(y) ; n <- dim(y)[1] ; p <- dim(y)[2]
	  S <- as.matrix(S) ; xi.mat <- matrix(xi,n,p,byrow=T)
        dif <- 1
        count <- 0
	solve.S <- solve(S) ; S.s <- sqrt.mt(S)
	y.xi <- y-xi.mat ; one <- matrix(1,p,1)
	del <- matrix(t(la1)%*%S.s,p,1)  ; Del <- diag(abs(la2)^(1/3)*sign(la2))%*%S.s
	del.y <-  t(del)%*%t(y.xi)
	Del.y <-  matrix(apply(y.xi,1,function(x) Del%*%x),p,n)
	Q <- apply(y.xi,1,function(x) t(x)%*%solve.S%*%(x))
	xi.w <- del.y + t(one)%*%Del.y^3
	LL <- 1
	while ((dif > tol) && (count <= iter.max)) {
	s1 <- (1-nu1+nu1*nu2^(1+p/2)*exp((1-nu2)*(Q/2)))/(1-nu1+nu1*nu2^(p/2)*exp((1-nu2)*(Q/2)))
	s2 <- rep(1,n)
	s3 <- as.numeric(xi.w + dnorm(xi.w)/pnorm(xi.w))
	if (CML){
	xi <- optim(xi,function(x){
			xi.mat <- matrix(x,n,p,byrow=T)
			y.xi <- y-xi.mat
			del.y <-  t(del)%*%t(y.xi)
			Del.y <-  matrix(apply(y.xi,1,function(x) Del%*%x),p,n)
			xi.w <- del.y + t(one)%*%Del.y^3
			-sum(log(2*(nu1*dmvnorm(y,x,S/nu2)+(1-nu1)*dmvnorm(y,x,S))*pnorm(xi.w)))
			},method="BFGS")$par
	} else {
	xi <-  optim(xi,function(x){
	xi.mat <- matrix(x,n,p,byrow=T)
	y.xi <- y-xi.mat
	del.y <-  t(del)%*%t(y.xi)
	Del.y <-  matrix(apply(y.xi,1,function(x) Del%*%x),p,n)
	Sum <- 0 ; 	one.D <- as.numeric(t(one)%*%Del.y^3)
	for(i in 1:n)
	Sum <- Sum + s2[i]*(del.y[i]*del+3*as.numeric(one.D[i]+del.y[i])*t(Del)%*%Del.y[,i]^2+del*one.D[i])-
			s3[i]*(del+3*t(Del)%*%Del.y[,i]^2)
		norm(solve.S%*%colSums(s1*y.xi)+Sum)
		},method="BFGS")$par 
		}
	xi.mat <- matrix(xi,n,p,byrow=T)
	y.xi <- y-xi.mat ; 	one.D <- as.numeric(t(one)%*%Del.y^3)
	S <- 1/n*matrix(rowSums(matrix(s1,p^2,n,byrow=T)*apply(y.xi,1,function(x) x%*%t(x))),p,p)
	solve.S <- solve(S) ;  S.s <- sqrt.mt(S)
	A <- matrix(rowSums(matrix(s2,p^2,n,byrow=T)*apply(y.xi,1,function(x) x%*%t(x))),p,p)
	del <- matrix(solve(A)%*%colSums(s3*y.xi-s2*one.D*y.xi),p,1)
	del.y <-  t(del)%*%t(y.xi)
	if (CML){
	Del <- optim(as.vector(Del),function(x){
		Del <- matrix (x , p, p)
		Del.y <-  matrix(apply(y.xi,1,function(x) Del%*%x),p,n)
		-sum(log(2*(nu1*dmvnorm(y,xi,S/nu2)+(1-nu1)*dmvnorm(y,xi,S))*pnorm(del.y + t(one)%*%Del.y^3)))
		},method="BFGS")$par
		} else {
	Del <-  optim(as.vector(Del),function(x){
	Del <- matrix (x , p, p)
	Del.y <-  matrix(apply(y.xi,1,function(x) Del%*%x),p,n)
	one.D <- as.numeric(t(one)%*%Del.y^3)
		eq <- 0
		for(i in 1:n)
		eq <- eq + ( s3[i]-s2[i]*(del.y[i]+one.D[i]) )*Del.y[,i]^2%*%t(y.xi[i,])
		norm( eq )
		},method="BFGS")$par
		}
	Del <- matrix(Del,p,p)
	Del.y <-  apply(y.xi,1,function(x) Del%*%x)
	Q <- apply(y.xi,1,function(x) t(x)%*%solve.S%*%(x))
	xi.w <- del.y + t(one)%*%Del.y^3
	if (equal){	
	nu1 <- optim(nu1,function(x){
		 -sum(log(2*(x*dmvnorm(y,xi,S/x)+(1-x)*dmvnorm(y,xi,S))*pnorm(xi.w)))
		},method="L-BFGS-B",lower=.01,upper=1)$par
	nu2 <- nu1
	LL.new <-   sum(log(2*(nu1*dmvnorm(y,xi,S/nu2)+(1-nu1)*dmvnorm(y,xi,S))*pnorm(xi.w))) # log-likelihood function
	aic <- -2 * LL.new + 2 * (3*p+p*(p+1)/2+1)
	bic <- -2 * LL.new + log(n) * (3*p+p*(p+1)/2+1)
		} else{
	nu <- optim(c(nu1,nu2),function(x){
		 -sum(log(2*(x[1]*dmvnorm(y,xi,S/x[2])+(1-x[1])*dmvnorm(y,xi,S))*pnorm(xi.w)))
		},method="L-BFGS-B",lower=.01,upper=1)$par
	LL.new <-   sum(log(2*(nu1*dmvnorm(y,xi,S/nu2)+(1-nu1)*dmvnorm(y,xi,S))*pnorm(xi.w))) # log-likelihood function
	nu1 <- nu[1] ; nu2 <- nu[2]
	aic <- -2 * LL.new + 2 * (3*p+p*(p+1)/2+2)
	bic <- -2 * LL.new + log(n) * (3*p+p*(p+1)/2+2)
	}
	count <- count +1 
	dif <- abs(LL.new/LL-1)
	LL <- LL.new
	}
	la1 <- sqrt.mt(S,inverse=F)%*%del ; la2=(Del%*%(sqrt.mt(S,inverse=F)))^3%*%one
	end <- proc.time()[3]
	time <- end-begin
	list(xi=xi, S=S, la1=la1 , la2=la2, nu1=nu1, nu2=nu2, loglik=LL.new, AIC=aic, BIC=bic, iter=count, elapsed=as.numeric(time))
	}
	EM.SST <- function(y, xi, S, la1, la2, nu, iter.max=100, tol=10^-6, CML=TRUE){  
   	 begin <- proc.time()[3] 
	sqrt.mt = function(S,inverse=T)	{
	p = ncol(S)
	if(p == 1) S.sqrt = as.matrix(sqrt(S))
	else{
	eig.S = eigen(S)
	if(inverse)
	S.sqrt = eig.S$ve %*% diag(1/sqrt(eig.S$va)) %*% t(eig.S$ve)
	else 	S.sqrt = eig.S$ve %*% diag(sqrt(eig.S$va)) %*% t(eig.S$ve)
 			}		}
	  y <- as.matrix(y) ; n <- dim(y)[1] ; p <- dim(y)[2]
	  S <- as.matrix(S) ; xi.mat <- matrix(xi,n,p,byrow=T)
        dif <- 1
        count <- 0
	solve.S <- solve(S) ; S.s <- sqrt.mt(S)
	y.xi <- y-xi.mat ; one <- matrix(1,p,1)
	del <- matrix(t(la1)%*%S.s,p,1)  ; Del <- diag(abs(la2)^(1/3)*sign(la2))%*%S.s
	del.y <-  t(del)%*%t(y.xi)
	Del.y <-  matrix(apply(y.xi,1,function(x) Del%*%x),p,n)
	LL <- 1
	while ((dif > tol) && (count <= iter.max)) {
	Q <- apply(y.xi,1,function(x) t(x)%*%solve.S%*%(x))
	xi.w <- del.y + t(one)%*%Del.y^3 ; M <- xi.w*sqrt((nu+p)/(Q+nu))
	s1 <-	s2 <- as.numeric(((nu+p)/(Q+nu))*pt(M*sqrt((nu+p+2)/(nu+p)), nu+p+2)/pt(M, nu+p))
	s3 <- as.numeric(1/pt(M, nu+p)*(M*sqrt((nu+p)/(nu+Q))*pt(M*sqrt((nu+p+2)/(nu+p)), nu+p+2)+
		(Q+nu)^(-1/2)*gamma((nu+p+1)/2)/sqrt(pi)/gamma((nu+p)/2)*(1+xi.w^2/(Q+nu))^(-(nu+p+1)/2)))
	if (CML){
	xi <- optim(xi,function(x){
			xi.mat <- matrix(x,n,p,byrow=T)
			y.xi <- y-xi.mat
			del.y <-  t(del)%*%t(y.xi)
			Del.y <-  matrix(apply(y.xi,1,function(x) Del%*%x),p,n)
			xi.w <- del.y + t(one)%*%Del.y^3
			M <- xi.w*sqrt((nu+p)/(Q+nu))
			-sum(log(2*dmvt(y,x,S,nu,log=F)*pt(M,nu+p))) 
			},method="BFGS")$par
	} else {
	xi <-  optim(xi,function(x){
	xi.mat <- matrix(x,n,p,byrow=T)
	y.xi <- y-xi.mat
	del.y <-  t(del)%*%t(y.xi)
	Del.y <-  matrix(apply(y.xi,1,function(x) Del%*%x),p,n)
	Sum <- 0 ; 	one.D <- as.numeric(t(one)%*%Del.y^3)
	for(i in 1:n)
	Sum <- Sum + s2[i]*(del.y[i]*del+3*as.numeric(one.D[i]+del.y[i])*t(Del)%*%Del.y[,i]^2+del*one.D[i])-
			s3[i]*(del+3*t(Del)%*%Del.y[,i]^2)
		norm(solve.S%*%colSums(s1*y.xi)+Sum)
		},method="BFGS")$par 
		}
	xi.mat <- matrix(xi,n,p,byrow=T)
	y.xi <- y-xi.mat ; 	one.D <- as.numeric(t(one)%*%Del.y^3)
	S <- 1/n*matrix(rowSums(matrix(s1,p^2,n,byrow=T)*apply(y.xi,1,function(x) x%*%t(x))),p,p)
	solve.S <- solve(S) ;  S.s <- sqrt.mt(S)
	A <- matrix(rowSums(matrix(s2,p^2,n,byrow=T)*apply(y.xi,1,function(x) x%*%t(x))),p,p)
	del <- matrix(solve(A)%*%colSums(s3*y.xi-s2*one.D*y.xi),p,1)
	del.y <-  t(del)%*%t(y.xi)
	if (CML){
	Del <- optim(as.vector(Del),function(x){
		Del <- matrix (x , p, p)
		Del.y <-  matrix(apply(y.xi,1,function(x) Del%*%x),p,n)
		-sum(log(2*dmvt(y,xi,S,nu,log=F)*pt((del.y + t(one)%*%Del.y^3)*sqrt((nu+p)/(Q+nu)), nu+p))) 
		},method="BFGS")$par
		} else {
	Del <-  optim(as.vector(Del),function(x){
	Del <- matrix (x , p, p)
	Del.y <-  matrix(apply(y.xi,1,function(x) Del%*%x),p,n)
	one.D <- as.numeric(t(one)%*%Del.y^3)
		eq <- 0
		for(i in 1:n)
		eq <- eq + ( s3[i]-s2[i]*(del.y[i]+one.D[i]) )*Del.y[,i]^2%*%t(y.xi[i,])
		norm( eq )
		},method="BFGS")$par
		}
	Del <- matrix(Del,p,p)
	Del.y <-  apply(y.xi,1,function(x) Del%*%x)
	nu <- optim(nu,function(x){
		-sum(log(2*dmvt(y,xi,S,x,log=F)*pt((del.y + t(one)%*%Del.y^3)*sqrt((x+p)/(Q+x)), x+p))) 
		},method="L-BFGS-B",lower=.01,upper=100)$par
	LL.new <-  sum(log(2*dmvt(y,xi,S,nu,log=F)*pt((del.y + t(one)%*%Del.y^3)*sqrt((nu+p)/(Q+nu)), nu+p))) # log-likelihood function
	count <- count +1 
	dif <- abs(LL.new/LL-1)
	LL <- LL.new
	}
	la1 <- sqrt.mt(S,inverse=F)%*%del ; la2=(Del%*%(sqrt.mt(S,inverse=F)))^3%*%one
	aic <- -2 * LL.new + 2 * (3*p+p*(p+1)/2+1)
	bic <- -2 * LL.new + log(n) * (3*p+p*(p+1)/2+1)
	end <- proc.time()[3]
	time <- end-begin
	list(xi=xi, S=S, la1=la1 , la2=la2, nu=nu, loglik=LL.new, AIC=aic, BIC=bic, iter=count,elapsed=as.numeric(time))
	}
	EM.SSTT <- function(y, xi, S, la1, la2, nu1, nu2, iter.max=100, tol=10^-6, CML=TRUE, equal=FALSE){  
	  begin <- proc.time()[3] 
	sqrt.mt = function(S,inverse=T)	{
	p = ncol(S)
	if(p == 1) S.sqrt = as.matrix(sqrt(S))
	else{
	eig.S = eigen(S)
	if(inverse)
	S.sqrt = eig.S$ve %*% diag(1/sqrt(eig.S$va)) %*% t(eig.S$ve)
	else 	S.sqrt = eig.S$ve %*% diag(sqrt(eig.S$va)) %*% t(eig.S$ve)
 			}		}
	  y <- as.matrix(y) ; n <- dim(y)[1] ; p <- dim(y)[2]
	  S <- as.matrix(S) ; xi.mat <- matrix(xi,n,p,byrow=T)
        dif <- 1
        count <- 0
	solve.S <- solve(S) ; S.s <- sqrt.mt(S)
	y.xi <- y-xi.mat ; one <- matrix(1,p,1)
	del <- matrix(t(la1)%*%S.s,p,1)  ; Del <- diag(abs(la2)^(1/3)*sign(la2))%*%S.s
	del.y <-  t(del)%*%t(y.xi)
	Del.y <-  matrix(apply(y.xi,1,function(x) Del%*%x),p,n)
	LL <- 1
	while ((dif > tol) && (count <= iter.max)) {
	Q <- apply(y.xi,1,function(x) t(x)%*%solve.S%*%(x))
	xi.w <- del.y + t(one)%*%Del.y^3
	s1 <- (nu1+p)/(nu1+Q) 
	s2 <- as.numeric(pt(xi.w*sqrt((nu2+1)/nu2), nu2+2)/pt(xi.w, nu2))
	s3 <- as.numeric(1/pt(xi.w, nu2)*(pt(xi.w*sqrt((nu2+1)/nu2), nu2+2)*xi.w+dt(xi.w, nu2)))
	if (CML){
	xi <- optim(xi,function(x){
			xi.mat <- matrix(x,n,p,byrow=T)
			y.xi <- y-xi.mat
			del.y <-  t(del)%*%t(y.xi)
			Del.y <-  matrix(apply(y.xi,1,function(x) Del%*%x),p,n)
			xi.w <- del.y + t(one)%*%Del.y^3
			-sum(log(2*dmvt(y,x,S,nu1,log=F)*pt(xi.w,nu2))) 
			},method="BFGS")$par
	} else {
	xi <-  optim(xi,function(x){
	xi.mat <- matrix(x,n,p,byrow=T)
	y.xi <- y-xi.mat
	del.y <-  t(del)%*%t(y.xi)
	Del.y <-  matrix(apply(y.xi,1,function(x) Del%*%x),p,n)
	Sum <- 0 ; 	one.D <- as.numeric(t(one)%*%Del.y^3)
	for(i in 1:n)
	Sum <- Sum + s2[i]*(del.y[i]*del+3*as.numeric(one.D[i]+del.y[i])*t(Del)%*%Del.y[,i]^2+del*one.D[i])-
			s3[i]*(del+3*t(Del)%*%Del.y[,i]^2)
		norm(solve.S%*%colSums(s1*y.xi)+Sum)
		},method="BFGS")$par 
		}
	xi.mat <- matrix(xi,n,p,byrow=T)
	y.xi <- y-xi.mat ; 	one.D <- as.numeric(t(one)%*%Del.y^3)
	S <- 1/n*matrix(rowSums(matrix(s1,p^2,n,byrow=T)*apply(y.xi,1,function(x) x%*%t(x))),p,p)
	solve.S <- solve(S) ;  S.s <- sqrt.mt(S)
	A <- matrix(rowSums(matrix(s2,p^2,n,byrow=T)*apply(y.xi,1,function(x) x%*%t(x))),p,p)
	del <- matrix(solve(A)%*%colSums(s3*y.xi-s2*one.D*y.xi),p,1)
	del.y <-  t(del)%*%t(y.xi)
	if (CML){
	Del <- optim(as.vector(Del),function(x){
		Del <- matrix (x , p, p)
		Del.y <-  matrix(apply(y.xi,1,function(x) Del%*%x),p,n)
		-sum(log(2*dmvt(y,xi,S,nu1,log=F)*pt(del.y + t(one)%*%Del.y^3, nu2))) 
		},method="BFGS")$par
		} else {
	Del <-  optim(as.vector(Del),function(x){
	Del <- matrix (x , p, p)
	Del.y <-  matrix(apply(y.xi,1,function(x) Del%*%x),p,n)
	one.D <- as.numeric(t(one)%*%Del.y^3)
		eq <- 0
		for(i in 1:n)
		eq <- eq + ( s3[i]-s2[i]*(del.y[i]+one.D[i]) )*Del.y[,i]^2%*%t(y.xi[i,])
		norm( eq )
		},method="BFGS")$par
		}
	Del <- matrix(Del,p,p)
	Del.y <-  apply(y.xi,1,function(x) Del%*%x)
	if (equal){
	nu <- optim(nu1,function(x){
		-sum(log(2*dmvt(y,xi,S,x,log=F)*pt(del.y + t(one)%*%Del.y^3, x))) 
		},method="L-BFGS-B",lower=.01,upper=100)$par
	nu1 <- nu2 <- nu
	}
	else{
	nu <- optim(c(nu1,nu2),function(x){
		-sum(log(2*dmvt(y,xi,S,x[1],log=F)*pt(del.y + t(one)%*%Del.y^3, x[2]))) 
		},method="L-BFGS-B",lower=.01,upper=100)$par
	nu1 <- nu[1] ; nu2 <- nu[2]
	}
	LL.new <-  sum(log(2*dmvt(y,xi,S,nu1,log=F)*pt(del.y + t(one)%*%Del.y^3, nu2))) # log-likelihood function
	count <- count +1 
	dif <- abs(LL.new/LL-1)
	LL <- LL.new
	}
	la1 <- sqrt.mt(S,inverse=F)%*%del ; la2=(Del%*%(sqrt.mt(S,inverse=F)))^3%*%one
	if(equal){
	aic <- -2 * LL.new + 2 * (3*p+p*(p+1)/2+1)
	bic <- -2 * LL.new + log(n) * (3*p+p*(p+1)/2+1)
	} else {
	aic <- -2 * LL.new + 2 * (3*p+p*(p+1)/2+2)
	bic <- -2 * LL.new + log(n) * (3*p+p*(p+1)/2+2)
	}
	end <- proc.time()[3]
	time <- end-begin
	list(xi=xi, S=S, la1=la1 , la2=la2, nu1=nu1, nu2=nu2, loglik=LL.new, AIC=aic, BIC=bic, iter=count,elapsed=as.numeric(time))
	}
	if (family=='MFSSN')
	fit <- EM.SSN(y, xi=xi, S=S, la1=la1 , la2=la2, iter.max=iter.max, tol=tol, CML=CML )
	if (family=='MFSSTN')
	fit <- EM.SSTN(y, xi=xi, S=S, la1=la1 , la2=la2, nu= nu1, iter.max=iter.max, tol=tol, CML=CML )
	if (family=='MFSSSLN')
	fit <- EM.SSSLN(y, xi=xi, S=S, la1=la1 , la2=la2, nu= nu1, iter.max=iter.max, tol=tol, CML=CML )
	if (family=='MFSSCN')
	fit <- EM.SSCN(y, xi=xi, S=S, la1=la1 , la2=la2, nu1= nu1, nu2=nu2, equal=equal, iter.max=iter.max, tol=tol, CML=CML )
	if (family=='MFSST')
	fit <- EM.SST(y, xi=xi, S=S, la1=la1 , la2=la2, nu= nu1, iter.max=iter.max, tol=tol, CML=CML )
	if (family=='MFSSTT')
	fit <- EM.SSTT(y, xi=xi, S=S, la1=la1 , la2=la2, nu1= nu1, nu2=nu2, equal=equal, iter.max=iter.max, tol=tol, CML=CML )
	return(list(family=family,fit=fit))
	}




