#Function that selects the bandwidth parameter of a Gaussian kernel
#K(i,j)=exp(-\theta (d_ij)^2))

theta_optim<-function(y,X)
{
  #y: response vector (n x 1)
  #X: matrix of markers (n x p)
  
  #check inputs
  if (is.null(rownames(X))) stop("X must have rownames\n")
  if (!is.character(rownames(X))) (rownames(X)<-as.character(rownames(X)))
  
  #Step 1:
  n <- nrow(X) #number of individuals
  p <- ncol(X) #number of markers
  
  #standardize the matrix of markers
  Z <- scale(X,center=TRUE,scale=TRUE) 
  #euclidian distance matrix
  D <- as.matrix(dist(Z, method = "euclidian"))/sqrt(p)
  #vector of values of theta
  theta <- setdiff(seq(0, max(D), length.out = 11),0)
  
  #Step 2:
  sol <- list()
  out <- list()
  nt  <- length(theta)
  
  #fit the model for each theta value
  for (i in 1:nt)
  {
    cat("i=",i,"\n")
    KG <- exp(-theta[i]*(D^2))	
    random <- list(kernel=list(K=KG, id=rownames(X)))
    sol[[i]] <- lmer_uvcov(y=y,fixed="1",random=random) 
  }
  
  #Step 3:
  LL <- rep(0, nt)
  
  for (i in 1:nt)
  {
    out[[i]] <- summary(sol[[i]])  
    LL[i] <- out[[i]]$logLik
  }
  
  LL <- as.numeric(LL)
  max.LL <- which.max(LL)
  theta.max <- theta[max.LL]
  K <- exp(-theta[max.LL]*(D^2))
  
  return(list(LL=LL, max.LL=LL[max.LL], theta=theta, 
              theta.max=theta.max, K=K))
}
