#Function that selects the bandwidth parameter of a Gaussian kernel
#K(i,j)=exp(-\theta (d_ij)^2))

theta_optim<-function(y, X, D = NULL, b = NULL)
{
  #y: response vector (n x 1)
  #X: matrix of markers (n x p)
  #D: distance matrix (p x p)
  #b: grid
  
  #check inputs
  if (is.null(rownames(X))) stop("X must have rownames\n")
  if (!is.character(rownames(X))) (rownames(X)<-as.character(rownames(X)))
  
  #Step 1:
  n <- nrow(X) #number of individuals
  p <- ncol(X) #number of markers
  
 #check the distance matrix
  if (is.null(D)) {
    #standardize the matrix of markers
    Z <- scale(X,center=TRUE,scale=TRUE) 
    #euclidian distance matrix
    D <- as.matrix(dist(Z, method = "euclidian"))/sqrt(p)
  }
  
  #grid of valures of theta
  if (is.null(b)){
    #vector of values of theta
    b <- setdiff(seq(0, max(D), length.out = 11),0)
  }
  
  #Step 2:
  sol <- list()
  out <- list()
  nt  <- length(b)
  
  #fit the model for each theta value
  for (i in 1:nt)
  {
    cat("i=",i,"\n")
    KG <- exp(-b[i]*(D^2))
    rownames(KG)<-colnames(KG)<-rownames(X)
    data<-data.frame(y=y,id=rownames(X))	
    sol[[i]] <- lmerUvcov(y~1+(1|id),data=data,Uvcov=list(id=KG))
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
  theta.max <- b[max.LL]
  K <- exp(-b[max.LL]*(D^2))
  
  return(list(LL=LL, max.LL=LL[max.LL], theta=b, 
              theta.max=theta.max, K=K))
}
