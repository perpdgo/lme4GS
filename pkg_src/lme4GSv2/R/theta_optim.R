
theta_optim <- function(formula, data = NULL, Uvcov = NULL, 
                        kernel = list(D = NULL, kernel_type = "gaussian", 
                                      theta_seq = NULL, MRK = NULL), 
                        verbose_lmer=0L, verbose_grid_search = 0L)
{
  #Step1: check the k_id list 
  if(is.null(kernel)) stop("'kernel' can not be NULL\n")
  if(!is.list(kernel)) stop("'kernel' must be a list\n")
  if(length(kernel)<1) stop("'kernel' is an empty list\n")
  
  #check the components of k_id
  #check the distance matrix
  if (!is.null(kernel$D)){
    if (!is.matrix(kernel$D)){
      stop("'D' must be a matrix\n")
    }
  }
  
  #compute the distance matrix if D is NULL
  else{
    X <- kernel$MRK
    
    if(is.null(X)){
      stop("The marker matrix 'MRK' can not be NULL\n")
    }
    
    if(!is.matrix(X)){
      stop("The marker matrix 'MRK' must be a matrix\n")
    }
    
    n <- nrow(X)
    p <- ncol(X)
    D <- as.matrix(dist(X, method = "euclidian"))/sqrt(p)
  }
  
  
  #check the grid for theta and compute if the grid is NULL
  if(is.null(kernel$theta_seq)){
    rho <- seq(0.15,0.90,length.out=20)
    theta<-(-log(rho))
    theta<-sort(theta) 
  } else
    theta <- kernel$theta_seq
  
  #Step2: compute the gaussian kernel for each value of theta
  
  sol <- list()
  out <- list()
  nt <- length(theta)
  formulanew <- formula  
  datanew <- data
  
  #fit the model for each value of theta
  #check the type of kernel
  
  if (kernel$kernel_type=="gaussian"){
    for (i in 1:nt){
      if (verbose_grid_search>0L){
        message("Case ",i,"/",length(theta),"\n")
        message("theta=",theta[i],"\n")
      }
      KG <- exp(-theta[i]*D^2)
      uvcovk <- list(k_id = list(K = KG))
      uvcovkk <- c(Uvcov, uvcovk)
      sol[[i]] <- lmerUvcov(formula = formulanew, data = datanew, 
                            Uvcov = uvcovkk,verbose=verbose_lmer)
    }  
  } else if (kernel$kernel_type=="exponential"){
    for (i in 1:nt){
      if(verbose_grid_search>0L){
        message("Case ",i,"/",length(theta),"\n")
        message("theta=",theta[i],"\n")
      }
      KE <- exp(-theta[i]*D)
      uvcovk <- list(k_id = list(K = KE))
      uvcovkk <- c(Uvcov, uvcovk)
      sol[[i]] <- lmerUvcov(formula = formulanew, data = datanew, 
                            Uvcov = uvcovkk,verbose=verbose_lmer)
    }  
  } else {
    stop("The type of kernel must be 'gaussian' or 'exponential'\n")
  }
  
  #step3: get the value of theta that maximizes the log-likelihood
  LL <- rep(0,nt)
  
  for (i in 1:nt){
    out[[i]] <- summary(sol[[i]])
    LL[i] <- out[[i]]$logLik
  }
  
  LL <- as.numeric(LL)
  max.LL <- which.max(LL)
  theta.max <- theta[max.LL]
  sol.max <- sol[[max.LL]] #Fitted model for optimal value of the bandwidth 
  
  if(kernel$kernel_type=="gaussian")
  {
  		K.opt <- exp(-theta.max*(D^2)) 
  }else if (kernel$kernel_type=="exponential"){
  		K.opt<- exp(-theta[i]*D)
  }
   
  #Return the goodies
  return(list(LL=LL, LL.max=LL[max.LL],theta=theta,
              theta.max=theta.max,fm=sol.max,K=K.opt))
}
