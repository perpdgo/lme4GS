#Function that selects the bandwidth parameter of a Gaussian kernel
#K(i,j)=exp(-\theta (d_ij)^2))

#Option 2: Create a function to get the optimal value of theta of the
#kernel gaussiano in the general model:
#y= X*beta + Z_1 u_1 + ... + Z_q u_q

#If you want to calculate the Gaussian or exponential kernel, 
#you must put the m_id in the formula for the id markers.

theta_optim <- function(formula, data = NULL, Uvcov = NULL, 
                        kernel = list(D = NULL, kernel_type = "gaussian", 
                                      theta_seq = NULL, MRK = NULL))
{
  #Step1: check the m_id list 
  if(is.null(kernel)) stop("'kernel' can not be NULL\n")
  if(!is.list(kernel)) stop("'kernel' must be a list\n")
  if(length(kernel)<1) stop("'kernel' is an empty list\n")
  
  #check the components of m_id
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
    Z <- scale(X, center = TRUE, scale = TRUE)
    D <- as.matrix(dist(Z, method = "euclidian"))/sqrt(p)
  }
  
  
  #check the grid for theta and compute if the grid is NULL
  if(is.null(kernel$theta_seq)){
    theta <- setdiff(seq(0, max(D), length.out=11),0)
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
      cat("i=",i,"\n")
      KG <- exp(-theta[i]*(D^2))
      uvcovk <- list(m_id = list(K = KG))
      uvcovkk <- c(Uvcov, uvcovk)
      sol[[i]] <- lmerUvcov(formula = formulanew, data = datanew, 
                            Uvcov = uvcovkk)
    }  
  } else if (kernel$kernel_type=="exponential"){
    for (i in 1:nt){
      cat("i=",i,"\n")
      KE <- exp(-theta[i]*D)
      uvcovk <- list(m_id = list(K = KE))
      uvcovkk <- c(Uvcov, uvcovk)
      sol[[i]] <- lmerUvcov(formula = formulanew, data = datanew, 
                            Uvcov = uvcovkk)
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
  sol.max <- sol[[max.LL]] #the results of the optimal value of
  K.opt <- exp(-theta.max*(D^2)) 
  
  #K <- exp(-theta.max*(D^2))
  
  #can be added to return max.LL, theta.max
  return(list(LL=LL, LL.max=LL[max.LL],theta=theta,
              theta.max=theta.max,out=sol.max, K=K.opt))
}

