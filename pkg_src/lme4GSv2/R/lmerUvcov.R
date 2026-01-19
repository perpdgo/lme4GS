#Definition of NOT IN 
'%notin%' <- Negate('%in%')

#Adapted from relfac.R in the lme4qtl package
#Removed relfac.svd
#Removed call to function llply from plyr package. llply was substituded by lapply
#Last Update Feb/04/2019

relfac <- function(mat, method.relfac = "auto", tol.relfac.evd = 1e-10) 
{

  # ?Matrix::chol
  # Returned value: a matrix of class 'Cholesky', i.e., upper triangular: R such that R'R = x.
  # Note that another notation is equivalent x = L L', where L is a lower triangular 
  # @ http://en.wikipedia.org/wiki/Cholesky_decomposition
  #
  # If the substitution is ZStar = Z L, then ZStar' = L' Z' = R Z'      
  
  # inc
  stopifnot(requireNamespace("Matrix", quietly = TRUE))
  
  rmat <- switch(method.relfac, 
    "chol" = relfac.chol(mat),
    "evd" = relfac.evd(mat, tol.relfac.evd),
    "auto" = {
      # filters
      mat.duplicated <- any(duplicated(lapply(1:ncol(mat), function(i) mat[, i])))
      
      # auto
      if(mat.duplicated) {
        rmat <- relfac.evd(mat, tol.relfac.evd)
        message("Computing relfac using evd")
      } else {
      
        rmat <- suppressWarnings(try({
          relfac.chol(mat)
        }, silent = TRUE))
        
        if(inherits(rmat, "try-error")) {
          rmat <- relfac.evd(mat, tol.relfac.evd)
          message("Computing relfac using evd")
        }else{
            message("Computing relfac using Cholesky")
        }
      }
      
      rmat
    },
    stop())
    
  rownames(rmat) <- rownames(mat)
  colnames(rmat) <- colnames(mat)
  
  return(rmat)
}

#Obtain relfactor using Cholesky factorization
relfac.chol <- function(mat)
{
  Matrix::chol(mat)
}


#Obtain relfactor using eigen-value decomposition
relfac.evd <- function(mat, tol)
{
  out <- eigen(mat, symmetric = TRUE)
  
  # clean eigen values
  ind <- (abs(out$values) < tol)
  if(any(ind)) {
    out$values[ind] <- 0
  }

  # clean eigen vectors only if eigen values were previously cleaned
  if(any(ind)) {
    ind_vec <- (abs(out$vectors) < tol)
    if(any(ind_vec)) {
      out$vectors[ind_vec] <- 0
    }
  }
  
  # return
  R <- diag(sqrt(out$values)) %*% t(out$vectors)
  return(as(R, "dgCMatrix"))
}

#Function to fit a mixed model:
#y=X*beta + Z_1 u_1 + Z_2 u_2 + ... + Z_q u_q + e 
#NA are not allowed in y

lmerUvcov <- function(formula, data = NULL, Uvcov = NULL, verbose=5L)
{
	#Match call
	mc<-match.call()
	
	#Checking inputs
	
	if(is.null(data)) stop("'data' can not be NULL\n")
	if(is.null(Uvcov)) stop("'Uvcov' can not be NULL\n")
	if(!is.list(Uvcov)) stop("'Uvcov' must be a list\n")
	if(length(Uvcov)<1) stop("'Uvcov' is an empty list\n")
	if(is.null(names(Uvcov)) | any(nchar(names(Uvcov))==0)) stop("Each element of the list 'Uvcov' must be named\n")
	
	#Step 1 parse the formula
	control<-lmerControl()
	control$checkControl$check.nobs.vs.rankZ <- "ignore"
	control$checkControl$check.nobs.vs.nlev <- "ignore"
	control$checkControl$check.nobs.vs.nRE <- "ignore"
  
  	# Parsed formula
  	parsedFormula <- lFormula(formula, data, control = control)
  	
  	#-------------------------------
  	# start of Uvcov-specific code
  	#-------------------------------
 
  	relnms <- names(Uvcov)
  	
  	#Sanity checks, assign factoring methods
	for(m in 1:length(relnms))
	{
		#Check that K is not null
		
		if(is.null(Uvcov[[m]]$K))
		{
			stop("K matrix associated with ", relnms[m], " can not be NULL\n")
		}
		
		if(!is.matrix(Uvcov[[m]]$K))
		{	
			stop("K object associated with ", relnms[m], " must be a matrix\n")
		}
		
		if(nrow(Uvcov[[m]]$K)!=ncol(Uvcov[[m]]$K))
		{
			stop("K matrix associated with ", relnms[m], " must be square\n")
		}
		
		#Check row and column names
		if(is.null(rownames(Uvcov[[m]]$K))) 
		{
			stop("rownames for matrix K associated with ", relnms[m], " can not be NULL\n")
		}
		
		if(is.null(colnames(Uvcov[[m]]$K))) 
		{
			stop("colnames for matrix K associated with ", relnms[m], " can not be NULL\n")
		}
		
		if(any(rownames(Uvcov[[m]]$K)!=colnames(Uvcov[[m]]$K))) 
		{
			stop("row and columns names for matrix K associated with ",relnms[m], " do not match\n")
		}

		#Check the factoring method and assign 'auto' if it was not assigned
		if(is.null(Uvcov[[m]]$f_method))
		{
			Uvcov[[m]]$f_method<-"auto"
			message("Factoring method for matrix K associated with ",relnms[m], " was set to 'auto'")
		}else{
			if(!(Uvcov[[m]]$f_method%in%c("auto","chol","evd"))) 
			{
				stop("Unsupported factoring method for matrix K associated with ",
				     relnms[m], ", supported factoring methods are: 'auto', 'chol' or 'evd'.\n")
			}
		}
	}
	
	#list to store relfactors
  	Tfac <- list()
  	
  	# list of factors
  	# The order match those in "Ztlist", see below 
  	flist <- parsedFormula$reTrms[["flist"]]   
  	
  	ind <- relnms %notin% names(flist)
  	
  	if(any(ind))
  	{
		stop(paste(relnms[ind],collapse=","), " not included in the formula\n")  		
  	}
  	
  	## random-effects design matrix components
    Ztlist <- parsedFormula$reTrms[["Ztlist"]]
    
    #Now we know that all the elements in relnms are included in the list of factors
    	
	for(m in 1:length(relnms))
	{
		relnms_m <- relnms[m]
				
		#Check the position of relnms_m in the list of factors
		j <- which(names(flist)==relnms_m)
		
		#This should be equal, if not stop
		if(names(flist)[j]!=relnms_m) stop("Ordering problem in Ztlist\n")
		
		z<-levels(parsedFormula$fr[,relnms_m])
		
		#Subset
		index<-rownames(Uvcov[[m]]$K)%in%z
		Tfac[[m]]<-Uvcov[[m]]$K[index,index]
		
		#Check unique
		if(any(duplicated(rownames(Tfac[[m]])))) stop("K matrix associated with ", relnms_m,
													   " has some duplicated ids\n")
													   
		if(length(z)!=nrow(Tfac[[m]])) stop(relnms_m, " has ", length(z), " levels ", 
											  " and the K matrix has ",nrow(Tfac[[m]]), " levels ")
		
		#Sorting
		index<-order(as.factor(rownames(Tfac[[m]])))
		Tfac[[m]]<-Tfac[[m]][index,index]
		
		if(any(rownames(Tfac[[m]])!=z)) stop("Ordering problem\n")
  
		Tfac[[m]] <- Matrix::Matrix(Tfac[[m]], sparse = TRUE)
		
		#R=L' in the case of Cholesky or
		#R=Lambda^0.5 * Gamma', in the case of eigen-value decomposition, K=Gamma*Lambda*Gamma'
		Tfac[[m]] <- relfac(Tfac[[m]],Uvcov[[m]]$f_method,1e-10)
		colnames(Ztlist[[j]])<-parsedFormula$fr[,relnms_m]
		
		#Compute ZStar', note ZStar = Z L, then ZStar' = L' Z' = R Z' in the case of Cholesky
		#note ZStar=Z Gamma Lamda^0.5 , then ZStar' = Lambda^0.5 Gamma' Z' = R Z' in the case of eigen-value decomposition  
		#In general ZStar' = R Z'
		
		Ztlist[[j]] <-  Tfac[[m]] %*% Ztlist[[j]]
		
	}
	
	names(Tfac)<-relnms
	
	parsedFormula$reTrms[["Ztlist"]] <- Ztlist
	
  	
  	parsedFormula$reTrms[["Zt"]] <- do.call(rbind, Ztlist)
  
  	#-------------------------------
  	# end of Uvcov-specific code
  	#-------------------------------
  	
  	#2 Deviance
  	
  	devianceFunction<-do.call(mkLmerDevfun, c(parsedFormula,
    					          list(start = NULL, verbose = TRUE, control = control)))
    					          
    #3 The optimization module
	optimizerOutput<-optimizeLmer(devianceFunction,verbose=verbose)

	#4 Output module 
	out<-mkMerMod(rho = environment(devianceFunction), opt = optimizerOutput,
                  reTrms = parsedFormula$reTrms, fr = parsedFormula$fr)
                  
    attr(out@frame,"relnms")<-relnms          #the name of the elements in the list 'Uvcov' provided by user
    attr(out@frame,"flist")<-flist            #list of factors
    
    cls <- "lmerUvcov"
    
    ans <- do.call(new, list(Class=cls, relfac = Tfac,
                             frame=out@frame, flist=out@flist, cnms=out@cnms, Gp=out@Gp,
                             theta=out@theta, beta=out@beta,u=out@u,lower=out@lower,
                             devcomp=out@devcomp, pp=out@pp,resp=out@resp,optinfo=out@optinfo))
    ans@call <- evalq(mc)
    ans
	
}


#Function to extract random effects for objects of class lmerUvcov
#for checking the arguments see the documentation of the ranef function

ranefUvcov<-function(object, postVar = FALSE, drop = FALSE, whichel = names(ans), ...)
{
	#cat("I am here\n")
	
	#Note that ranef provides the BLUPs for the model y=X*beta + ZStar*uStar+e
	#with ZStar = Z L in the case of Cholesky and ZStar= Z Lambda Gamma^0.5 in the case of eigen-value decomposition
	#ZStar' = L' Z' = R Z' in the case of Cholesky
	#ZStar' = Lambda^0.5 Gamma' Z' = R Z' in the case of eigen-value decomposition
	
	#Note u = L*uStar=R' * uStar, in the case of the Cholesky decomposition
	#u =  Lambda Gamma^0.5 uStar = R' * uStar in the case of eigen-value decomposition
	#In general we have u=R' * uStar
	
	ans <- ranef(as(object, "merMod"), postVar, drop = FALSE)
	ans <- ans[whichel]

	#R factors
	rf <- object@relfac
  
	for (nm in names(rf)) 
	{
		#blups in the transformed model (uStar)
		dm <- data.matrix(ans[[nm]])
		cn <- colnames(dm)
		rn <- rownames(dm)
		if(any(colnames(rf[[nm]])!=rownames(nm))) stop("Ordering problem in ")
		
		#blups in the original model (u)
		#Note that we need to take the transpose of R
		#u=R' * uStar
		
		dm <- as.matrix(t(rf[[nm]]) %*% dm)
		colnames(dm) <- cn
		rownames(dm) <- rn
		ans[[nm]] <- data.frame(dm, check.names = FALSE)
		
	}#End for

	if (drop)
	{
	  ans <- lapply(ans, function(el)
					     {
						   if (ncol(el) > 1) return(el)
						   pv <- drop(attr(el, "postVar"))
						   el <- drop(as.matrix(el))
						   if (!is.null(pv))
							  attr(el, "postVar") <- pv
						   el
					      } #End of function definition
					)#End of lapply
	  }#End if
	  
	  ans
}

#Overwrite generic ranef in S4

setMethod("ranef", signature(object = "lmerUvcov"),ranefUvcov)

#Obtain BLUPs for new levels of random effects
#object is an object returned by lmerUvcov
#newrandom two level list with ids to be predicted and variance covariance matrix that contains information
#of these ids and the ids used to fit the model. 

ranefUvcovNew<-function(object,Uvcov)
{

	#Check inputs
	if(!inherits(object, "lmerUvcov")) stop("object must be of class 'lmerUvcov'\n")
	
	if(is.null(Uvcov)) stop("'Uvcov' can not  be NULL\n")
	if(!is.list(Uvcov)) stop("'Uvcov' must be a list\n")
	if(length(Uvcov)<1) stop("'Uvcov' is an empty list\n")
	if(is.null(names(Uvcov)) | any(nchar(names(Uvcov))==0)) stop("Each element of the list 'Uvcov' must be named\n")
	
	relnms <- names(Uvcov)
  	
  	#Sanity checks, assign factoring methods
	
	for(m in 1:length(relnms))
	{
		#Check that K is not null
		
		if(is.null(Uvcov[[m]]$K))
		{
			stop("K matrix associated with ", relnms[m], " can not be NULL\n")
		}
		
		if(!is.matrix(Uvcov[[m]]$K))
		{	
			stop("K object associated with ", relnms[m], " must be a matrix\n")
		}
		
		if(nrow(Uvcov[[m]]$K)!=ncol(Uvcov[[m]]$K))
		{
			stop("K matrix associated with ", relnms[m], " must be square\n")
		}
		
		#Check row and column names
		if(is.null(rownames(Uvcov[[m]]$K))) 
		{
			stop("rownames for matrix K associated with ", relnms[m], " can not be NULL\n")
		}
		
		if(is.null(colnames(Uvcov[[m]]$K))) 
		{
			stop("colnames for matrix K associated with ", relnms[m], " can not be NULL\n")
		}
		
		if(any(rownames(Uvcov[[m]]$K)!=colnames(Uvcov[[m]]$K))) 
		{
			stop("row and columns names for matrix K associated with ",relnms[m], " do not match\n")
		}
	}

	tol.evd <- 1e-10
	blupsTrn <- ranefUvcov(object, postVar = FALSE, drop = FALSE)
	
	ind <- relnms %notin% names(blupsTrn)
  	
  	if(any(ind))
  	{
		stop(paste(relnms[ind],collapse=","), " not included in the formula used to fit the model \n")  		
  	}
  	
  	#At this point we are sure that all the elements in Uvcov all have associated blups 
	
	out<-list()
	
	for(m in 1:length(relnms))
	{
		relnms_m<-relnms[m]
		
		j<-which(names(blupsTrn)==relnms_m)
		blupTrn<-blupsTrn[[j]]
		
		if(ncol(blupTrn)>1) stop(relnms_m," only intercept models are supported\n")
				
		K<-Uvcov[[m]]$K
			
		if(any(!rownames(blupTrn)%in%rownames(K))) 
		{
			stop("not all levels for ",relnms_m, " in training set are present in provided K matrix\n")
		}
			
		trn<-rownames(blupTrn)
		all_ids<-rownames(K)
		tst<-setdiff(rownames(K),trn)
			
		if(length(tst)>0)
		{
			trn<-rownames(K)%in%trn
			K11<-K[trn,trn]
			K21<-K[!trn,trn]
				
			#Sort elements
			index<-order(rownames(blupTrn))
			blupTrn<-blupTrn[index,]
	
			index<-order(rownames(K11))
			K11<-K11[index,index]
		
			#Eigen_value decomposition
			outEigen<-eigen(K11)
	
			index<-outEigen$values>tol.evd
			outEigen$vectors<-outEigen$vectors[,index]
			outEigen$values<-outEigen$values[index]
			K11inv<-outEigen$vectors%*%diag(1.0/outEigen$values)%*%t(outEigen$vectors)
			rownames(K11inv)<-rownames(K11)
			colnames(K11inv)<-colnames(K11)
	
			#Sort the columns of K21, so that the column names match with that in K11inv
			index<-order(colnames(K21))
			K21<-K21[,index]
	
			if(any(colnames(K21)!=rownames(K11inv))) stop("Ordering problem!\n")
			if(any(colnames(K11inv)!=rownames(blupTrn))) stop("Ordering problem!\n")
	
			blupTst<-K21%*%K11inv%*%blupTrn
	
			blupTst<-blupTst[match(tst,rownames(blupTst)),,drop=FALSE]
			out[[m]]<-blupTst
			
		}else{
				stop("Testing set for ",relnms_m," has 0 elements\n")
		}						
	}
	
	names(out)<-relnms
	
	#Return the goodies
	out
}

