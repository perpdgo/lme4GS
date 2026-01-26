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

setMethod("ranef",
          signature(object = "merUvcov"),
          ranefUvcov)
          
#Obtain BLUPs for new levels of random effects
#object is an object returned by lmerUvcov
#newrandom two level list with ids to be predicted and variance covariance matrix that contains information
#of these ids and the ids used to fit the model. 

ranefUvcovNew<-function(object,Uvcov)
{

	#Check inputs
	if(!inherits(object, "merUvcov")) stop("object must be of class 'merUvcov'\n")
	
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

