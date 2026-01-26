#Function to fit a mixed model:
#y=X*beta + Z_1 u_1 + Z_2 u_2 + ... + Z_q u_q + e 
#NA are not allowed in y

lmerUvcov <- function(formula, data = NULL, Uvcov = NULL, verbose=0L)
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
			if(verbose>0L)
			{
				message("Factoring method for matrix K associated with ",relnms[m], " was set to 'auto'")
			}
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
		
		#Sorting must match levels obtained in z
     	index <- match(z, rownames(Tfac[[m]]))
		Tfac[[m]]<-Tfac[[m]][index,index]
		
		if(any(rownames(Tfac[[m]])!=z)) stop("Ordering problem\n")
  
		Tfac[[m]] <- Matrix::Matrix(Tfac[[m]], sparse = TRUE)
		
		#R=L' in the case of Cholesky or
		#R=Lambda^0.5 * Gamma', in the case of eigen-value decomposition, K=Gamma*Lambda*Gamma'
		Tfac[[m]] <- relfac(Tfac[[m]],Uvcov[[m]]$f_method,1e-10,verbose=verbose)
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
    					          list(start = NULL, verbose = verbose, control = control)))
    					          
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