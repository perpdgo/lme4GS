#Function to fit a mixed model:
#y=X*beta + Z_1 u_1 + Z_2 u_2 + ... + Z_k u_k + e 
#NA are not allowed in y

lmer_uvcov<-function(y,fixed="1",random,verbose=5L)
{
	
	#Check inputs
	if(class(random)!="list") stop("random should be a list\n")
	if(is.null(names(random)) | any(nchar(names(random))==0)) stop("Each element of the list 'random' must be named\n")

	#Check the factoring method and assign 'auto' if it was not assigned
	for(k in 1:length(random))
	{
		if(is.null(random[[k]]$f_method))
		{
			random[[k]]$f_method<-"auto"
			message("Factoring method for vcov matrix set to 'auto' for ",names(random)[k])
		}else{
			if(!(random[[k]]$f_method%in%c("auto","chol","evd"))) stop("Supported factoring methods are: 'auto', 'chol' or 'evd'.\n")
		}
	}
	
	#Step 1 parse the formula
	control<-lmerControl()
	control$checkControl$check.nobs.vs.rankZ <- "ignore"
	control$checkControl$check.nobs.vs.nlev <- "ignore"
	control$checkControl$check.nobs.vs.nRE <- "ignore"
	
	f<-"y~"
	f<-paste0(f,fixed)
	
	for(k in 1:length(random))
	{
		f<-paste0(f,"+ (1|random[[",k,"]]$id)")
	}
	
	parsedFormula<-lFormula(f,control=control)
	
	#Step 2, deviance function
	
	for(k in 1:length(random))
	{
	    #If K is NULL then we assign the Identity matrix
	    if(!is.null(random[[k]]$K))
	    {
			random[[k]]$id<-as.factor(random[[k]]$id)
			ids.unique<-levels(random[[k]]$id)
			random[[k]]$K<-Matrix::Matrix(random[[k]]$K[ids.unique,ids.unique],sparse=TRUE)
	    }
	}
	
	#SOME THING VERY STRANGE HERE, IT IS POSSIBLE THAT THE ORDER OF ELEMENTS IN Ztlist 
	#ARE NOT THE SAME THAT IN THE LIST random
	
	Ztlist<-parsedFormula$reTrms[["Ztlist"]]
	fnmns<-names(parsedFormula$reTrms[["flist"]])
	
	#The names as I was expecting
	expected_names<-paste0("random[[",1:length(random),"]]$id")
	
	for(k in 1:length(random))
	{
		#In which position of the list fnmns is random[[k]]$id
		j<-which(fnmns[k]==expected_names)
		
		#If K is NULL then we assign the Identity matrix
		if(!is.null(random[[j]]$K))
		{
			Ztlist[[k]]<-relfac(random[[j]]$K,random[[j]]$f_method,1e-10)%*%Ztlist[[k]]
			
			#CHECK HERE, confirm if the names of columns are assigned correctly when 
			#factors are present
			colnames(Ztlist[[k]])<-random[[j]]$id
			
			#Ztlist[[k]]<-relfac(random[[j]]$K)%*%Ztlist[[k]]
		}
	}
	
	parsedFormula$reTrms[["Ztlist"]]<-Ztlist
	parsedFormula$reTrms[["Zt"]]<-do.call(rbind, Ztlist)  
	
	devianceFunction<-do.call(mkLmerDevfun, c(parsedFormula,
    					          list(start = NULL, verbose = TRUE, control = control)))
	
	#3 The optimization module
	optimizerOutput<-optimizeLmer(devianceFunction,verbose=verbose)

	#4 Output module 
	out<-mkMerMod(rho = environment(devianceFunction), opt = optimizerOutput,
                  reTrms = parsedFormula$reTrms, fr = parsedFormula$fr)
                  
    
    #Add attributes
    attr(out@frame,"enmns")<-expected_names   #Expected names
    attr(out@frame,"rnmns")<-names(random)    #the name of the elements in the list 'random' provided by user
    attr(out@frame,"fnmns")<-fnmns            #the name of the elements in the list created by lme4
    
    return(out)               
}
