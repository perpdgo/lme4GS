#Obtain BLUPs for new levels of random effects
#object is an object returned by lmer_uvcov
#newrandom two level list with ids to be predicted and variance covariance matrix that contains information
#of these ids and the ids used to fit the model. 

predict_uvcov<-function(object,newrandom)
{
	#Check inputs
	if(class(newrandom)!="list") stop("newrandom should be a list\n")
	if(is.null(names(newrandom)) | any(nchar(names(newrandom))==0)) stop("Each element of the list 'newrandom' must be named\n")
	
	#Check that the names of newrandom are present in object
	
	if(any(!names(newrandom)%in%attr(object@frame,"rnmns"))) 
	{
	    stop("All named elements of the list 'newrandom' must be in the list 'random' provided to lmer_custom_varcov\n")
	}
	
	tol.evd <- 1e-10
	
	blup<-list()
	
	Ztlist<-getME(object,"Ztlist")
	
	for(k in 1:length(newrandom))
	{
		j<-which(attr(object@frame,"rnmns")==names(newrandom)[k])
		label<-attr(object@frame,"enmns")[j]
		
		#Note that ranef provides the BLUPs for the model y=X*beta + ZStar*uStar+e
	
		out_ranef<-ranef(object,whichel=label)
		out_ranef<-out_ranef[[1]]
		uStar_trn<-out_ranef[,1]
		names(uStar_trn)<-rownames(out_ranef)
		
		#The BLUPs are obtained when multiplying ZStar= Z*relfac(K) by uStar_trn
		#u=Z*relfac(K)*uStar_trn. Note that ZStar is given in Ztlist
		
		#In which position in Ztlist the Zt that we need
		j<-which(attr(object@frame,"fnmns")==label)
		
		#Extract Z
		Z<-t(Ztlist[[j]])
		
		if(any(colnames(Z)!=names(uStar_trn))) stop("Ordering problem\n")
		u_trn<-as.vector(Z%*%uStar_trn)
		names(u_trn)<-rownames(Z)
		
		#Remove duplicated
		index<-duplicated(names(u_trn))
		u_trn<-u_trn[!index]
	
		K<-newrandom[[k]]$K
	
		trn<-names(u_trn)
		tst<-newrandom[[k]]$id
		all_ids<-unique(c(trn,tst))
	
		#Subset K
		subset<-rownames(K)%in%all_ids
		K<-K[subset,subset]
	
		trn<-rownames(K)%in%names(u_trn)
		K11<-K[trn,trn]
		K21<-K[!trn,trn]
	
		#Sort elements
		index<-order(names(u_trn))
		u_trn<-u_trn[index]
	
		index<-order(rownames(K11))
		K11<-K11[index,index]
		
		#Eigen_value decomposition
		out_eigen<-eigen(K11)
	
		index<-out_eigen$values>tol.evd
		out_eigen$vectors<-out_eigen$vectors[,index]
		out_eigen$values<-out_eigen$values[index]
		K11inv<-out_eigen$vectors%*%diag(1.0/out_eigen$values)%*%t(out_eigen$vectors)
		rownames(K11inv)<-rownames(K11)
		colnames(K11inv)<-colnames(K11)
	
		#Sort the columns of K21, so that the column names match with that in K11inv
		index<-order(colnames(K21))
		K21<-K21[,index]
	
		if(any(colnames(K21)!=rownames(K11inv))) stop("Ordering problem!\n")
		if(any(colnames(K11inv)!=names(u_trn))) stop("Ordering problem!\n")
	
		u_tst<-as.vector(K21%*%K11inv%*%u_trn)
		names(u_tst)<-rownames(K21)
	
	
		#FIXME: check that the manes of the output object are in the same order than id
		u_tst<-u_tst[match(tst,names(u_tst))]
		blup[[k]]<-u_tst
		
	}
	names(blup)<-names(newrandom)
	blup
}

ranef_uvcov<-function(object)
{
	 blup<-list()

         Ztlist<-getME(object,"Ztlist")

         labels<-attr(object@frame,"enmns")

         for(k in 1:length(labels))
         {
                #Note that ranef provides the BLUPs for the model y=X*beta + ZStar*uStar+e

                out_ranef<-ranef(object,whichel=labels[k])
                out_ranef<-out_ranef[[1]]
                uStar_trn<-out_ranef[,1]
                names(uStar_trn)<-rownames(out_ranef)

                #The BLUPs are obtained when multiplying ZStar= Z*relfac(K) by uStar_trn
                #u=Z*relfac(K)*uStar_trn. Note that ZStar is given in Ztlist

                #In which position in Ztlist the Zt that we need
                j<-which(attr(object@frame,"fnmns")==labels[k])

                #Extract Z
                Z<-t(Ztlist[[j]])

                if(any(colnames(Z)!=names(uStar_trn))) stop("Ordering problem\n")
                u_trn<-as.vector(Z%*%uStar_trn)
                names(u_trn)<-rownames(Z)

                #Remove duplicated
                index<-duplicated(names(u_trn))
                u_trn<-u_trn[!index]
		blup[[j]]<-u_trn
	}
	#Add the names
	names(blup)<-attr(object@frame,"rnmns")

	#Return the goodies
	return(blup)
}
