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
#g(y|u)=X*beta + Z_1 u_1 + Z_2 u_2 + ... + Z_q u_q + e 
#u = [u1|u2|u3|...|uq]
#NA are not allowed in y

glmerUvcov <- function(formula, data = NULL, Uvcov = NULL, family = gaussian, verbose = 1L, 
                       start = NULL, nAGQ = 1L, contrasts = NULL, devFunOnly = FALSE){
  #subset, weights, na.action, offset,mustart, etastart,)
  
  # Match call
  mc <- match.call()
  
  # Checking inputs
  if(is.null(data)) stop("'data' can not be NULL\n")
  if(is.null(Uvcov)) stop("'Uvcov' can not be NULL\n")
  if(!is.list(Uvcov)) stop("'Uvcov' must be a list\n")
  if(length(Uvcov)<1) stop("'Uvcov' is an empty list\n")
  if(is.null(names(Uvcov)) | any(nchar(names(Uvcov))==0)) stop("Each element of the list 'Uvcov' must be named\n")
  
  # family-checking code duplicated here and in glFormula (for now)
  # full stop is gaussian default is not changed, ask later how to 
  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame(2))
  if( is.function(family)) family <- family()
  if (isTRUE(all.equal(family, gaussian()))) {
    ## redirect to lmerUvcov (with warning)
    #warning("calling glmer() with family=gaussian (identity link) as a shortcut to lmer() is deprecated;",
    #        " please call lmer() directly")
    #mc[[1]] <- quote(lme4GS::lmerUvcov)
    #mc["family"] <- NULL            # to avoid an infinite loop
    #return(eval(mc, parent.frame()))
    stop("Use lmeUvcov when calling family = gaussian (identity link)")
  }
  
  
  # Step 1: parse the formula
  control <- glmerControl()
  control$checkControl$check.nobs.vs.rankZ <- "ignore"
  control$checkControl$check.nobs.vs.nlev <- "ignore"
  control$checkControl$check.nobs.vs.nRE <- "ignore"
  
  # Parsed Formula 
  parsedFormula <- glFormula(formula = formula, data = data, family = family,
                             contrasts = contrasts, control = control)
  
  # Check the response variable for correct input
  
  if (is.matrix(y <- model.response(parsedFormula$fr))
      && ((family$family != "binomial" && ncol(y) > 1) ||
          (ncol(y) >2))) {
    stop("can't handle matrix-valued responses: consider using refit()")
  }
  
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
  
  # list to store relfactors
  Tfac <- list()
  
  # list of factors
  # The order match those in "Ztlist", see below
  # flist in parsedFormula is the list of grouping factors used in the random effects terms
  flist <- parsedFormula$reTrms[["flist"]]
  
  ind <- relnms %notin% names(flist)
  
  # Checks to see if there are any random effects that are included in Uvcov but not in the formula
  # this ignores any random effects with no defined Uvcov structure included in formula
  if(any(ind))
  {
    stop(paste(relnms[ind],collapse=","), " not included in the formula\n")  		
  }
  
  ## random-effects design matrix components
  # List of components of the transpose of the random effects model matrix separated by random effects
  Ztlist <- parsedFormula$reTrms[["Ztlist"]]
  
  # Now we know that all the elements in relnms are included in the list of factors
  # we can then compute Zstar considering that the function glFormula returns Z'
  for(m in 1:length(relnms))
  {
    relnms_m <- relnms[m]
    
    #Check the position of relnms_m in the list of factors
    j <- which(names(flist)==relnms_m)
    
    #This should be equal, if not stop
    if(names(flist)[j]!=relnms_m) stop("Ordering problem in Ztlist\n")
    
    z <- levels(parsedFormula$fr[,relnms_m])
    
    # Subset
    index <- rownames(Uvcov[[m]]$K)%in%z
    Tfac[[m]]<-Uvcov[[m]]$K[index,index]
    
    # Check unique
    if(any(duplicated(rownames(Tfac[[m]])))) stop("K matrix associated with ", relnms_m,
                                                  " has some duplicated ids\n")
    
    if(length(z)!=nrow(Tfac[[m]])) stop(relnms_m, " has ", length(z), " levels ", 
                                        " and the K matrix has ",nrow(Tfac[[m]]), " levels ")
    
    #Sorting must match levels obtained in z
    index <- match(z, rownames(Tfac[[m]]))
    Tfac[[m]] <- Tfac[[m]][index,index]
    
    if(any(rownames(Tfac[[m]])!=z)) stop("Ordering problem\n")
    
    Tfac[[m]] <- Matrix::Matrix(Tfac[[m]], sparse = TRUE)
    
    # R=L' in the case of Cholesky or
    # R=Lambda^0.5 * Gamma', in the case of eigen-value decomposition, K=Gamma*Lambda*Gamma'
    Tfac[[m]] <- relfac(Tfac[[m]],Uvcov[[m]]$f_method,1e-10)
    colnames(Ztlist[[j]]) <- parsedFormula$fr[,relnms_m]
    
    # Compute ZStar', note ZStar = Z L, then ZStar' = L' Z' = R Z' in the case of Cholesky
    # note ZStar=Z Gamma Lamda^0.5 , then ZStar' = Lambda^0.5 Gamma' Z' = R Z' in the case of eigen-value decomposition  
    # In general ZStar' = R Z'
    
    Ztlist[[j]] <-  Tfac[[m]] %*% Ztlist[[j]]
    
  }
  
  names(Tfac) <- relnms
  
  # We change the values of the model matrices obtained originally and substitute the Zstars
  parsedFormula$reTrms[["Ztlist"]] <- Ztlist 
  
  parsedFormula$reTrms[["Zt"]] <- do.call(rbind, Ztlist)
  
  #-------------------------------
  # end of Uvcov-specific code
  #-------------------------------
  
  ## create deviance function for covariance parameters (theta)
  nAGQinit <- if(control$nAGQ0initStep){0L} else{1L}
  devfun <- do.call(mkGlmerDevfun, c(parsedFormula, list(verbose = verbose,
                                                         control = control,
                                                         nAGQ = nAGQinit)))
  
  if (nAGQ==0 && devFunOnly) return(devfun)
  
  ## optimize deviance function over covariance parameters
  ## FIXME: perhaps should be in glFormula instead??
  if (is.list(start)) {
    start.bad <- setdiff(names(start),c("theta","fixef"))
    if (length(start.bad)>0) {
      stop(sprintf("bad name(s) for start vector (%s); should be %s and/or %s",
                   paste(start.bad,collapse=", "),
                   shQuote("theta"),
                   shQuote("fixef")),call.=FALSE)
    }
    if (!is.null(start$fixef) && nAGQ==0)
      stop("should not specify both start$fixef and nAGQ==0")
  }
  
  ## FIX ME: allow calc.derivs, use.last.params etc. if nAGQ=0
  if(control$nAGQ0initStep) {
    opt <- optimizeGlmer(devfun,
                         optimizer = control$optimizer[[1]],
                         ## DON'T try fancy edge tricks unless nAGQ=0 explicitly set
                         restart_edge=if (nAGQ==0) control$restart_edge else FALSE,
                         boundary.tol=if (nAGQ==0) control$boundary.tol else 0,
                         control = control$optCtrl,
                         start=start,
                         nAGQ = 0,
                         verbose=verbose,
                         calc.derivs=FALSE)
  }
  
  # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## -------------------- Function updatStart  ---------------------
  # used later in optimizing the deviance function
  updateStart <- function(start, theta) {
    if (is.numeric(start)) {
      theta
    } else {
      if (!is.null(start$theta))
        start$theta <- theta
      start
    }
  }
  
  # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  if(nAGQ > 0L) {
    ## update deviance function to include fixed effects as inputs
    devfun <- updateGlmerDevfun(devfun, parsedFormula$reTrms, nAGQ = nAGQ)
    str(devfun)
    if (control$nAGQ0initStep) {
      start <- updateStart(start,theta=opt$par)
    }
    ## if nAGQ0 was skipped
    ## we don't actually need to do anything here, it seems --
    ## getStart gets called again in optimizeGlmer
    
    if (devFunOnly) return(devfun)
    ## reoptimize deviance function over covariance parameters and fixed effects
    opt <- optimizeGlmer(devfun,
                         optimizer = control$optimizer[[2]],
                         restart_edge=control$restart_edge,
                         boundary.tol=control$boundary.tol,
                         control = control$optCtrl,
                         start=start,
                         nAGQ=nAGQ,
                         verbose = verbose,
                         stage=2,
                         calc.derivs=control$calc.derivs,
                         use.last.params=control$use.last.params)
  } 
  
  #tmp<<-control$calc.derivs
  #cc <- if (!control$calc.derivs){NULL} else {
  #  if (verbose > 10) cat("checking convergence\n")
  #  checkConv(attr(opt,"derivs"),opt$par,
  #            ctrl = control$checkConv,
  #            lbound=environment(devfun)$lower)
  #}
  cc <- cat("checking convergence\n")
      checkConv(attr(opt,"derivs"),opt$par,
                ctrl = control$checkConv,
                lbound=environment(devfun)$lower)
  
  ## prepare output
  out <- mkMerMod(rho = environment(devfun), opt = opt, reTrms = parsedFormula$reTrms, 
                  fr = parsedFormula$fr, lme4conv = cc)
  
  attr(out@frame,"relnms") <- relnms          #the name of the elements in the list 'Uvcov' provided by user
  attr(out@frame,"flist") <- flist            #list of factors
  
  cls <- "glmerUvcov"
  
  ans <- do.call(new, list(Class = cls, relfac = Tfac,
                           frame = out@frame, flist = out@flist, cnms = out@cnms, Gp = out@Gp,
                           theta = out@theta, beta = out@beta,u = out@u,lower = out@lower,
                           devcomp = out@devcomp, pp = out@pp,resp = out@resp,optinfo = out@optinfo))
  ans@call <- evalq(mc)
  ans
  
}

#Function to extract random effects for objects of class glmerUvcov
#for checking the arguments see the documentation of the ranef function

ranef_glmer_Uvcov<-function(object, postVar = FALSE, drop = FALSE, whichel = names(ans), ...)
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

setMethod("ranef", signature(object = "glmerUvcov"),ranef_glmer_Uvcov)




