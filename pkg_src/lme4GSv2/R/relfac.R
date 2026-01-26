#Definition of NOT IN 
'%notin%' <- Negate('%in%')

#Adapted from relfac.R in the lme4qtl package
#Removed relfac.svd
#Removed call to function llply from plyr package. llply was substituded by lapply
#Last Update Feb/04/2019

relfac <- function(mat, method.relfac = "auto", tol.relfac.evd = 1e-10,verbose=0L) 
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
        if(verbose>0L){
        	message("Computing relfac using evd")
	}
      } else {
      
        rmat <- suppressWarnings(try({
          relfac.chol(mat)
        }, silent = TRUE))
        
        if(inherits(rmat, "try-error")) {
          rmat <- relfac.evd(mat, tol.relfac.evd)
	  if(verbose>0L)
	  {
          	message("Computing relfac using evd")
	  }
        }else{
	    if(verbose>0L)
	    {
            	message("Computing relfac using Cholesky")
	    }
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
