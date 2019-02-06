\name{lmer_uvcov}
\alias{lmer_uvcov}
\title{Fits a linear mixed model with user specified variance covariance-matrices.}

\description{
  Fits a linear mixed model with user specified variance covariance-matrices.
}


\usage{
  lmer_uvcov(y,fixed="1",random,verbose=5L)
}

\arguments{
   \item{y}{(numeric, \eqn{n}{n}) the data-vector (NAs not allowed).}
   \item{fixed}{A formula that specifies fixed effects.}
   \item{random}{(list) This is a two-level list used to specify the variance covariance matrixes for random effects.}
   \item{verbose}{verbose output from optimizeLmer function?.}	
}

\details{

The routine fits the linear mixed model:

\ifelse{html}{\figure{lmm.png}}{{\deqn{\boldsymbol y= \boldsymbol X \boldsymbol \beta+ \boldsymbol Z_1 \boldsymbol u_1 + \dots + \boldsymbol Z_q \boldsymbol u_q + \boldsymbol e}}}

where \eqn{\boldsymbol y} is the response vector, \eqn{\boldsymbol X} is the matrix for
fixed effects,  \eqn{\boldsymbol\beta} is the vector of fixed effects, 
\eqn{\boldsymbol Z_j} is a design matrix for random effects, \eqn{\boldsymbol u_j} is
a vector of random effects, \eqn{j=1,\dots,q}. We assume that 
\eqn{\boldsymbol u_j \sim N(\boldsymbol 0, \sigma^2_j \boldsymbol K_j)}, 
\eqn{j=1,\dots,q} and \eqn{\boldsymbol e \sim N(\boldsymbol 0, \sigma^2_e \boldsymbol I)}.

The liner mixed model can be re-written as:

\deqn{\boldsymbol y= \boldsymbol X \boldsymbol \beta+ \boldsymbol Z_1^* \boldsymbol u_1^* + \dots + \boldsymbol Z_q^* \boldsymbol u_q^* + \boldsymbol e}

where \eqn{\boldsymbol Z_j^*= \boldsymbol Z_j \times \boldsymbol L_j}, with \eqn{\boldsymbol L_j} from Cholesky 
factorization for \eqn{\boldsymbol K_j}. Alternatively, 
\eqn{\boldsymbol Z_j^*= \boldsymbol Z_j \times \boldsymbol \Gamma_j \boldsymbol \Lambda_j^{1/2}}, with 
\eqn{\boldsymbol \Gamma_j} and \eqn{\boldsymbol \Lambda_j} the matrix of eigen-vectors and 
eigen-values obtained from the eigen-value decomposition for \eqn{\boldsymbol K_j}.

}

\value{
  An object of class merMod (more specifically,
  an object of \emph{subclass} \code{lmerMod}), for which many methods
  are available (e.g. \code{methods(class="merMod")})
}


\author{
Paulino Perez-Rodriguez
}

\examples{

\dontrun{
  library(BGLR)
  library(lme4GS)

  #Example 1, wheat 
  data(wheat)
  X<-wheat.X
  Z<-scale(X,center=TRUE,scale=TRUE)
  G<-tcrossprod(Z)/ncol(Z)
  A<-wheat.A
  rownames(G)<-colnames(G)<-rownames(A)
  y<-wheat.Y[,1]

  #id a vector with ids
  random<-list(mrk=list(K=G,id=rownames(G)),
               ped=list(K=A,id=rownames(A)))

  out<-lmer_uvcov(y,fixed="1",random=random)

  summary(out)

  plot(y,predict(out))

  #Example 2, mice

  data(mice)

  y<-mice.pheno$Obesity.BMI
  y<-scale(y,center=TRUE,scale=TRUE)
  y<-as.vector(y)

  X<-mice.X
  Z<-scale(X,center=TRUE,scale=TRUE)
  G<-tcrossprod(Z)/ncol(Z)

  gender<-mice.pheno$GENDER
  litter<-mice.pheno$Litter
  cage<-mice.pheno$cage

  fixed<-"1+gender+litter"

  #K=NULL means that we are using the identity matrix
  random<-list(cage=list(K=NULL,id=cage),
	           mrk=list(K=G,id=rownames(G)))

  out<-lmer_uvcov(y,fixed=fixed,random=random)

  plot(y,predict(out))
}

}

