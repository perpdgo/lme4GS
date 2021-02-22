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

\ifelse{html}{\out{y=X&beta; + Z<sub>1</sub> u<sub>1</sub> + ... + Z<sub>q</sub> u<sub>q</sub> + e,}}{\deqn{\boldsymbol y= \boldsymbol X \boldsymbol \beta+ \boldsymbol Z_1 \boldsymbol u_1 + \dots + \boldsymbol Z_q \boldsymbol u_q + \boldsymbol e,}{y=X beta + Z[1] u[1] +...+ Z[q] u[q] + e,}}

where \eqn{\boldsymbol y}{y} is the response vector, \eqn{\boldsymbol X}{X} is the matrix for
fixed effects,  \ifelse{html}{\out{&beta;}}{\eqn{\boldsymbol\beta}{beta}} is the vector of fixed effects, 
\ifelse{html}{\out{Z<sub>j</sub>}}{\eqn{\boldsymbol Z_j}{Z[j]}} is a design matrix for random effects, 
\ifelse{html}{\out{u<sub>j</sub>}}{\eqn{\boldsymbol u_j}{u[j]}} is a vector of random effects, 
\eqn{j=1,\dots,q}{j=1,...,q}. We assume that 
\ifelse{html}{\out{u<sub>j</sub>&#126;N(0,&sigma;<sup>2</sup><sub>j</sub> K <sub>j</sub>)}}{\eqn{\boldsymbol u_j \sim N(\boldsymbol 0, \sigma^2_j \boldsymbol K_j)}{u_j ~ N(0, sigma[j]^2] K[j])}},
\eqn{j=1,\dots,q}{j=1,...,q} and 
\ifelse{html}{\out{e&#126;N(0,&sigma;<sup>2</sup><sub>e</sub>I)}}{\eqn{\boldsymbol e \sim N(\boldsymbol 0, \sigma^2_e \boldsymbol I)}{e ~ N(0, sigma[e]^2 I)}}.

The linear mixed model can be re-written as:

\ifelse{html}{\out{y=y=X&beta; + Z<sub>1</sub><sup>&ast;</sup> u<sub>1</sub><sup>&ast;</sup>+...+Z<sub>q</sub><sup>&ast;</sup> u<sub>q</sub><sup>&ast;</sup>+e,}}{\deqn{\boldsymbol y= \boldsymbol X \boldsymbol \beta+ \boldsymbol Z_1^* \boldsymbol u_1^* + \dots + \boldsymbol Z_q^* \boldsymbol u_q^* + \boldsymbol e,}{y= Xbeta + Z[1]^* u[1]^* + ... + Z[q]^* u[q]^* +e,}}

where \ifelse{html}{\out{Z<sub>j</sub><sup>&ast;</sup>=Z<sub>j</sub> L<sub>j</sub>}}{\eqn{\boldsymbol Z_j^*= \boldsymbol Z_j \times \boldsymbol L_j}{Z[j]^* = Z[j] L[j]}}, with \ifelse{html}{\out{L<sub>j</sub>}}{\eqn{\boldsymbol L_j}{L[j]}} from Cholesky 
factorization for \ifelse{html}{\out{K<sub>j</sub>}}{\eqn{\boldsymbol K_j}{K[j]}}. Alternatively,
\ifelse{html}{\out{Z<sub>j</sub><sup>&ast;</sup>=Z<sub>j</sub>&Gamma;<sub>j</sub>&Lambda;<sup>1/2</sup>}}{\eqn{\boldsymbol Z_j^*= \boldsymbol Z_j \times \boldsymbol \Gamma_j \boldsymbol \Lambda_j^{1/2}}{Z[j]^*=Z[j] Gamma[j] Lambda[j]^(1/2)}}, with 
\ifelse{html}{\out{&Gamma;<sub>j</sub>}}{\eqn{\boldsymbol \Gamma_j}{Gamma[j]}} and 
\ifelse{html}{\out{&Lambda;<sub>j</sub>}}{\eqn{\boldsymbol \Lambda_j}{Lambda[j]}} the matrix of 
eigen-vectors and eigen-values obtained from the eigen-value decomposition 
for \ifelse{html}{\out{K<sub>j</sub>}}{\eqn{\boldsymbol K_j}{K[j]}}. The factorization method for
\ifelse{html}{\out{K<sub>j</sub>}}{\eqn{\boldsymbol K_j}{K[j]}} is selected automatically at runtime. 

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

