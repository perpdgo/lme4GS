\name{ranef_uvcov}
\alias{ranef_uvcov}
\title{Extract the conditional means of the random effects}

\description{
	A function to extract the conditional means of the random effects from a fitted model object.
}

\usage{
  ranef_uvcov(object)
}

\arguments{
 \item{object}{is an object returned by lmer_uvcov.}  
}

\details{

The function ranef extract the conditional means for the liner mixed effects model:

\ifelse{html}{\out{y=X&beta; + Z<sub>1</sub><sup>&ast;</sup> u<sub>1</sub><sup>&ast;</sup>+...+Z<sub>q</sub><sup>&ast;</sup> u<sub>q</sub><sup>&ast;</sup>+e,}}{\deqn{\boldsymbol y= \boldsymbol X \boldsymbol \beta+ \boldsymbol Z_1^* \boldsymbol u_1^* + \dots + \boldsymbol Z_q^* \boldsymbol u_q^* + \boldsymbol e,}{y= X beta + Z[1]^* u[1]^* + ... + Z[q]^* u[q]^* +e,}}

where  \ifelse{html}{\out{Z<sub>j</sub><sup>&ast;</sup>=Z<sub>j</sub> L<sub>j</sub>}}{\eqn{\boldsymbol Z_j^*= \boldsymbol Z_j \times \boldsymbol L_j}{Z[j]^* = Z[j] L[j]}}, with \ifelse{html}{\out{L<sub>j</sub>}}{\eqn{\boldsymbol L_j}{L[j]}} from Cholesky
factorization for \ifelse{html}{\out{K<sub>j</sub>}}{\eqn{\boldsymbol K_j}{K[j]}}. Alternatively,
\ifelse{html}{\out{Z<sub>j</sub><sup>&ast;</sup>=Z<sub>j</sub>&Gamma;<sub>j</sub>&Lambda;<sup>1/2</sup>}}{\eqn{\boldsymbol Z_j^*= \boldsymbol Z_j \times \boldsymbol \Gamma_j \boldsymbol \Lambda_j^{1/2}}{Z[j]^*=Z[j] Gamma[j] Lambda[j]^(1/2)}}, with
\ifelse{html}{\out{&Gamma;<sub>j</sub>}}{\eqn{\boldsymbol \Gamma_j}{Gamma[j]}} and 
\ifelse{html}{\out{&Lambda;<sub>j</sub>}}{\eqn{\boldsymbol \Lambda_j}{Lambda[j]}} 
the matrix of eigen-vectors and eigen-values obtained from the eigen-value decomposition 
for \ifelse{html}{\out{K<sub>j</sub>}}{\eqn{\boldsymbol K_j}{K[j]}}. So, 
the conditional means of the random effects in the linear mixed effects model:

\ifelse{html}{\out{y=X&beta; + Z<sub>1</sub> u<sub>1</sub> + ... + Z<sub>q</sub> u<sub>q</sub> + e,}}{\deqn{\boldsymbol y= \boldsymbol X \boldsymbol \beta+ \boldsymbol Z_1 \boldsymbol u_1 + \dots + \boldsymbol Z_q \boldsymbol u_q + \boldsymbol e,}{y=X beta + Z[1] u[1]  +...+ Z[q] u[q] + e,}}

are obtained as follows: \ifelse{html}{\out{&ucirc;<sub>j</sub>=L<sub>j</sub>&ucirc;<sub>j</sub><sup>&ast;</sup>}}{\eqn{\hat{\boldsymbol{u}}_j=\boldsymbol L_j \hat{\boldsymbol u}_j^*}{hat(u)[j]=L[j] hat(u)[j]^*}} if the Cholesky 
factorization is used or
\ifelse{html}{\out{&ucirc;<sub>j</sub>=&Gamma;<sub>j</sub>&Lambda;<sub>j</sub><sup>1/2</sup>&ucirc;<sub>j</sub><sup>&ast;</sup>}}{\eqn{\hat{\boldsymbol{u}}_j=\boldsymbol \Gamma_j \boldsymbol \Lambda_j^{1/2} \hat{\boldsymbol u}_j^*}{hat(u)[j]= Gamma[j] Lambda[j]^(1/2) hat(u)[j]^*}} if 
the the eigen-value decomposition is used.

}

\value{
A list of data frames, one for each grouping factor.
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

  out2<-ranef_uvcov(out)

  str(out2)

}

}

