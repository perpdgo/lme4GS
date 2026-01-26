\name{ranefUvcov}
\alias{ranefUvcov}
\title{Extract the conditional means of the random effects}

\description{
	A function to extract the conditional means of the random effects from a fitted model object.
}

\usage{
  ranefUvcov(object, postVar = FALSE, drop = FALSE, whichel = names(ans), ...)
}

\arguments{
 \item{object}{is an object returned by lmerUvcov.}  
 \item{postVar}{a logical argument indicating if the conditional variance-covariance matrices of the random effects should be added as an attribute.}
 \item{drop}{should components of the return value that would be data frames with a single column, usually a column called ``(Intercept)'', be returned as named vectors instead?.}
 \item{whichel}{character vector of names of grouping factors for which the random effects should be returned.}
 \item{...}{some methods for these generic functions require additional arguments.}
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

\references{

Caamal-Pat D., P. Perez-Rodriguez, J. Crossa, C. Velasco-Cruz, S. Perez-Elizalde, M. Vazquez-Pena. 2021.
lme4GS: An R-Package for Genomic Selection. \emph{Front. Genet.} \bold{12}:680569. doi: 10.3389/fgene.2021.680569
doi: 10.3389/fgene.2021.680569

}

\author{
Paulino Perez-Rodriguez
}

\examples{

\donttest{

library(BGLR)
library(lme4GS)

########################################################################
#Example wheat
########################################################################
data(wheat)
X<-wheat.X
Z<-scale(X,center=TRUE,scale=TRUE)
G<-tcrossprod(Z)/ncol(Z)
A<-wheat.A
rownames(G)<-colnames(G)<-rownames(A)
y<-wheat.Y[,1]

data<-data.frame(y=y,m_id=rownames(G),a_id=rownames(A))

fm1<-lmerUvcov(y~(1|m_id)+(1|a_id),data=data,
               Uvcov=list(m_id=list(K=G),a_id=list(K=A)))

summary(fm1)

#Predictions
plot(y,predict(fm1))

#Random effects
ranef(fm1)

#Equivalently
ranefUvcov(fm1)

}

}

