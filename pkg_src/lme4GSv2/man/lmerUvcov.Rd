\name{lmerUvcov}
\alias{lmerUvcov}
\title{Fits a linear mixed model with user specified variance covariance-matrices.}

\description{
  Fits a linear mixed model with user specified variance covariance-matrices.
}


\usage{
  lmerUvcov(formula, data = NULL, Uvcov = NULL,verbose=5L)
}

\arguments{
   \item{formula}{a two-sided linear formula object describing both the fixed-effects and random-effects part of the model, with the
                  response on the left of a `~' operator and the terms, separated by `+' operators, on the right.  Random-effects
                  terms are distinguished by vertical bars (`|') separating expressions for design matrices from grouping factors.}
   \item{data}{an optional data frame containing the variables named in `formula'.} 
   \item{Uvcov}{list.}
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

\ifelse{html}{\out{y=X&beta; + Z<sub>1</sub><sup>&ast;</sup> u<sub>1</sub><sup>&ast;</sup>+...+Z<sub>q</sub><sup>&ast;</sup> u<sub>q</sub><sup>&ast;</sup>+e,}}{\deqn{\boldsymbol y= \boldsymbol X \boldsymbol \beta+ \boldsymbol Z_1^* \boldsymbol u_1^* + \dots + \boldsymbol Z_q^* \boldsymbol u_q^* + \boldsymbol e,}{y= Xbeta + Z[1]^* u[1]^* + ... + Z[q]^* u[q]^* +e,}}

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

\references{

Caamal-Pat D., P. Perez-Rodriguez, J. Crossa, C. Velasco-Cruz, S. Perez-Elizalde, M. Vazquez-Pena. 2021.
lme4GS: An R-Package for Genomic Selection. \emph{Front. Genet.} \bold{12}:680569. doi: 10.3389/fgene.2021.680569
doi: 10.3389/fgene.2021.680569

}

\examples{

\donttest{

library(BGLR)
library(lme4GS)

########################################################################
#Example 1, wheat 
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

########################################################################
#Example 2, mice 
########################################################################

data(mice)

X<-mice.X
Z<-scale(X,center=TRUE,scale=TRUE)
G<-tcrossprod(Z)/ncol(Z)

mice.pheno$m_id<-rownames(G)
mice.pheno$Obesity.BMI<-scale(mice.pheno$Obesity.BMI)


fm2<-lmerUvcov(Obesity.BMI~GENDER+Litter+(1|cage)+(1|m_id),
               data=mice.pheno,Uvcov=list(m_id=list(K=G)))
summary(fm2)

plot(mice.pheno$Obesity.BMI,predict(fm2))

ranef(fm2)

}

}

