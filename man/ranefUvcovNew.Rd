\name{ranefUvcovNew}
\alias{ranefUvcovNew}
\title{Obtain BLUPs for new levels of random effects with user specified variance covariance-matrices.}

\description{
  Obtain BLUPs for new levels of random effects with user specified variance covariance-matrices.
}

\usage{
  ranefUvcovNew(object,Uvcov)
}

\arguments{
	\item{object}{is an object returned by lmerUvcov.}
	\item{Uvcov}{two level list with ids to be predicted and variance covariance matrix that contains information of 
                     these ids and the ids used to fit the model.}
}

\details{

Assume that the random effect 
\ifelse{html}{\out{u<sub>j</sub> &#126; N(0, &sigma;<sup>2</sup><sub>j</sub>K<sub>j</sub>)}}{\eqn{\boldsymbol u_j \sim N(\boldsymbol 0,\sigma_j^2 \boldsymbol K_j)}{u[j] ~N(0,sigma[j]^2 K[j])}} and the matrix \ifelse{html}{\out{K<sub>j</sub>}}{\eqn{\boldsymbol K_j}{K[j]}} is partitioned as follows:

\ifelse{html}{\out{u<sub>j</sub>=(u<sub>j1</sub> u<sub>j2</sub>)'}}{\deqn{\boldsymbol u_j= \left( {\begin{array}{c} \boldsymbol u_{j1} \\ \boldsymbol u_{j2} \end{array}}\right)}{u[j]=(u[j1] u[j2])'}}

and

\ifelse{html}{\figure{Kj.png}}{\deqn{\boldsymbol K_j= \left[ {\begin{array}{cc} \boldsymbol K_{j11} & \boldsymbol K_{j12} \\ \boldsymbol K_{j21} & \boldsymbol K_{j22} \end{array}}\right]}{K[j]=(K[j11] K[j12]; K[j21] K[j22])}}

The BLUP for \eqn{\boldsymbol u_{j2}}{u[j2]} can be obtained as:

\ifelse{html}{\out{E(u<sub>j2</sub>|y<sub>1</sub>)=K<sub>j21</sub>K<sub>j11</sub><sup>-1</sup>&ucirc;<sub>j1</sub>}}{\deqn{\boldsymbol E(\boldsymbol u_{j2}|\boldsymbol y_1)=\boldsymbol K_{j21} \boldsymbol K_{j11}^{-1} \hat {\boldsymbol u}_{j1} }{E(u[j2] | y[1])=K[j21] K[j11]^(-1) hat(u)[j1]}}

}


\author{
Paulino Perez-Rodriguez
}

\examples{

\dontrun{
 
library(BGLR)
library(lme4GS)

data(wheat)
X<-wheat.X
Z<-scale(X,center=TRUE,scale=TRUE)
G<-tcrossprod(Z)/ncol(Z)
A<-wheat.A
rownames(G)<-colnames(G)<-rownames(A)
y<-wheat.Y[,1]

#Predict 10/100 of records selected at random. 
#The data were partitioned in 10 groups at random
#and we predict individuals in group 2.

fold<-2
y_trn<-y[wheat.sets!=fold]
y_tst<-y[wheat.sets==fold]

A_trn=A[wheat.sets!=fold,wheat.sets!=fold]
G_trn=G[wheat.sets!=fold,wheat.sets!=fold]

pheno_trn=data.frame(y_trn=y_trn,m_id=rownames(A_trn),a_id=rownames(G_trn))

#######################################################################################
#Marker based prediction
#######################################################################################
	
fm1<-lmerUvcov(y_trn~1+(1|m_id),data=pheno_trn,Uvcov=list(m_id=list(K=G_trn)))

plot(pheno_trn$y_trn,predict(fm1),xlab="Phenotype",ylab="Pred. Gen. Value")

#BLUP for individuals in the testing set
blup_tst<-ranefUvcovNew(fm1,Uvcov=list(m_id=list(K=G)))
blup_tst<-blup_tst$m_id[,1]

#Comparison
#Check the names
names(y_tst)<-rownames(G)[wheat.sets==fold]
blup_tst<-blup_tst[match(names(y_tst),names(blup_tst))]

yHat_tst<-fixef(fm1)[1]+blup_tst
points(y_tst,yHat_tst,col="red",pch=19)

#Correlation in testing set
cor(y_tst,yHat_tst)

#######################################################################################
#Pedigree based prediction
#######################################################################################

fm2<-lmerUvcov(y_trn~1+(1|a_id),data=pheno_trn,Uvcov=list(a_id=list(K=A_trn)))

plot(pheno_trn$y_trn,predict(fm2),xlab="Phenotype",ylab="Pred. Gen. Value")

#BLUP for individuals in the testing set
blup_tst<-ranefUvcovNew(fm2,Uvcov=list(a_id=list(K=A)))
blup_tst<-blup_tst$a_id[,1]

#Comparison
#Check the names
names(y_tst)<-rownames(A)[wheat.sets==fold]
blup_tst<-blup_tst[match(names(y_tst),names(blup_tst))]

yHat_tst<-fixef(fm2)[1]+blup_tst
points(y_tst,yHat_tst,col="red",pch=19)

#Correlation in testing set
cor(y_tst,yHat_tst)

#######################################################################################
#Markers + Pedigree based prediction

fm3<-lmerUvcov(y_trn~1+(1|m_id)+(1|a_id),data=pheno_trn,
               Uvcov=list(m_id=list(K=G_trn),a_id=list(K=A_trn)))

plot(pheno_trn$y_trn,predict(fm3),xlab="Phenotype",ylab="Pred. Gen. Value")

#BLUP for individuals in the testing set
blup_tst<-ranefUvcovNew(fm3,Uvcov=list(m_id=list(K=G),a_id=list(K=A)))

blup_tst_m<-blup_tst$m_id[,1]
blup_tst_a<-blup_tst$a_id[,1]

#Comparison
#Check the names
names(y_tst)<-rownames(A)[wheat.sets==fold]
blup_tst_m<-blup_tst_m[match(names(y_tst),names(blup_tst_m))]
blup_tst_a<-blup_tst_a[match(names(y_tst),names(blup_tst_a))]

yHat_tst<-fixef(fm3)[1] + blup_tst_m + blup_tst_a
points(y_tst,yHat_tst,col="red",pch=19)

#Correlation in testing set
cor(y_tst,yHat_tst)

}

}

