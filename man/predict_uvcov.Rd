\name{predict_uvcov}
\alias{predict_uvcov}
\title{Obtain BLUPs for new levels of random effects with user specified variance covariance-matrices.}

\description{
  Obtain BLUPs for new levels of random effects with user specified variance covariance-matrices.
}

\usage{
  predict_uvcov(object,newrandom)
}

\arguments{
	\item{object}{is an object returned by lmer_uvcov.}
	\item{newrandom}{newrandom two level list with ids to be predicted and variance covariance matrix that contains information of 
                     these ids and the ids used to fit the model.}
}

\details{

Assume that the random effect \eqn{\boldsymbol u_j \sim N(\boldsymbol 0, \boldsymbol K_j)}
and the matrix \eqn{\boldsymbol K_j} is partitioned as follows:

\deqn{\boldsymbol u_j= \left[ {\begin{array}{c} \boldsymbol u_{j1} \\ \boldsymbol u_{j2} \end{array}}\right]}

and

\deqn{\boldsymbol K_j= \left[ {\begin{array}{cc} \boldsymbol K_{j11} & \boldsymbol K_{j21} \\ \boldsymbol K_{j22} & \boldsymbol K_{j12} \end{array}}\right]}

The BLUP for \eqn{\boldsymbol u_{j2}} can be obtained as:

\deqn{\boldsymbol E(\boldsymbol u_{j2}|\boldsymbol y_2)=\boldsymbol K_{j21} \boldsymbol K_{j11}^{-1} \hat {\boldsymbol u}_{j1} }

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

  #Predict 10/100 of records selected at random. 
  #The data were partitioned in 10 groups at random
  #and we predict individuals in group 2.

  fold<-2
  ytrn<-y[wheat.sets!=fold]
  ytst<-y[wheat.sets==fold]

  #######################################################################################
  #Marker based prediction
  #######################################################################################

  random<-list(mrk=list(K=G,id=names(ytrn)))
	
  out<-lmer_uvcov(ytrn,fixed="1",random=random)

  plot(ytrn,predict(out),xlab="Phenotype",ylab="Pred. Gen. Value")

  #Random effect list for prediction
  newrandom<-list(mrk=list(K=G,id=names(y)[wheat.sets==fold]))

  blup_tst<-predict_uvcov(out,newrandom)
  blup_tst<-blup_tst[[1]]
  points(ytst,blup_tst,col="red",pch=19)

  #Correlation in testing set
  cor(ytst,blup_tst)

  #######################################################################################
  #Pedigree based prediction
  #######################################################################################

  random<-list(ped=list(K=A,id=names(ytrn)))
	
  out<-lmer_uvcov(ytrn,fixed="1",random=random)

  plot(ytrn,predict(out),xlab="Phenotype",ylab="Pred. Gen. Value")

  #Random effect list for prediction
  newrandom<-list(ped=list(K=A,id=names(y)[wheat.sets==fold]))

  blup_tst<-predict_uvcov(out,newrandom)
  blup_tst<-blup_tst[[1]]
  points(ytst,blup_tst,col="red",pch=19)

  #Correlation in testing set
  cor(ytst,blup_tst)


  #######################################################################################
  #Markers + Pedigree based prediction
  #######################################################################################
 
  random<-list(mrk=list(K=G,id=names(ytrn)),
               ped=list(K=A,id=names(ytrn)))
	
  out<-lmer_uvcov(ytrn,fixed="1",random=random)

  plot(ytrn,predict(out),xlab="Phenotype",ylab="Pred. Gen. Value")

  #Random effect list for prediction
  newrandom<-list(mrk=list(K=G,id=names(y)[wheat.sets==fold]),
                  ped=list(K=A,id=names(y)[wheat.sets==fold]))

  blup_tst<-predict_uvcov(out,newrandom)
  blup_tst<-blup_tst[[1]]+blup_tst[[2]]
  points(ytst,blup_tst,col="red",pch=19)

  #Correlation in testing set
  cor(ytst,blup_tst)

}

}

