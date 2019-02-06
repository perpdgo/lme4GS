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

  summary(out)

  plot(y,predict(out))

}

}

