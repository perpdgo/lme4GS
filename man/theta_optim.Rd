\name{theta_optim}
\alias{theta_optim}
\title{Selection of bandwidth parameter of Gaussian kernel.}

\description{
  Obtain the optimal value of the bandwidth parameter for the Gaussian kernel by cross-validation.
}


\usage{
  theta_optim<-function(y,X)
}

\arguments{
   \item{y}{(numeric, \eqn{n}{n}) the data-vector (NAs not allowed).}
   \item{X}{A matrix, \eqn{n}{p} that specifies the markers.}
}

\value{A list that contains
   \item{theta}{A vector that specifies the different values of \theta.}
   \item{theta.max}{The optimal value of \theta that maximizes the log-likelihood.}
   \item{LL}{A vector that specifies the values of the likelihood with its corresponding values of \theta.}
   \item{max.LL}{The maximum value of the log-likelihood.}
   \item{K}{Tha Gassuian Kernel compute with the optimal value of \theta.}  
}


\author{
Paulino Perez-Rodriguez
}

\examples{

\dontrun{
  library(BGLR)
  library(lme4GS)

#Example wheat

data(wheat)
X <- wheat.X 
A <- wheat.A 
y <- wheat.Y[,1]

rownames(X) <- rownames(A)

out_res <- theta_optim(y,X)
}

}
