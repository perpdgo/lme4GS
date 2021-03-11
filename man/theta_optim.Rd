\name{theta_optim}
\alias{theta_optim}
\title{Selection of bandwidth parameter of Gaussian kernel.}

\description{
  Obtain the optimal value of the bandwidth parameter for the Gaussian and Exponential kernel by cross-validation.
}


\usage{
 theta_optim <- function(formula, data = NULL, Uvcov = NULL, 
                     kernel = list(D = NULL, kernel_type = "gaussian", 
                     theta_seq = NULL, MRK = NULL))
}

\arguments{
   \item{formula}{(a two-sided linear formula object describing both the fixed-effects and random-effects part of the model, with the
                  response on the left of a `~' operator and the terms, separated by `+' operators, on the right.  Random-effects
                  terms are distinguished by vertical bars (`|') separating expressions for design matrices from grouping factors.
                  The random-effect term associated with the markers with gaussian kernel and exponential kernel like variance-covariance 
                  matrix must be named as `k_id' in the formula and the data.}
   \item{data}{an optional data frame containing the variables named in `formula'.}
   \item{Uvcov}{a list to specify the variance-covariance matrix for random effects.}
   \item{kernel}{a list to specify the distance matrix of markers (D) this can be calculated or provided by the user, the kernel type can be gaussian or exponential, 
                 an optional grid values of the parameter bandwidth and the markers matrix.}
}

\value{A list that contains
   \item{theta}{A vector that specifies the different values of theta.}
   \item{theta.max}{The optimal value of theta that maximizes the log-likelihood.}
   \item{LL}{A vector that specifies the values of the likelihood with its corresponding values of theta.}
   \item{max.LL}{The maximum value of the log-likelihood.}
   \item{out}{The solution with the optimal value of theta.}
   \item{K}{The Kernel obtained with the optimal value of theta.}  
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

data=data.frame(y=y,k_id=rownames(X))

res<-theta.optim(y~(1|k_id),data=data,
                 kernel = list(D = NULL, kernel_type = "gaussian", 
                               theta.seq = NULL, MRK = X))

}

}
