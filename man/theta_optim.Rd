\name{theta_optim}
\alias{theta_optim}
\title{Selection of bandwidth parameter for Gaussian and exponential kernels.}

\description{
  Obtain the optimal value of the bandwidth parameter for the Gaussian and exponential kernels.
}


\usage{
  theta_optim(formula, data = NULL, Uvcov = NULL,
                       kernel = list(D = NULL, kernel_type = "gaussian", 
                                           theta_seq = NULL, MRK = NULL))
}

\arguments{
  \item{formula}{a two-sided linear formula object describing both the fixed-effects and random-effects part of the model, with the
                  response on the left of a `~' operator and the terms, separated by `+' operators, on the right.  Random-effects
                  terms are distinguished by vertical bars (`|') separating expressions for design matrices from grouping factors.}
   \item{data}{an optional data frame containing the variables named in `formula'.} 
   \item{Uvcov}{list.}
   \item{kernel}{list with the following elements, i)D: Distance matrix (can be NULL), 
   	                    ii) kernel_type: character, can be either "gaussian" or "exponential",
   	                    ii)theta_seq: sequence of values for theta from which we select the optimum (can be NULL), 
   	                    iv) MRK: marker matrix from wich Euclidean distance is computed (can be NULL).}
}

\value{A list that contains: 
	\item{LL}{Log-likelihood.}
	\item{LL.max}{Maximum of likelihood.}
	\item{theta}{Sequence of values for the bandwidth.}
	\item{theta.max}{Value of bandwidth when log-likelihood attains the maximum.}
	\item{fm}{Fitted model with the optimum bandwidth parameter.}
	\item{K.opt}{The kernel for the optimum bandwith parameter.}  
}


\author{
Paulino Perez-Rodriguez, Diana Caamal-Pat
}

\examples{
	
\dontrun{
library(BGLR)
library(lme4GS)

data(wheat)

l = 1
y = wheat.Y[,l]
X = wheat.X
A = wheat.A

rownames(X) <- rownames(A)

#model y=1*mu+Z_1*u_1+e, u_1~NM(0, \sigma_1*KG), KG: Gaussian kernel
wheat = data.frame(y=y, k_id=rownames(X))

fm1 <- thetaoptimmod(y~(1|k_id), data = wheat, Uvcov = NULL,
                    kernel = list(D = NULL, kernel_type = "gaussian", 
                                  theta_seq = NULL, MRK = X))

fm2 <- thetaoptimmod(y~(1|k_id), data = wheat, Uvcov = NULL,
                    kernel = list(D = NULL, kernel_type = "exponential", 
                                  theta_seq = NULL, MRK = X))
                                  
########################################################################
#model y=1*mu+Z_1*u_1+Z_2*u_2+e, u_1~NM(0, \sigma_1*KG), KG: Gaussian kernel
#u_2~NM(0, \sigma_2*A)

wheat = data.frame(y=y, a_id=rownames(A), k_id=rownames(X))

fm3 <- thetaoptimmod(y~(1|a_id)+(1|k_id), data = wheat, Uvcov = list(a_id=list(K=A)), 
                    kernel = list(D = NULL, kernel_type = "gaussian", 
                                  theta_seq = NULL, MRK = X))

########################################################################
data(mice)

X<-mice.X
A<-mice.A

mice.pheno$a_id<-rownames(A)
mice.pheno$k_id<-rownames(X)

out<-thetaoptimmod(Obesity.BMI~GENDER+Litter+(1|cage)+(1|k_id)+(1|a_id),data=mice.pheno,
               Uvcov=list(a_id=list(K=A)),kernel = list(D = NULL, kernel_type = "gaussian", 
                                                        theta_seq = NULL, MRK = X))
}
}
