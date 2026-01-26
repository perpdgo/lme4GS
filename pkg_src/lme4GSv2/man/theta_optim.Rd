\name{theta_optim}
\alias{theta_optim}
\title{Selection of bandwidth parameter for Gaussian and exponential kernels.}

\description{
  Obtain the optimal value of the bandwidth parameter for the Gaussian and exponential kernels.
}


\usage{
  theta_optim(formula, data = NULL, Uvcov = NULL,
                       kernel = list(D = NULL, kernel_type = "gaussian", 
                                     theta_seq = NULL, MRK = NULL), 
                       verbose_lmer= 0L, verbose_grid_search=0L)
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
   \item{verbose_lmer}{integer scalar, verbose output from optimizeLmer function?. If '> 0' verbose output is generated during the optimization of the parameter estimates, 
                       default value is 0L.}
   \item{verbose_grid_search}{integer scalar, if '>0' verbose output is generated, default value is 0L.}
}

\value{A list that contains: 
	\item{LL}{Log-likelihood.}
	\item{LL.max}{Maximum of likelihood.}
	\item{theta}{Sequence of values for the bandwidth.}
	\item{theta.max}{Value of bandwidth when log-likelihood attains the maximum.}
	\item{fm}{Fitted model with the optimum bandwidth parameter.}
	\item{K.opt}{The kernel for the optimum bandwith parameter.}  
}

\references{

Caamal-Pat D., P. Perez-Rodriguez, J. Crossa, C. Velasco-Cruz, S. Perez-Elizalde, M. Vazquez-Pena. 2021.
lme4GS: An R-Package for Genomic Selection. \emph{Front. Genet.} \bold{12}:680569. doi: 10.3389/fgene.2021.680569
doi: 10.3389/fgene.2021.680569

}


\author{
Paulino Perez-Rodriguez, Diana Caamal-Pat
}

\examples{
	
\donttest{
library(BGLR)
library(lme4GS)

data(wheat)

y = wheat.Y[,1]
X = wheat.X
A = wheat.A

rownames(X) <- rownames(A)

#model y=1*mu+Z_1*u_1+e, u_1~NM(0, \sigma_1*KG), KG: Gaussian kernel
wheat = data.frame(y=y, k_id=rownames(X))

fm1 <- theta_optim(y~(1|k_id), data = wheat, Uvcov = NULL,
                   kernel = list(D = NULL, kernel_type = "gaussian", 
                                 theta_seq = seq(3,8,length.out=10), MRK = X),
                   verbose_lmer=0L,verbose_grid_search=1L)

plot(fm1$theta,fm1$LL,xlab=expression(theta),ylab="Log-Likelihood")
fm1$theta.max
fm1$LL.max


}
}
