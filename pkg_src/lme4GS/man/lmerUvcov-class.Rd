\name{lmerUvcov-class}

\docType{class}

\alias{lmerUvcov-class}
\alias{ranef,lmerUvcov-method}

\title{User defined variance covariance  mixed-effects model fits}

\description{A mixed-effects model fit by \code{\link{lmerUvcov}}.
  This class extends class \code{"\link[lme4:merMod-class]{merMod}"} class and includes one
  additional slot, \code{relfac}, which is a list of (left) Cholesky
  factors of the relationship matrices derived from 
  user provided variance covariance matrixes.}

\section{Objects from the Class}{Objects are created by calls to the
  \code{\link{lmerUvcov}} function.}

\section{Extends}{
Class \code{"\link[lme4:merMod-class]{merMod}"}, directly.
}

\section{Methods}{
  \describe{
   \item{ranef}{\code{signature(object = "lmerUvcov")}: incorporates
      the user defined variance covariance matrix into the random effects as returned for the object
      viewed as a \code{"\link[lme4:merMod-class]{merMod}"} object.}
  }
}

\seealso{\code{\link{lmerUvcov}}
}

\examples{
showClass("lmerUvcov")
}

\keyword{classes}

