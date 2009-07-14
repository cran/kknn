\name{train.kknn}
\alias{train.kknn}
\alias{plot.train.kknn}
\alias{print.train.kknn}
\alias{predict.train.kknn}
\alias{summary.train.kknn}
\title{Training kknn}
\description{
  Training of kknn method via leave-one-out crossvalidation. 
}
\usage{
train.kknn(formula, data, kmax = 11, distance = 2, kernel = NULL,
	ykernel = NULL, contrasts = c('unordered' = "contr.dummy",
	ordered = "contr.ordinal"), ...)
}
\arguments{
\item{formula}{A formula object.}
\item{data}{Matrix or data frame. }
\item{kmax}{Maximum number of k. }
\item{distance}{Parameter of Minkowski distance.}
\item{kernel}{Kernel to use. Possible choices are "rectangular" 
	(which is standard unweighted knn), "triangular", "epanechnikov" 
	(or beta(2,2)), "biweight" (or beta(3,3)), "triweight" (or beta(4,4)), 
	"cos", "inv" and "gaussian".}
\item{ykernel}{Window width of an y-kernel, especially for prediction 
	of ordinal classes.}  	  
\item{contrasts}{A vector containing the 'unordered' and 'ordered' contrasts to use.}	
\item{\dots}{Further arguments passed to or from other methods.}
}
\value{
\code{train.kknn} returns a list-object of class \code{train.kknn} including
	the components 
\item{MISCLASS}{Matrix of misclassification errors.}
\item{MEAN.ABS}{Matrix of mean absolute errors.}
\item{MEAN.SQU}{Matrix of mean squared errors.}
\item{fitted.values}{List of predictions for all combinations of kernel and k.}
\item{best.parameters}{List containing the best parameter value for kernel and k.}
\item{response}{Type of response variable, one of \emph{continuous}, \emph{nominal} or \emph{ordinal}.}
\item{distance}{Parameter of Minkowski distance.}
\item{call}{The matched call.}
\item{terms}{The 'terms' object used.}
}


\references{Hechenbichler K. and Schliep K.P. (2004)  \emph{Weighted k-Nearest-Neighbor Techniques
 	and Ordinal Classification}, Discussion Paper 399, SFB 386, Ludwig-Maximilians University Munich
(\url{http://www.stat.uni-muenchen.de/sfb386/papers/dsp/paper399.ps})

Hechenbichler K. (2005)  \emph{Ensemble-Techniken und ordinale Klassifikation}, PhD-thesis
}
\author{ Klaus P. Schliep \email{K.P.Schliep@massey.ac.nz} }
\seealso{ \code{\link[kknn]{kknn}} and \code{\link[kknn]{simulation}} }
\examples{
library(kknn)

data(miete)
(train.con <- train.kknn(nmqm ~ wfl + bjkat + zh, data = miete, 
	kmax = 25, kernel = c("rectangular", "triangular", "epanechnikov",
	"gaussian", "rank")))
plot(train.con)
(train.ord <- train.kknn(wflkat ~ nm + bjkat + zh, miete, kmax = 25,
 	kernel = c("rectangular", "triangular", "epanechnikov", "gaussian", 
 	"rank")))
plot(train.ord)
(train.nom <- train.kknn(zh ~ wfl + bjkat + nmqm, miete, kmax = 25, 
	kernel = c("rectangular", "triangular", "epanechnikov", "gaussian", 
	"rank")))
plot(train.nom)

data(glass)
glass <- glass[,-1]
(fit.glass1 <- train.kknn(Type ~ ., glass, kmax = 15, kernel = 
	c("triangular", "rectangular", "epanechnikov"), distance = 1))
(fit.glass2 <- train.kknn(Type ~ ., glass, kmax = 15, kernel = 
	c("triangular", "rectangular", "epanechnikov"), distance = 2))
plot(fit.glass1)
plot(fit.glass2)

}
\keyword{ classif }% at least one, from doc/KEYWORDS
