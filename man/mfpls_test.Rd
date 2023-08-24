\name{fplsdc_test}
\alias{fplsdc_test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Test of functional partial least squares through distance covariance.
}
\description{
Test of independence for functional partial least squares through distance covariance(FPLS-DC). The null distribution of statistic is approximated by gamma distribution and the scalar and shape parameter is approximated by permutation procedure.
}
\usage{
fplsdc_test(X,Y,index.num = NULL,per.num = 50)
}
\arguments{
  \item{X}{A numeric matrix, or dataframe. Each column is an observation vector.}
  \item{Y}{A numeric matrix, or dataframe.}
  \item{index.num}{The index vector of variable for testing. If \code{index.num} = NULL, the function return the p.value of all variables.}
  \item{per.num}{Number of permutation replications for the estimation of scalar and shape parameter of gamma distribution.}
}

\value{
\item{p.value}{The p.value vector for the test.}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
library(MASS)
library(FPLSDC)
  n <- 200
  feature_num <- 50
  sigma_sqrt <- matrix(NA,feature_num,feature_num)
  for (i in seq(feature_num)) {
    for (j in seq(feature_num)) {
      if (i == j){
        sigma_sqrt[i,j] = 1
      }
      else
        sigma_sqrt[i,j] = 0.5^abs(i-j)
    }
  }
  X <- mvrnorm(n = n, mu = rep(0,feature_num), Sigma = sigma_sqrt)
  Y <- matrix(rnorm(n*200),n,200)
### Independent case
fplsdc_test1 <- fplsdc_test(X,Y,index.num = 1)
### Dependent case
  for (j in 51:150) {
        Y[,j] <-  Y[,j] + ((2/sqrt(12))*((0.002*j)^(-1))*X[,1] + (2/sqrt(12))*((0.002*j)^(-1))*(X[,2])+(2/sqrt(12))*((0.002*j)^(-1))*X[,3])^2
    }
fplsdc_test2 <- fplsdc_test(X,Y,index.num = c(1,2,3))
}
