\name{ifpls_dc}
\alias{ifpls_dc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Iterative functional partial least squares through
distance covariance.
}
\description{
Iterative functional partial least squares through
distance covariance(IFPLS-DC). IFPLS-DC can be used to presents a variable selection methodology by adding L0 constraint on \code{X}.
}
\usage{
ifpls_dc(X, Y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{A numeric matrix, or dataframe. Each column is an observation vector.}
  \item{Y}{A numeric matrix, or dataframe. The number of columns of Y is greater than the number of basis functions, the defult number is 5. }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{index}{Index vector with selected variables.}
\item{rank}{The index of all predictors ranked by IFPLS-DC.}
\item{dcovU}{The maximum of unbiased distance covariance correspond to selected variables.}
\item{par}{The parameters correpsond to the selected variables.}
\item{bhat}{The coefficient function for Y if the number of columns of Y is greater than the number of basis functions.}
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
  for (i in 1:n){
    for (j in 51:150) {
        Y[i,j] <-  Y[i,j] + ((2/sqrt(13))*((0.002*j)^(-1))*X[i,1] + (2/sqrt(13))*((0.002*j)^(-1))*(X[i,2])+(2/sqrt(13))*((0.002*j)^(-1))*X[i,3] + (1/sqrt(13))*((0.002*j)^(-1))*X[i,30])^2
    }
  }
  ## high dimensional case
  ifpls_select1 <- ifpls_dc(X,Y)
  ## one dimensional case
  ifpls_select2 <- ifpls_dc(X,Y[,52])
}

