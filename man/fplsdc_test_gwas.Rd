\name{fplsdc_test_gwas}
\alias{fplsdc_test_gwas}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
GWAS test of functional partial least squares through distance covariance.
}
\description{
GWAS test of independence for functional partial least squares through distance covariance(FPLS-DC). The null distribution of statistic is approximated by gamma distribution and the scalar and shape parameters is approximated by modified parameters, which are obtained through permutation procedure on a refered data (\code{X1,Y}). The probability of \code{X1} = 0,1,2 is 1/3.
}
\usage{
fplsdc_test_gwas(X,Y,index.num = NULL,per.num = 50)
}
\arguments{
  \item{X}{A numeric matrix, or dataframe. Each column is an observation vector.}
  \item{Y}{A numeric matrix, or dataframe.}
  \item{index.num}{The index vector of variable for testing. If \code{index.num} = NULL, the function return the p.value of all variables.}
  \item{per.num}{Number of permutation replications for the estimation of scalar and shape parameter of gamma distribution.}
}

\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{p.value}{The p.value vector for the GWAS test.}
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
### indepedent case
X <- matrix(sample(0:2, size = 200*10, replace = TRUE, prob = c(1/4, 1/4, 1/2)),200,10)
Y <- matrix(rnorm(200*200),200,200)
p.value1 <- fplsdc_test_gwas(X,Y)
### dependent case
Y <- matrix(rnorm(200*200),200,200)
for (j in 51:150) {
Y[,j] <- ((1/sqrt(2))*((0.002*j)^(-1))*(X[,1] - X[,5]))^2 + rnorm(200)
}
p.value2 <- fplsdc_test_gwas(X,Y)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
