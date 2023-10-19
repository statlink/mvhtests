\name{Log-likelihood ratio test for equality of two or more covariance matrices}
\alias{likel.cov}
\title{
Log-likelihood ratio test for equality of two or more covariance matrices
}

\description{
Log-likelihood ratio test for equality of two or more covariance matrices.
}

\usage{
likel.cov(x, ina, a = 0.05)
}

\arguments{
\item{x}{
A matrix containing Euclidean data.
}
\item{ina}{
A vector denoting the groups of the data.
}
\item{a}{
The significance level, set to 0.05 by default.
}
}

\details{
Tthe hypothesis test is that of the equality of at least two covariance matrices: \eqn{H_0:\pmb{\Sigma}_1=\ldots=\pmb{\Sigma}_k}. The algorithm is taken from Mardia, Bibby and Kent (1979, pg. 140). The log-likelihood ratio test is the multivariate generalization of Bartlett's test of homogeneity of variances. The test statistic takes the following form
\deqn{
-2log{\lambda}=n\log{\left|{\bf S}\right|}-\sum_{i=1}^kn_i\log{\left|{\bf S_i}\right|}=\sum_{i=1}^kn_i\log{\left|{\bf S}_i^{-1}{\bf S}\right|},
}
where \eqn{{\bf S}_i} is the \eqn{i}-th sample biased covariance matrix and \eqn{{\bf S}=n^{-1}\sum_{i=1}^kn_i{\bf S}_i} is the maximum likelihood estimate of the common covariance matrix (under the null hypothesis) with \eqn{n=\sum_{i=1}^kn_i}. The degrees of freedom of the asymptotic chi-square distribution are \eqn{\frac{1}{2}\left(p+1\right)\left(k-1\right)}.
}

\value{
A vector with the the test statistic, the p-value, the degrees of freedom and the critical value of the test.
}

\references{
Mardia K. V., Kent J. T., and Bibby J. M. (1979). Multivariate Analysis, pg. 140. London: Academic
Press.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{equal.cov}, \link{Mtest.cov}
}
}

\examples{
x <- as.matrix( iris[, 1:4] )
ina <- iris[, 5]
likel.cov(x, ina)
}

