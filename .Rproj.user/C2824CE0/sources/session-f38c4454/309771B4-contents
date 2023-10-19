\name{Log-likelihood ratio test for equality of one covariance matrix}
\alias{equal.cov}
\title{
Log-likelihood ratio test for equality of one covariance matrix
}

\description{
Log-likelihood ratio test for equality of one covariance matrix.
}

\usage{
equal.cov(x, Sigma, a = 0.05)
}

\arguments{
\item{x}{
A matrix containing Euclidean data.
}
\item{Sigma}{
The hypothesis covariance matrix.
}
\item{a}{
The significance level, set to 0.05 by default.
}
}

\details{
The hypothesis test is that the the sample covariance is equal to some specified covariance matrix: \eqn{H_0:\pmb{\Sigma}=\pmb{\Sigma}_0}, with \eqn{\pmb{\mu}} unknown. The algorithm for this test is taken from Mardia, Bibby and Kent (1979, pg. 126-127).
The test is based upon the log-likelihood ratio test. The form of the test is
\deqn{
-2\log{\lambda}=n \text{tr}\left\lbrace \pmb{\Sigma}_0^{-1}{\bf S}\right\rbrace-n\log{\left|\pmb{\Sigma}_0^{-1}{\bf S} \right|}-np,
}
where \eqn{n} is the sample size, \eqn{\pmb{\Sigma}_0} is the specified covariance matrix under the null hypothesis, \eqn{{\bf S}} is the sample covariance matrix and \eqn{p} is the dimensionality of the data (or the number of variables). Let \eqn{\alpha} and \eqn{g} denote the arithmetic mean and the geometric mean respectively of the eigenvalues of \eqn{\pmb{\Sigma}_0^{-1}{\bf S}}, so that \eqn{tr\left\lbrace \pmb{\Sigma}_0^{-1}{\bf S}\right\rbrace=p\alpha} and
\eqn{\left|\pmb{\Sigma}_0^{-1}{\bf S} \right|=g^p}, then the test statistic becomes
\deqn{
-2\log{\lambda}=np\left(\alpha-log{(g)}-1 \right).
}
The degrees of freedom of the \eqn{\chi^2} distribution are \eqn{\frac{1}{2}p\left(p+1\right)}.
}

\value{
A vector with the the test statistic, the p-value, the degrees of freedom and the critical value of the test.
}

\references{
Mardia K. V., Kent J. T., and Bibby J. M. (1979). Multivariate Analysis, pg. 126--127. London: Academic
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
\code{\link{likel.cov}, \link{Mtest.cov}
}
}

\examples{
x <- as.matrix( iris[, 1:4] )
s <- cov(x) * 1.5
equal.cov(x, s)
}

