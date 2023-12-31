\name{Box's M test for equality of two or more covariance matrices}
\alias{Mtest.cov}
\title{
Box's M test for equality of two or more covariance matrices
}

\description{
Box's M test for equality of two or more covariance matrices.
}

\usage{
Mtest.cov(x, ina, a = 0.05)
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
According to Mardia, Kent and Bibby (1979, pg. 140), it may be argued that if \eqn{n_i} is small, then the log-likelihood ratio test (function \code{\link{likel.cov}}) gives too much weight to the contribution of \eqn{{\bf S}}. This consideration led Box (1949) to propose another test statistic in place of that seen in \code{\link{likel.cov}} . Box's \eqn{M} is given by
\deqn{
M=\gamma\sum_{i=1}^k\left(n_i-1\right)\log{\left|{\bf S}_{i}^{-1}{\bf S}_p \right|},
}
where \eqn{\gamma=1-\frac{2p^2+3p-1}{6\left(p+1\right)\left(k-1\right)}\left(\sum_{i=1}^k\frac{1}{n_i-1}-\frac{1}{n-k}\right)} and \eqn{{\bf S}_{i}} and \eqn{{\bf S}_{p}} are the \eqn{i}-th unbiased covariance estimator and the pooled covariance matrix, respectively with \eqn{{\bf S}_p=\frac{\sum_{i=1}^k\left(n_i-1\right){\bf S}_i}{n-k}}. Box's \eqn{M} also has an asymptotic \eqn{\chi^2} distribution with \eqn{\frac{1}{2}\left(p+1\right)\left(k-1\right)} degrees of freedom. Box's approximation seems to be good if each \eqn{n_i} exceeds 20 and if \eqn{k} and \eqn{p} do not exceed 5 (Bibby and Kent (1979) pg. 140).
}

\value{
A vector with the the test statistic, the p-value, the degrees of freedom and the critical value of the test.
}

\references{
Mardia K.V., Kent J.T. and Bibby J.M. (1979). Multivariate Analysis. London: Academic Press.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{equal.cov}, \link{likel.cov}
}
}

\examples{
x <- as.matrix( iris[, 1:4] )
ina <- iris[, 5]
Mtest.cov(x, ina)
}


