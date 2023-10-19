\name{Hypothesis test for two high-dimensional mean vectors}
\alias{sarabai}
\title{
Hypothesis test for two high-dimensional mean vectors
}

\description{
Hypothesis test for two high-dimensional mean vectors.
}

\usage{
sarabai(x1, x2)
}

\arguments{
\item{x1}{
A matrix containing the Euclidean data of the first group.
}
\item{x2}{
A matrix containing the Euclidean data of the second group.
}
}

\details{
High dimensional data are the multivariate data which have many variables (\eqn{p}) and usually a small number of observations (\eqn{n}). It also happens that \eqn{p>n} and this is the case here in this Section. We will see a simple test for the case of \eqn{p>n}. In this case, the covariance matrix is not invertible and in addition it can have a lot of zero eigenvalues.

The test we will see was proposed by Bai and Saranadasa (1996). Ever since, there have been some more suggestions but I chose this one for its simplicity. There are two datasets, \eqn{{\bf X}_1} and \eqn{{\bf X}_2} of sample sizes \eqn{n_1} and \eqn{n_2}, respectively. Their corresponding sample mean vectors and covariance matrices are \eqn{\bar{{\bf X}}_1}, \eqn{\bar{{\bf X}}_2} and \eqn{{\bf S}_1}, \eqn{{\bf S}_2} respectively. The assumption here is the same as that of the Hotelling's test we saw before.

Let us define the pooled covariance matrix at first, calculated under the assumption of equal covariance matrices
\eqn{
{\bf S}_n=\frac{\left(n_1-1\right){\bf S}_1+\left(n_2-1\right){\bf S}_2}{n}},
where \eqn{n=n_1+n_2}. Then define \eqn{B_n=\sqrt{ \frac{n^2}{\left(n+2\right)\left(n-1\right)}\left\lbrace\text{tr}\left({\bf S}_n^2\right)-
\frac{1}{n}\left[\text{tr}\left({\bf S}_n\right)\right]^2 \right\rbrace }}.
The test statistic is
\deqn{
Z=\frac{\frac{n_1n_2}{n_1+n_2}\left(\bar{{\bf X}}_1-\bar{{\bf X}}_2\right)^T\left(\bar{{\bf X}}_1-\bar{{\bf X}}_2\right)
-\text{tr}\left({\bf S}_n\right)}{\sqrt{\frac{2\left(n+1\right)}{n}}B_n}.
}

Under the null hypothesis (equality of the two mean vectors) the test statistic follows the standard normal distribution. Bai and Saranadasa (1996) established the asymptotic normality of the test statistics and showed that it has attractive power property when \eqn{p/n \rightarrow c < \infty} and under some restriction on the maximum eigenvalue of the common population covariance matrix. However, the requirement of \eqn{p} and \eqn{n} being of the same order is too restrictive to be used in the "large \eqn{p} small \eqn{n}" situation.
}

\value{
A vector with the test statistic and the p-value.
}

\references{
Bai Z. D. and Saranadasa H. (1996). Effect of high dimension: by an example of a two
sample problem. Statistica Sinica, 6(2): 311--329.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{hotel2T2}, \link{maov}, \link{el.test2}, \link{eel.test2}
}
}

\examples{
x1 <- matrix( rnorm(40 * 100), ncol = 100 )
x2 <- matrix( rnorm(50 * 100), ncol = 100 )
sarabai(x1, x2)
}

