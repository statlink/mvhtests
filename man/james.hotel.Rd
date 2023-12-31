\name{Relationship between the Hotelling's and James' tests}
\alias{james.hotel}
\title{
Relationship between the Hotelling's \eqn{T^2} and James' test
}

\description{
Relationship between the Hotelling's \eqn{T^2} and James' test.
}

\usage{
james.hotel(x1, x2)
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
Emerson (2009, pg. 76--81) mentioned a very nice result between the Hotelling's one sample \eqn{T^2} and James test for two mean vectors
\deqn{
J\left(\pmb{\mu}\right) = T_1^2\left(\pmb{\mu}\right) + T_2^2\left(\pmb{\mu}\right),
}
where \eqn{J\left(\pmb{\mu}\right)} is the James test statistic (James, 1954) and \eqn{T_1^2\left(\pmb{\mu}\right)} and \eqn{T_1^2\left(\pmb{\mu}\right)} are the two one sample Hotelling's \eqn{T^2} test statistic values (see function \code{\link{hotel1T2}}) for each sample from their common mean vector \eqn{\pmb{\mu}_c} (see the help file of \code{\link{james}}). In fact, James test statistic is found from minimizing the right hand side of the above expression with respect to \eqn{\pmb{\mu}}. The sum is mimized when \eqn{\pmb{\mu}} takes the form of the common mean vector \eqn{\pmb{\mu}_c}. The same is true for the t-test in the univariate case.

I have created this function illustrating this result, so this one is for educational purposes. It calculates the James test statistic, the sum of the two \eqn{T^2} test statistics, the common mean vector and the one found via numerical optimization. In the univariate case, the common mean vector is a weighted linear combination of the two sample means. So, if we take a segment connecting the two means, the common mean is somewhere on that segment.
}

\value{
A list including:
\item{tests}{
A vector with two values, the James test statistic value and the sum of the two Hotelling's test statistic using
the common mean.
}
\item{mathematics.mean}{
The common mean computed the closed form expression seen in the help file of \code{\link{james}}.
}
\item{optimised.mean}{
The common mean vector obtained from the minimisation process.
}
}

\references{
Emerson S. (2009). Small sample performance and calibration of the Empirical Likelihood method.
PhD thesis, Stanford university.

James G.S. (1954). Tests of Linear Hypothese in Univariate and Multivariate Analysis
when the Ratios of the Population Variances are Unknown. Biometrika, 41(1/2): 19--43.

}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{hotel2T2}, \link{maovjames}, \link{el.test2}, \link{eel.test2}
}
}

\examples{
james.hotel( as.matrix(iris[1:50, 1:4]), as.matrix(iris[51:100, 1:4]) )
james( as.matrix(iris[1:50, 1:4]), as.matrix(iris[51:100, 1:4]), R = 1 )
}

