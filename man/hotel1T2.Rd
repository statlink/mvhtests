\name{Hotelling's multivariate version of the 1 sample t-test for Euclidean data}
\alias{hotel1T2}
\title{
Hotelling's multivariate version of the 1 sample t-test for Euclidean data
}

\description{
Hotelling's test for testing one Euclidean population mean vector.
}

\usage{
hotel1T2(x, M, a = 0.05, R = 999, graph = FALSE)
}

\arguments{
\item{x}{
A matrix containing Euclidean data.
}
\item{a}{
The significance level, set to 0.05 by default.
}
\item{M}{
The hypothesized mean vector.
}
\item{R}{
If R is 1 no bootstrap calibration is performed and the classical p-value via the F distribution is returned. If R is greater than 1, the bootstrap p-value is returned.
}
\item{graph}{
A boolean variable which is taken into consideration only when bootstrap calibration is performed. IF TRUE the histogram of the bootstrap test statistic values is plotted.
}
}

\details{
The hypothesis test is that a mean vector is equal to some specified vector \eqn{H_0:\pmb{\mu}=\pmb{\mu}_0}. We assume that \eqn{\pmb{\Sigma}} is unknown. The first approach to this hypothesis test is parametrically, using the Hotelling's \eqn{T^2} test Mardia, Bibby and Kent (1979, pg. 125-126). The test statistic is given by
\deqn{
T^2=\frac{\left(n-p\right)n}{\left(n-1\right)p}\left(\bar{{\bf X}}-\pmb{\mu}\right)^T{\bf S}^{-1}\left(\bar{{\bf X}}-\pmb{\mu} \right).
}
Under the null hypothesis, the above test statistic follows the \eqn{F_{p,n-p}} distribution. The bootstrap version of the one-sample multivariate generalization of the simple t-test is also included in the function. An extra argument (R) indicates whether bootstrap calibration should be used or not. If R=1, then the asymptotic theory applies, if R>1, then the bootstrap p-value will be applied and the number of re-samples is equal to R.
}

\value{
A list including:
\item{m}{
The sample mean vector.
}
\item{info}{
The test statistic, the p-value, the critical value and the degrees of freedom of the F distribution (numerator and denominator).
This is given if no bootstrap calibration is employed.
}
\item{pvalue}{
The bootstrap p-value is bootstrap is employed.
}
\item{runtime}{
The runtime of the bootstrap calibration.
}
}

\references{
Mardia K.V., Kent J.T. and Bibby J.M. (1979). Multivariate analysis. London: Academic Press.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{eel.test1}, \link{el.test1}, \link{james}, \link{hotel2T2}, \link{maov}, \link{el.test2}
}
}

\examples{
x <- matrix( rnorm( 100 * 4), ncol = 4)
hotel1T2(x, numeric(4), R = 1)
hotel1T2(x, numeric(4), R = 999, graph = TRUE)
}

\keyword{ Hypothesis test }
\keyword{ mean vector }
