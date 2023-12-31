\name{Exponential empirical likelihood for a one sample mean vector hypothesis testing}
\alias{eel.test1}
\title{
Exponential empirical likelihood for a one sample mean vector hypothesis testing
}

\description{
Exponential empirical likelihood for a one sample mean vector hypothesis testing.
}

\usage{
eel.test1(x, mu, tol = 1e-06, R = 1)
}

\arguments{
\item{x}{
A matrix containing Euclidean data.
}
\item{mu}{
The hypothesized mean vector.
}
\item{tol}{
The tolerance value used to stop the Newton-Raphson algorithm.
}
\item{R}{
The number of bootstrap samples used to calculate the p-value. If R = 1 (default value), no bootstrap calibration is performed
}
}

\details{
Exponential empirical likelihood or exponential tilting was first introduced by Efron (1981) as a way to perform a "tilted" version of the bootstrap for the one sample mean hypothesis testing. Similarly to the empirical likelihood, positive weights \eqn{p_i}, which sum to one, are allocated to the observations, such that the weighted sample mean \eqn{{\bf \bar{x}}} is equal to some population mean \eqn{\pmb{\mu}_0}, under the \eqn{H_0}. Under \eqn{H_1} the weights are equal to \eqn{\frac{1}{n}}, where \eqn{n} is the sample size. Following Efron (1981), the choice of \eqn{p_is} will minimize the Kullback-Leibler distance from \eqn{H_0} to \eqn{H_1}
\deqn{
D\left(L_0,L_1\right)=\sum_{i=1}^np_i\log\left(np_i\right),
}
subject to the constraint \eqn{\sum_{i=1}^np_i{\bf x}_i=\pmb{\mu}_0}. The probabilities take the form
\deqn{
p_i=\frac{e^{\pmb{\lambda}^T{\bf x}_i}}{\sum_{j=1}^ne^{\pmb{\lambda}^T{\bf x}_j}}
}
and the constraint becomes
\deqn{
\frac{\sum_{i=1}^ne^{\pmb{\lambda}^T{\bf x}_i}\left({\bf x}_i-\pmb{\mu}_0\right)}{\sum_{j=1}^ne^{\pmb{\lambda}^T{\bf x}_j}}=0 \Rightarrow \frac{\sum_{i=1}^n{\bf x}_ie^{\pmb{\lambda}^T{\bf x}_i}}{\sum_{j=1}^ne^{\pmb{\lambda}^T{\bf x}_j}}-\pmb{\mu}_0={\bf 0}.
}
A numerical search over \eqn{\pmb{\lambda}} is required. Under \eqn{H_0} \eqn{\Lambda \sim \chi^2_d}, where \eqn{d} denotes the number of variables. Alternatively the bootstrap p-value may be computed.
}

\value{
A list including:
\item{p}{
The estimated probabilities.
}
\item{lambda}{
The value of the Lagrangian parameter \eqn{\lambda}.
}
\item{iter}{
The number of iterations required by the newton-Raphson algorithm.
}
\item{info}{
The value of the log-likelihood ratio test statistic along with its corresponding p-value.
}
\item{runtime}{
The runtime of the process.
}
}

\references{
Efron B. (1981) Nonparametric standard errors and confidence intervals. Canadian Journal of
Statistics, 9(2): 139--158.

Jing B.Y. and Wood A.T.A. (1996). Exponential empirical likelihood is not
Bartlett correctable. Annals of Statistics, 24(1): 365--369.

Owen A. B. (2001). Empirical likelihood. Chapman and Hall/CRC Press.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{el.test1}, \link{hotel1T2}, \link{james}, \link{hotel2T2}, \link{maov}, \link{el.test2}
}
}

\examples{
x <- as.matrix( iris[, 1:4] )
eel.test1(x, numeric(4) )
el.test1(x, numeric(4) )
}

\keyword{ Hypothesis test }
\keyword{ mean vector }
