\name{Empirical likelihood for a one sample mean vector hypothesis testing}
\alias{el.test1}
\title{
Empirical likelihood for a one sample mean vector hypothesis testing
}

\description{
Empirical likelihood for a one sample mean vector hypothesis testing.
}

\usage{
el.test1(x, mu, R = 1, ncores = 1, graph = FALSE)
}

\arguments{
\item{x}{
A matrix containing Euclidean data.
}
\item{mu}{
The hypothesized mean vector.
}
\item{R}{
If R is 1 no bootstrap calibration is performed and the classical p-value via the \eqn{\chi^2} distribution is returned. If R is greater than 1, the bootstrap p-value is returned.
}
\item{ncores}{
The number of cores to use, set to 1 by default.
}
\item{graph}{
A boolean variable which is taken into consideration only when bootstrap calibration is performed. IF TRUE the histogram of the bootstrap test statistic values is plotted.
}
}

\details{
The \eqn{H_0} is that \eqn{\pmb{\mu} = \pmb{\mu}_0} and the constraint imposed by EL is
\deqn{
\frac{1}{n}\sum_{i=1}^{n}\left\lbrace\left[1+\pmb{\lambda}^T\left({\bf x}_i-\pmb{\mu}_0 \right)\right]^{-1}\left({\bf x}_i-\pmb{\mu}_0\right)\right\rbrace={\bf 0},
}
where the \eqn{\pmb{\lambda}} is the Lagrangian parameter introduced to maximize the above expression. Note that the maximization of is with respect to the \eqn{\pmb{\lambda}}. The probabilities have the following form
\deqn{
p_i=\frac{1}{n}\left[1+\pmb{\lambda}^T \left({\bf x}_i-\pmb{\mu}_0 \right)\right]^{-1}.
}

The log-likelihood ratio test statistic can be written as
\deqn{
\Lambda=\sum_{i=1}^{n}\log{np_i}.
}
where \eqn{d} denotes the number of variables. Under \eqn{H_0} \eqn{\Lambda \sim \chi^2_d}, asymptotically. Alternatively the bootstrap p-value may be computed.
}

\value{
A list with the outcome of the function \bold{el.test}() from the package \bold{emplik} which includes
the -2 log-likelihood ratio, the observed P-value by chi-square approximation, the final value of Lagrange multiplier \eqn{\lambda}, the gradient at the maximum, the Hessian matrix, the weights on the observations (probabilities multiplied by the sample size) and the number of iteration performed.
In addition the runtime of the procedure is reported. In the case of bootstrap, the bootstrap p-value is also returned.
}

\references{
Jing B.Y. and Wood A.T.A. (1996). Exponential empirical likelihood is
not Bartlett correctable. Annals of Statistics, 24(1): 365--369.

Owen A. (1990). Empirical likelihood ratio confidence regions.
Annals of Statistics, 18(1): 90--120.

Owen A.B. (2001). Empirical likelihood. Chapman and Hall/CRC Press.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{eel.test1}, \link{hotel1T2}, \link{james}, \link{hotel2T2}, \link{maov}, \link{el.test2}
}
}

\examples{
x <- as.matrix(iris[, 1:4])
el.test1(x, mu = numeric(4) )
eel.test1(x, mu = numeric(4) )
}

