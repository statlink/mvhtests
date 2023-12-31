\name{Repeated measures ANOVA (univariate data) using Hotelling's T2 test}
\alias{rm.hotel}
\title{
Repeated measures ANOVA (univariate data) using Hotelling's \eqn{T^2} test
}

\description{
Repeated measures ANOVA (univariate data) using Hotelling's \eqn{T^2} test.
}

\usage{
rm.hotel(x, a = 0.05)
}

\arguments{
\item{x}{
A numerical matrix with the repeated measurements. Each column contains the
values of the repeated measurements.
}
\item{a}{
The level of significance, default value is equal to 0.05.
}
}

\details{
We now show how can one use Hotelling's \eqn{T^2} test to analyse univariate repeated measures. Univariate analysis of variance for repeated measures is the classical way, but we can use this multivariate test as well. In the repeated measures ANOVA case, we have many repeated observations from the same \eqn{n} subjects, usually at different time points and the interest is to see whether the means of the samples are equal or not \eqn{\mu_1=\mu_2=\ldots=\mu_k} assuming \eqn{k} repeated measurements. We can of course change this null hypothesis and test many combinations of means. The idea in any case is to construct a matrix of contrasts. I will focus here in the first case only and in particular the null hypothesis and the matrix of contrasts \eqn{\bf C} are
\deqn{
\left( {\begin{array}{c}
\mu_1=\mu_2 \\
\mu_2=\mu_3 \\
\vdots  \\
\mu_{k-1}=\mu_k \end{array}} \right)=
\left( {\begin{array}{ccccc}
1 & -1 & 0 & \ldots & 0 \\
1 &  0 & -1 & \dots & 0 \\
\vdots & \vdots & \vdots & \vdots & \vdots \\
1 & 0 & 0 & \ldots & -1 \\
\end{array}} \right)\pmb{\mu}={\bf C}\pmb{\mu}.
}

The contrast matrix \eqn{\bf C} has \eqn{k-1} independent rows and if there is no treatment effect, \eqn{{\bf C}\pmb{\mu}={\bf 0}}.

The test statistic is
\deqn{
T_r^2=\frac{\left(n-k+1\right)}{\left(n-1\right)\left(k-1\right)}n\left({\bf C}\bar{\bf x}\right)^T
\left({\bf CSC}^T\right)^{-1}\left({\bf C}\bar{\bf x}\right) \sim F_{k-1,n-k+1}.
}
}

\value{
A list including:
\item{m}{
The mean vector.
}
\item{result}{
A vector with the test statistic value, it's associated p-value, the numerator and denominator
degrees of freedom and the critical value.
}
}

%\references{
%}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{maov}, \link{hotel2T2}, \link{james}
}
}

\examples{
x <- as.matrix(iris[, 1:4]) ## assume they are repeated measurements
rm.hotel(x)
}

