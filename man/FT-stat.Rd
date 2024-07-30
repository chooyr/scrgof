\name{FT_stat}
\alias{scrgof}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Freeman-Tukey test
}
\description{
	A function for calculating Freeman-Tukey fit statistics from observed and expected values.
}
\usage{
FT_stat(observed, expected)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{observed}{
		observed values from data or simulations
}
  \item{expected}{
		expected values from the fitted model
}
}
\details{
	{FT_stat} calculates the Freeman-Tukey statistic: \eqn{FT-stat = \displaystyle\sum\limits_{n=1}^{N} \left(\displaystyle\sqrt{O_n} - \displaystyle\sqrt{E_n} )^2}. 
}

\author{
Yan Ru Choo
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
