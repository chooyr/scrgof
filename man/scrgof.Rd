\name{scrgof}
\alias{scrgof}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Goodness-of-Fit Tests For Spatial Capture-Recapture
}
\description{
	Tests goodness-of-fit for SECR models, using Monte Carlo resampling methods to obtain the empirical distributions of fit statistics conditioned on realised activity centres.
}
\usage{
scrgof(object, nsamp = 1000, maxperpoly = 100, renumber = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
		a fitted SECR model
}
  \item{nsamp}{
		number of Monte Carlo resamples to be drawn from the model
}
}
\details{
	{scrgof} calculates fit using summary statistics suggested in Royle et al. (2014).
	FT-ind-trap is a fit statistic based on the number of encounters of each individual at each trap.  FT-individuals is based on the total encounters of each individual. FT-traps is based on the total encounters at each trap.  These summary statistics are conditioned on the predicted activity centre locations of observed and hypothetical animals in the population.  Fit statistics are then calculated from the summary statistics in the form of Freeman-Tukey statistics: \eqn{FT-stat = \displaystyle\sum\limits_{n=1}^{N} \left(\displaystyle\sqrt{O_n} - \displaystyle\sqrt{E_n} )^2}.  The uncertainty of the model parameters and the activity centre locations are propagated to the fit statistics using Monte Carlo resampling methods.
	
	{scrgof} is also a wrapper for {secr::sim.detect}, which is used to calculate expected values and simulate new capture histories when resampling from the fitted model, and is subjected to the same limitations i.e. modles fitted using conditional likelihood, telemetry models and exotic behavioural response models.
}
\value{
  \item{scrgof_pval}{
		list of p-values for each GoF test
	}
	\item{gof_ik}{
		data frame of FT-ind-trap fit statistics for simulated and actual observations
	}
	\item{gof_i}{
		data frame of FT-individuals fit statistics for simulated and actual observations
	}
	\item{gof_j}{
		data frame of FT-traps fit statistics for simulated and actual observations
	}
	
}
\references{
Royle, J. A., Chandler, R. B., Sollmann, R., and Gardner, B. (2014b). Chapter 8 - Model Selection and Assessment. In Royle, J. A., Chandler, R. B., Sollmann, R., and Gardner, B., editors, Spatial Capture-recapture, pages 219–243. Academic Press, Boston.
}
\author{
Yan Ru Choo
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# Analysis of wolverine data from Royle et al (2014)

library(secr)
library(scrgof)

data(bears)

# Fitting capture history to a null model
bear0 <- secr.fit(bearCH, buffer=1e5, detectfn=0)
summary(bear0)
# test GoF of bear0
bear0_GoF <- scrgof(bear0, nsamp = 1000)
# return p-values
bear0_GoF$scrgof_pval

# Fitting capture history to a Mb model
bearh2 <- secr.fit(bearCH, buffer=1e5, model=list(g0 ~ h2, sigma ~ h2))
summary(bearh2)
# test GoF of bearh2
bearh2_gof <- scrgof(bearh2, nsamp = 1000)
# return p-values
bearh2_GoF$scrgof_pval

