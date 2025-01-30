\name{summary.gammi}
\alias{summary.gammi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summary Method for gammi Fits
}
\description{
Obtain summary statistics from a fit generalized additive mixed model (gammi) object.
}
\usage{
\method{summary}{gammi}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
Object of class "gammi"
}
  \item{...}{
Additional arguments (currently ignored)  
}
}
\details{
Produces significance testing and model diagnostic information. The significance tests use the Bayesian interpretation of a smoothing spline. The variable importance indices sum to 100 but can be negative for some terms. The variance inflation factors should ideally be 1 for all terms; values greater than 5 or 10 can indicate noteworthy multicollinearity.
}
\value{
An object of class "summary.gammi", which is a list with components:
  \item{call }{the model call, i.e., \code{object$call}}
  \item{term.labels }{the model term labels (character vector)}
  \item{family }{the exponential \code{\link[stats]{family}} object}
  \item{logLik}{log-likelihood for the solution}
  \item{aic}{AIC for the solution}
  \item{deviance }{the model deviance (numeric)}
  \item{deviance.resid }{the deviance residuals}
  \item{r.squared }{the model R-squared (numeric); see Note}
  \item{df}{the total degrees of freedom = \code{object$edf + object$df.random}}
  \item{significance }{the signififance testing information (matrix)}
  \item{importance }{the variable importance information (numeric)}
  \item{vif }{the variance inflation factors (numeric)}
}
\references{
Helwig, N. E. (2024). Precise tensor product smoothing via spectral splines. \emph{Stats, 7}(1), 34-53, \doi{10.3390/stats7010003}
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\note{
The model R-squared is the proportion of the null deviance that is explained by the model, i.e.,

\code{r.squared = 1 - deviance / null.deviance}

where \code{deviance} is the deviance of the model, and \code{null.deviance} is the deviance of the null model. 

When the \code{random} argument is used, \code{null.deviance} and \code{r.squared} will be \code{NA}. This is because there is not an obvious null model when random effects are included, e.g., should the null model include or exclude the random effects? Assuming that is it possible to define a reasonable \code{null.deviance} in such cases, the above formula can be applied to calculate the model R-squared for models that contain random effects.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{gammi}} for fitting generalized additive mixed models

\code{\link{plot.gammi}} for plotting effects from \code{gammi} objects

\code{\link{predict.gammi}} for predicting from \code{gammi} objects
}
\examples{
# load 'gammi' package
library(gammi)

# load data
data(exam)

# header of data
head(exam)

# fit model
mod <- gammi(Exam.score ~ VRQ.score, data = exam,
             random = ~ (1 | Primary.school) + (1 | Secondary.school))
       
# summarize results
summary(mod)

# refit model with Secondary.school as penalized nominal effect
mod <- gammi(Exam.score ~ Secondary.school + VRQ.score, data = exam,
             random = ~ (1 | Primary.school))
             
# summarize results
summary(mod)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{ htest }
\keyword{ regression }
\keyword{ smooth }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
