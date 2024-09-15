\name{gammi}
\alias{gammi}
\alias{gammi.default}
\alias{gammi.formula}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit a Generalized Additive Mixed Model
}
\description{
Fits generalized additive models (GAMs) and generalized additive mixed model (GAMMs) using \bold{lme4} as the tuning engine. Predictor groups can be manually input (default S3 method) or inferred from the model (S3 "formula" method). Smoothing parameters are treated as variance components and estimated using REML/ML (gaussian) or Laplace approximation to ML (others).
}
\usage{
gammi(x, \dots)

\method{gammi}{default}(x,
      y,
      group,
      family = gaussian,
      fixed = NULL,
      random = NULL,
      data = NULL,
      REML = TRUE,
      control = NULL,
      start = NULL,
      verbose = 0L,
      nAGQ = 10L,
      subset, 
      weights, 
      na.action, 
      offset, 
      mustart, 
      etastart,
      ...)

\method{gammi}{formula}(formula, 
      data, 
      family = gaussian,
      fixed = NULL, 
      random = NULL, 
      REML = TRUE,
      control = NULL,
      start = NULL,
      verbose = 0L,
      nAGQ = 10L,
      subset, 
      weights, 
      na.action, 
      offset, 
      mustart, 
      etastart,
      ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Model (design) matrix of dimension \code{nobs} by \code{nvars} (\eqn{n \times p}).
}
  \item{y}{
Response vector of length \eqn{n}. 
}
  \item{group}{
Group label vector (factor, character, or integer) of length \eqn{p}. Predictors with the same label are assumed to have the same variance parameter.
}
\item{formula}{
Model formula: a symbolic description of the model to be fitted. Uses the same syntax as \code{\link{lm}} and \code{\link{glm}}. 
}
  \item{family}{
Assumed exponential \code{\link{family}} (and link function) for the response variable.  
}
  \item{fixed}{
For default method: a character vector specifying which \code{group} labels should be treated as fixed effects. For formula method: a one-sided formula specifying the fixed effects model structure. 
}
  \item{random}{
A one-sided formula specifying the random effects structure using \bold{lme4} syntax. See Note.  
}
  \item{data}{
Optional data frame containing the variables referenced in \code{formula}, \code{fixed}, and/or \code{random}.
}
  \item{REML}{
Logical indicating whether REML versus ML should be used to tune the smoothing parameters and variance components.
}
  \item{control}{
List containing the control parameters (output from \code{\link[lme4]{lmerControl}} or \code{\link[lme4]{glmerControl}}).
}
  \item{start}{
List (with names) of starting parameter values for model parameters.
}
  \item{verbose}{
Postive integer that controls the level of output displayed during optimization.
}
  \item{nAGQ}{
Numer of adaptive Gaussian quadrature points. Only used for non-Gaussian responses with a single variance component.
}
  \item{subset}{
Optional expression indicating the subset of rows to use for the fitting (defaults to all rows).
}
  \item{weights}{
Optional vector indicating prior observations weights for the fitting (defaults to all ones).
}
  \item{na.action}{
Function that indicates how \code{NA} data should be dealt with. Default (of \code{na.omit}) will omit any observations with missing data on any variable.
}
  \item{offset}{
Optional vector indicating each observation's offset for the fitting (defaults to all zeros).
}
  \item{mustart}{
Optional starting values for the mean (fitted values).
}
  \item{etastart}{
Optional starting values for the linear predictors.
}
  \item{...}{
Optional arguments passed to the \code{\link{spline.model.matrix}} function, e.g., spline \code{knots} or \code{df} for each term.
}
}
\details{
Fits a generalized additive mixed model (GAMM) of the form
\deqn{
g(\mu) = f(\mathbf{X}, \mathbf{Z}) + \mathbf{X}^\top \boldsymbol\beta + \mathbf{Z}^\top \boldsymbol\alpha
}
where 
\itemize{
  \item \eqn{\mu = E(Y | \mathbf{X}, \mathbf{Z})} is the conditional expectation of the response \eqn{Y} given the predictor vectors \eqn{\mathbf{X} = (X_1, \ldots, X_p)^\top} and \eqn{\mathbf{Z} = (Z_1, \ldots, Z_q)^\top}
  
  \item the function \eqn{g(\cdot)} is a user-specified (invertible) link function
  
  \item the function \eqn{f(\cdot)} is an unknown smooth function of the predictors (specified by \code{formula})
  
  \item the vector \eqn{\mathbf{X}} is the fixed effects component of the design (specified by \code{fixed})
  
  \item the vector \eqn{\mathbf{Z}} is the random effects component of the design (specified by \code{random})
  
  \item the vector \eqn{\boldsymbol\beta} contains the unknown fixed effects coefficients
  
  \item the vector \eqn{\boldsymbol\alpha} contains the unknown Gaussian random effects
  
}

Note that the mean function \eqn{f(\cdot)} can include main and/or interaction effects between any number of predictors. Furthermore, note that the fixed effects in \eqn{\mathbf{X}^\top \boldsymbol\beta} and the random effects in \eqn{\mathbf{Z}^\top \boldsymbol\alpha} are both optional.

}
\value{
An object of class \code{"gammi"} with the following elements:
\item{fitted.values }{model predictions on the data scale}
\item{linear.predictors }{model predictions on the link scale}
\item{coefficients}{coefficients used to make the predictions}
\item{random.coefficients}{coefficients corresponding to the \code{random} argument, i.e., the BLUPs.}
\item{term.labels}{labels for the terms included in the \code{coefficients}}
\item{dispersion}{estimated dispersion parameter = \code{deviance/df.residual} when \code{is.null(random)}}
\item{vcovchol}{Cholesky factor of covariance matrix such that \code{tcrossprod(vcovchol)} gives the covariance matrix for the combined coefficient vector  \code{c(coefficients, random.coefficients)}}
\item{family}{exponential family distribution (same as input)}
\item{logLik}{log-likelihood for the solution}
\item{aic}{AIC for the solution}
\item{deviance}{model deviance, i.e., two times the negative log-likelihood}
\item{null.deviance}{deviance of the null model, i.e., intercept only. Will be \code{NA} if the \code{random} argument is used.}
\item{r.squared}{proportion of null deviance explained = \code{1 - deviance/null.deviance}. Will be \code{NA} if the \code{random} argument is used; see Note.}
\item{nobs}{number of observations used in fit}
\item{leverages}{leverage scores for each observation}
\item{edf}{effective degrees of freedom = \code{sum(leverages)}}
\item{df.random}{degress of freedom corresponding to \code{random} formula, i.e., number of co/variance parameters}
\item{df.residual}{residual degrees of freedom = \code{nobs - edf}}
\item{x}{input \code{x} matrix (default method only)}
\item{group}{character vector indicating which columns of \code{x} belong to which model terms}
\item{scale}{numeric vector giving the scale parameter used to z-score each term's data}
\item{fixed}{fixed effects terms (default method) or formula (formula method); will be \code{NULL} if no fixed terms are included}
\item{random}{random effects formula}
\item{mer}{object of class \code{"merMod"}, such as output by \link[lme4]{lmer}, with model fit information on a standardized scale}
\item{VarCorr}{data frame with variance and covariance parameter estimates from \code{mer} transformed back to the original scale}
\item{call}{function call}
\item{data}{input data}
\item{contrasts}{list of contrasts applied to \code{fixed} terms; will be \code{NULL} if no fixed terms are included}
\item{spline.info}{list of spline parameters for terms in \code{x} or \code{formula}}
\item{formula}{input model formula}
}
\references{
Bates, D., Maechler, M., Bolker, B., & Walker, S. (2015). Fitting linear mixed-effects models using lme4. \emph{Journal of Statistical Software, 67}(1), 1–48. \doi{10.18637/jss.v067.i01}

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

\section{Random Syntax}{
The \code{random} argument uses standard \link[lme4]{lmer} syntax:
\itemize{
\item \code{(1 | g)} for a random intercept for each level of \code{g}
\item \code{(1 | g1) + (1 | g2)} for random intercepts for g1 and g2
\item \code{(1 | g1/g2) = (1 | g1) + (1 | g1:g2)} for random intercepts for g1 and g2 nested within g1
\item \code{(x | g) = (1 + x | g)} for a correlated random intercept and slope of x for each level of \code{g}
\item \code{(x || g) = (1 | g) + (0 + x | g)} for an uncorrelated random intercept and slope of x for each level of \code{g}
}
}

\section{Warning}{
For stable computation, any terms entered through \code{x} (default method) or \code{formula} and/or \code{fixed} (formula method) are z-scored prior to fitting the model. Note that terms entered through \code{random} are not standardized. 

The \code{"mer"} component of the output contains the model fitting results for a z-scored version of the original data (i.e., this fit is on a different scale). Consequently, the \code{"mer"} component should \bold{not} be used for prediction and/or inference purposes. All prediction and inference should be conducted using the plot, predict, and summary methods mentioned in the `See Also' section.

The \code{"VarCorr"} component contains the estimated variance/covariance parameters transformed back to the original scale.
}

\seealso{
\code{\link{plot.gammi}} for plotting effects from \code{gammi} objects

\code{\link{predict.gammi}} for predicting from \code{gammi} objects

\code{\link{summary.gammi}} for summarizing results from \code{gammi} objects
}
\examples{
##############***##############   EXAM EXAMPLE   ##############***##############

# load 'exam' help file
?exam

# load data
data(exam)

# header of data
head(exam)

# fit model
mod <- gammi(Exam.score ~ VRQ.score, data = exam,
             random = ~ (1 | Primary.school) + (1 | Secondary.school))
             
# plot results
plot(mod)

# summarize results
summary(mod)




\donttest{
#############***#############   GAUSSIAN EXAMPLE   #############***#############

#~~~Example 1:  Single Predictor ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# generate data
n <- 1000
x <- seq(0, 1, length.out = n)
fx <- sin(2 * pi * x)
set.seed(1)
y <- fx + rnorm(n)

# fit model via formula method
mod <- gammi(y ~ x)
mod

# fit model via default method
modmat <- spline.model.matrix(y ~ 0 + x)
tlabels <- attr(modmat, "term.labels")
tassign <- attr(modmat, "assign")
g <- factor(tlabels[tassign], levels = tlabels)
mod0 <- gammi(modmat, y, g)
mod0

# summarize fit model
summary(mod)

# plot function estimate
plot(mod)



#~~~Example 2:  Additive Model ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# mean function
eta <- function(x, z, additive = TRUE){
  mx1 <- cos(2 * pi * (x - pi))
  mx2 <- 30 * (z - 0.6)^5
  mx12 <- 0
  if(!additive) mx12 <- sin(pi * (x - z))
  mx1 + mx2 + mx12
}

# generate data
set.seed(1)
n <- 1000
x <- runif(n)
z <- runif(n)
fx <- eta(x, z)
y <- fx + rnorm(n)

# fit model via formula method
mod <- gammi(y ~ x + z)
mod

# fit model via default method
modmat <- spline.model.matrix(y ~ 0 + x + z)
tlabels <- attr(modmat, "term.labels")
tassign <- attr(modmat, "assign")
g <- factor(tlabels[tassign], levels = tlabels)
mod0 <- gammi(modmat, y, g)
mod0

# summarize fit model
summary(mod)

# plot function estimate
plot(mod)



#~~~Example 3:  Interaction Model ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# mean function
eta <- function(x, z, additive = TRUE){
  mx1 <- cos(2 * pi * (x - pi))
  mx2 <- 30 * (z - 0.6)^5
  mx12 <- 0
  if(!additive) mx12 <- sin(pi * (x - z))
  mx1 + mx2 + mx12
}

# generate data
set.seed(1)
n <- 1000
x <- runif(n)
z <- runif(n)
fx <- eta(x, z, additive = FALSE)
y <- fx + rnorm(n)

# fit model via formula method
mod <- gammi(y ~ x * z)
mod

# fit model via default method
modmat <- spline.model.matrix(y ~ 0 + x * z)
tlabels <- attr(modmat, "term.labels")
tassign <- attr(modmat, "assign")
g <- factor(tlabels[tassign], levels = tlabels)
mod0 <- gammi(modmat, y, g)
mod0

# summarize fit model
summary(mod)

# plot function estimate
plot(mod)



#~~~Example 4:  Random Intercept ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# mean function
eta <- function(x, z, additive = TRUE){
  mx1 <- cos(2 * pi * (x - pi))
  mx2 <- 30 * (z - 0.6)^5
  mx12 <- 0
  if(!additive) mx12 <- sin(pi * (x - z))
  mx1 + mx2 + mx12
}

# generate mean function
set.seed(1)
n <- 1000
nsub <- 50
x <- runif(n)
z <- runif(n)
fx <- eta(x, z)

# generate random intercepts
subid <- factor(rep(paste0("sub", 1:nsub), n / nsub),
                levels = paste0("sub", 1:nsub))
u <- rnorm(nsub, sd = sqrt(1/2))

# generate responses
y <- fx + u[subid] + rnorm(n, sd = sqrt(1/2))

# fit model via formula method
mod <- gammi(y ~ x + z, random = ~ (1 | subid))
mod

# fit model via default method
modmat <- spline.model.matrix(y ~ 0 + x + z)
tlabels <- attr(modmat, "term.labels")
tassign <- attr(modmat, "assign")
g <- factor(tlabels[tassign], levels = tlabels)
mod0 <- gammi(modmat, y, g, random = ~ (1 | subid))
mod0

# summarize fit model
summary(mod)

# plot function estimate
plot(mod)




#############***#############   BINOMIAL EXAMPLE   #############***#############

#~~~Example 1:  Single Predictor ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# generate data
n <- 1000
x <- seq(0, 1, length.out = n)
fx <- sin(2 * pi * x)
set.seed(1)
y <- rbinom(n = n, size = 1, prob = 1 / (1 + exp(-fx)))

# fit model
mod <- gammi(y ~ x, family = binomial)
mod

# summarize fit model
summary(mod)

# plot function estimate
plot(mod)



#~~~Example 2:  Additive Model ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# mean function
eta <- function(x, z, additive = TRUE){
  mx1 <- cos(2 * pi * (x - pi))
  mx2 <- 30 * (z - 0.6)^5
  mx12 <- 0
  if(!additive) mx12 <- sin(pi * (x - z))
  mx1 + mx2 + mx12
}

# generate data
set.seed(1)
n <- 1000
x <- runif(n)
z <- runif(n)
fx <- 1 + eta(x, z)
y <- rbinom(n = n, size = 1, prob = 1 / (1 + exp(-fx)))

# fit model
mod <- gammi(y ~ x + z, family = binomial)
mod

# summarize fit model
summary(mod)

# plot function estimate
plot(mod)



#~~~Example 3:  Interaction Model ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# mean function
eta <- function(x, z, additive = TRUE){
  mx1 <- cos(2 * pi * (x - pi))
  mx2 <- 30 * (z - 0.6)^5
  mx12 <- 0
  if(!additive) mx12 <- sin(pi * (x - z))
  mx1 + mx2 + mx12
}

# generate data
set.seed(1)
n <- 1000
x <- runif(n)
z <- runif(n)
fx <- eta(x, z, additive = FALSE)
y <- rbinom(n = n, size = 1, prob = 1 / (1 + exp(-fx)))

# fit model
mod <- gammi(y ~ x * z, family = binomial)
mod

# summarize fit model
summary(mod)

# plot function estimate
plot(mod)



#~~~Example 4:  Random Intercept ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# mean function
eta <- function(x, z, additive = TRUE){
  mx1 <- cos(2 * pi * (x - pi))
  mx2 <- 30 * (z - 0.6)^5
  mx12 <- 0
  if(!additive) mx12 <- sin(pi * (x - z))
  mx1 + mx2 + mx12
}

# generate mean function
set.seed(1)
n <- 1000
nsub <- 50
x <- runif(n)
z <- runif(n)
fx <- 1 + eta(x, z)

# generate random intercepts
subid <- factor(rep(paste0("sub", 1:nsub), n / nsub),
                levels = paste0("sub", 1:nsub))
u <- rnorm(nsub, sd = sqrt(1/2))

# generate responses
y <- rbinom(n = n, size = 1, prob = 1 / (1 + exp(-(fx+u[subid]))))

# fit model
mod <- gammi(y ~ x + z, random = ~ (1 | subid), family = binomial)
mod

# summarize fit model
summary(mod)

# plot function estimate
plot(mod)
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{ regression }
\keyword{ smooth }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
