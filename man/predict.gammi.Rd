\name{predict.gammi}
\alias{predict.gammi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predict Method for gammi Fits
}
\description{
Obtain predictions from a fit generalized additive mixed model (gammi) object. 
}
\usage{
\method{predict}{gammi}(object, 
        newx,
        newdata, 
        se.fit = FALSE,
        type = c("link", "response", "terms"),
        conf.int = FALSE, 
        conf.level = 0.95,
        ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
Object of class "gammi"
}
  \item{newx}{
  Matrix of new \code{x} scores for prediction (default S3 method). Must have \eqn{p} columns arranged in the same order as the \code{x} matrix used to fit the model.
}
  \item{newdata}{
Data frame of new data scores for prediction (S3 "formula" method). Must contain all variables in the \code{formula} (and \code{fixed} formula if applicable) used to fit the model.
}
  \item{se.fit}{
Logical indicating whether standard errors of predictions should be returned.
}
  \item{type}{
Type of prediction to return: link = linear prediction, response = fitted value, and terms = matrix where each columns contains each term's linear predictor contribution.
}
  \item{conf.int}{
Logical indicating whether confidence intervals for predictions should be returned.
}
  \item{conf.level}{
Scalar between 0 and 1 controlling the confidence level for the interval. Ignored if \code{conf.int = FALSE}.
}
  \item{...}{
Additional arugments (ignored).  
}
}
\details{
The default of \code{type = "link"} returns the model implied linear predictor corresponding to \code{newx} or \code{newdata}, i.e., 
\deqn{g(\hat{\boldsymbol\mu}_{\theta (\mathrm{new}) }) = \hat{f}_\theta(\mathbf{X}_\mathrm{new}, \mathbf{Z}_\mathrm{new}) + \mathbf{X}_\mathrm{new}^\top \hat{\boldsymbol\beta}_\theta}
where \eqn{\hat{f}_\theta(\cdot)} is the estimated smooth function (with the subscript of \eqn{\theta} denoting the dependence on the variance parameters), and \eqn{\hat{\boldsymbol\beta}_\theta} are the fixed effect estimates (if applicable). Note that \eqn{\mathbf{X}_\mathrm{new}} and \eqn{\mathbf{Z}_\mathrm{new}} denote the new data at which the predictions will be formed.

Using \code{type = "response"} returns the predictions on the fitted value scale, i.e.,
\deqn{\hat{\boldsymbol\mu}_{\theta (\mathrm{new})} = g^{-1} \left( \hat{f}_\theta(\mathbf{X}_\mathrm{new}, \mathbf{Z}_\mathrm{new}) + \mathbf{X}_\mathrm{new}^\top \hat{\boldsymbol\beta}_\theta \right)}
where \eqn{g^{-1}(\cdot)} denotes the inverse of the chosen link function.

Using \code{type = "terms"} returns a matrix where each column contains the linear predictor contribution for a different model term, i.e., the \eqn{k}-th column contains
\deqn{\hat{f}_{\theta k}(\mathbf{X}_\mathrm{new}, \mathbf{Z}_\mathrm{new}) + \mathbf{X}_{\mathrm{new} k }^\top \hat{\boldsymbol\beta}_{\theta k}}
where \eqn{\hat{f}_{\theta k}} is the \eqn{k}-th additive function, i.e., \eqn{\hat{f}_\theta(\mathbf{X}_\mathrm{new}, \mathbf{Z}_\mathrm{new}) = \sum_{k=1}^K \hat{f}_{\theta k}(\mathbf{X}_\mathrm{new}, \mathbf{Z}_\mathrm{new})} and the second term denotes the (optional) fixed-effect contribution for the \eqn{k}-th term, i.e., \eqn{\mathbf{X}_\mathrm{new}^\top \hat{\boldsymbol\beta}_\theta = \sum_{k=1}^K \mathbf{X}_{\mathrm{new} k }^\top \hat{\boldsymbol\beta}_{\theta k}}
}

\value{
If \code{type = "link"} or \code{type = "response"}, returns either a vector (of predictions corresponding to the new data) or a data frame that contains the predictions, along with their standard errors and/or confidence interval endpoints (as controlled by \code{se.fit} and \code{conf.int} arguments).

If \code{type = "terms"}, returns either a matrix (with columns containing predictions for each term) or a list that contains the term-wise predictions, along with their standard errors and/or confidence interval endpoints (as controlled by \code{se.fit} and \code{conf.int} arguments).
}
\references{
Helwig, N. E. (2024). Precise tensor product smoothing via spectral splines. \emph{Stats, 7}(1), 34-53, \doi{10.3390/stats7010003}
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\note{
Terms entered through the \code{random} argument of the \code{\link{gammi}} function are \bold{not} included as a part of predictions.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{gammi}} for fitting generalized additive mixed models

\code{\link{plot.gammi}} for plotting effects from \code{gammi} objects

\code{\link{summary.gammi}} for summarizing results from \code{gammi} objects
}
\examples{
# load 'gammi' package
library(gammi)

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

# get fitted values via predict
fit <- predict(mod, newdata = data.frame(x = x, z = z))
max(abs(fit - mod$fitted.values))

# get fitted values with SE and CI
fit <- predict(mod, newdata = data.frame(x = x, z = z), conf.int = TRUE)
head(fit)

# get fitted values with SE and CI for each term
fit <- predict(mod, newdata = data.frame(x = x, z = z), 
               type = "terms", conf.int = TRUE)
str(fit)                                  # list with 4 components
head(sapply(fit, function(x) x[,1]))      # for x effect
head(sapply(fit, function(x) x[,2]))      # for z effect
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