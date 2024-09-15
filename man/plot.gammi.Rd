\name{plot.gammi}
\alias{plot.gammi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot Method for gammi Fits
}
\description{
Plots main and interaction effects from a fit gammi object.
}
\usage{
\method{plot}{gammi}(x, terms = x$term.labels, conf.int = TRUE, n = 400,
     intercept = FALSE, random = TRUE, ask = dev.interactive(), 
     xlab = NULL, ylab = NULL, zlab = NULL, main = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Object of class "gammi"
}
  \item{terms}{
Which model term(s) should be plotted? Default plots all terms.
}
  \item{conf.int}{
Should a 95\% confidence interval be added to the plot(s)?
}
  \item{n}{
Number of points used to plot each of the (continuous) terms.
}
  \item{intercept}{
Should the intercept be added to the y-axis of the plot(s)?
}
  \item{random}{
Should Q-Q plots of the random coefficients be produced?  
}
  \item{ask}{
Should the user be asked before each plot is produced?  
}
  \item{xlab}{
Optional x-axis label for plot(s).
}
  \item{ylab}{
Optional y-axis label for plot(s).
}
  \item{zlab}{
Optional z-axis label for plot(s).
}
  \item{main}{
Optional title for plot(s).
}
  \item{...}{
Additional arguments passed to internal plotting functions.  
}
}
\details{
Default use plots each effect function along with a 95\% confidence interval (if applicable). Line plots are used for continuous predictors, bar plots are used for categorical predictors, Q-Q plots are used for random effects, and image plots are used for two-way interactions. The \code{visualizer1} and \code{visualizer2} functions are used to plot main and interaction effects, respectively.
}
\value{
A plot is produced and nothing is returned.
}
\references{
Helwig, N. E. (2024). Precise tensor product smoothing via spectral splines. \emph{Stats, 7}(1), 34-53, \doi{10.3390/stats7010003}
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\note{
Three-way and higher-order interactions are not currently supported. 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{gammi}} for fitting generalized additive mixed models

\code{\link{predict.gammi}} for predicting from \code{gammi} objects

\code{\link{summary.gammi}} for summarizing results from \code{gammi} objects
}
\examples{
# load data
data(exam)

# header of data
head(exam)

# fit model
mod <- gammi(Exam.score ~ VRQ.score, data = exam,
             random = ~ (1 | Primary.school) + (1 | Secondary.school))
       
# plot terms
plot(mod)

# refit model with Secondary.school as penalized nominal effect
mod <- gammi(Exam.score ~ Secondary.school + VRQ.score, data = exam,
             random = ~ (1 | Primary.school))
             
# plot terms
plot(mod)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{ hplot }
\keyword{ regression }
\keyword{ smooth }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.