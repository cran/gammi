\name{visualizers}
\alias{visualizer1}
\alias{visualizer2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Internal Functions for Plot Method
}
\description{
Internal functions used by the \code{\link{plot.gammi}} function to visualize main effects and two-way interaction effects in fit \code{gammi} objects.
}
\usage{
visualizer1(x, y, bars = FALSE, bw = 0.02, lty = 1, lwd = 2, col = "black", 
            lwr = NULL, upr = NULL, ci.lty = 2, ci.lwd = 1.25, ci.col = "black",
            zero = TRUE, zero.lty = 3, xlim = NULL, ylim = NULL, 
            xlab = NULL, ylab = NULL, main = NULL, add = FALSE, ...)
            
visualizer2(x, y, z, col = NULL, ncolor = 21,
            xlim = NULL, ylim = NULL, zlim = NULL, zline = 1.5,
            xlab = NULL, ylab = NULL, zlab = NULL, main = NULL, 
            xticks = NULL, xlabels = NULL, yticks = NULL, ylabels = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, y, z}{
For 1D plots: \code{x} and \code{y} are the primary inputs to the \code{\link{plot}} function. For 2D plots: these are the primary inputs to the \code{\link{image}} function.
}
  \item{bars}{
For 1D plots: logical indicating whether to create a line plot (default) or a bar plot (\code{bars = TRUE}).
}
  \item{bw}{
For 1D plots: width of the bars relative to range of \code{x} (ignored if \code{bars = FALSE}).
}
  \item{lty, lwd}{
For 1D plots: line type and width for 1D plots.
}
  \item{col}{
For 1D plots: single color for line/bar plot. For 2D plots: vector of colors for image plot.
}
  \item{ncolor}{
For 2D plots: number of colors used for image plot and color legend, see Note.
}
  \item{lwr, upr}{
For 1D plots: number vectors defining the lower and upper bounds to plot for a confidence interval. Must be the same length as \code{x} and \code{y}.
}
  \item{ci.lty, ci.lwd, ci.col}{
For 1D plots: the type, width, and color for the confidence interval lines drawn from the \code{lwr} and \code{upr} arguments.
}
  \item{zero, zero.lty}{
For 1D plots: \code{zero} is a logicical indicating whether a horizontal line at y = 0 should be included, and \code{zero.lty} controls the line type
}
  \item{xlim, ylim, zlim}{
For 1D plots: \code{xlim} and \code{ylim} are the axis limits input to the \code{\link{plot}} function. For 2D plots: these are the axis limits input to the \code{\link{image}} function (note: \code{zlim} controls range for color legend).
}
  \item{xlab, ylab, zlab}{
For 1D plots: \code{xlab} and \code{ylab} are the axis labels input to the \code{\link{plot}} function. For 2D plots: these are the axis labels input to the \code{\link{image}} function (note: \code{zlab} controls label for color legend).
}
  \item{main}{
Title of the plot.
}
  \item{add}{
Should lines/bars be added to current plot?
}
  \item{zline}{
For 2D plots: margin line for the z-axis label.
}
  \item{xticks, yticks}{
For 2D plots: tick marks for x-axis and y-axis grid lines.
}  
  \item{xlabels, ylabels}{
For 2D plots: labels corresponding to the input tick marks that define the grid lines.
}
  \item{...}{
Additional arguments passed to the \code{\link{plot}} and \code{\link{image}} functions.  
}
}
\details{
The \code{visualizer1} function is used to plot 1D (line/bar) plots, and the \code{visaulizer2} function is used to plot 2D (image) plots. These functions are not intended to be called by the user, but they may be useful for producing customized visualizations that are beyond the scope of the \code{\link{plot.gammi}} function.
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
The vector of colors used to construct the plots is defined as \code{colorRampPalette(col)(ncolor)}, which interpolates a color palette of length \code{ncolor} from the input colors in the vector \code{col}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{plot.gammi}} for plotting effects from \code{gammi} objects
}
\examples{
# load 'gammi' package
library(gammi)

# load 'exam' help file
?exam

# load data
data(exam)

# header of data
head(exam)

# fit model
mod <- gammi(Exam.score ~ VRQ.score, data = exam,
             random = ~ (1 | Primary.school) + (1 | Secondary.school))
             
# plot results (using S3 method)
plot(mod, include.random = FALSE)

# plot results (using visualizer)
xnew <- seq(min(exam$VRQ.score), max(exam$VRQ.score), length.out = 400)
pred <- predict(mod, newdata = data.frame(VRQ.score = xnew), 
                type = "terms", conf.int = TRUE)
visualizer1(x = xnew, y = pred$fit, lwr = pred$lwr, upr = pred$upr,
            xlab = "VRQ.score", ylab = "Exam.score", main = "VRQ.score effect")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{ hplot }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
