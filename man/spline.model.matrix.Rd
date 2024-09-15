\name{spline.model.matrix}
\alias{spline.model.matrix}
\alias{rowKronecker}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Construct Design Matrices via Spectral Splines
}
\description{
Creates a design (or model) matrix using the \code{\link{spline.basis}} function to expand variables via a spectral spline basis. 
}
\usage{
spline.model.matrix(object, data, ...)

rowKronecker(X, Y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
a \code{\link{formula}} or \code{\link{terms}} object describing the fit model
}
  \item{data}{
a data frame containing the variables referenced in \code{object}
}
  \item{...}{
additional arguments passed to the \code{\link{spline.basis}} function, e.g., \code{df}, \code{knots}, \code{m}, etc. Arguments must be passed as a named list, see Examples.
}

  \item{X}{
matrix of dimension \eqn{n \times p}  
}
  \item{Y}{
matrix of dimension \eqn{n \times q}  
}
}
\details{
This is a reproduction of the \code{\link[grpnet]{rk.model.matrix}} function in the \bold{grpnet} package (Helwig, 2024b). 

Designed to be a more flexible alternative to the \code{\link{model.matrix}} function. The \code{\link{spline.basis}} function is used to construct a marginal basis for each variable that appears in the input \code{object}. Tensor product interactions are formed by taking a rowwise Kronecker product of marginal basis matrices. Interactions of any order are supported using standard formulaic conventions, see Note.
}
\value{
The design matrix corresponding to the input formula and data, which has the following attributes:

\item{assign}{an integer vector with an entry for each column in the matrix giving the term in the formula which gave rise to the column}
\item{term.labels}{a character vector containing the labels for each of the terms in the model}
\item{knots}{a named list giving the knots used for each variable in the formula}
\item{m}{a named list giving the penalty order used for each variable in the formula}
\item{periodic}{a named list giving the periodicity used for each variable in the formula}
\item{xlev}{a named list giving the factor levels used for each variable in the formula}
}
\references{
Helwig, N. E. (2021). Spectrally sparse nonparametric regression via elastic net regularized smoothers. \emph{Journal of Computational and Graphical Statistics, 30}(1), 182-191. \doi{10.1080/10618600.2020.1806855}

Helwig, N. E. (2024a). Precise tensor product smoothing via spectral splines. \emph{Stats, 7}(1), 34-53, \doi{10.3390/stats7010003}

Helwig, N. E. (2024b). grpnet: Group Elastic Net Regularized GLMs and GAMs. R package version 0.4. \doi{https://doi.org/10.32614/CRAN.package.grpnet}
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\note{
For formulas of the form \code{y ~ x + z}, the constructed model matrix has the form \code{cbind(spline.basis(x), spline.basis(z))}, which simply concatenates the two marginal basis matrices. For formulas of the form \code{y ~ x : z}, the constructed model matrix has the form \code{rowKronecker(spline.basis(x), spline.basis(z))}, where \code{\link{rowKronecker}} denotes the row-wise kronecker product. The formula \code{y ~ x * z} is a shorthand for \code{y ~ x + z + x : z}, which concatenates the two previous results. Unless it is suppressed (using \code{0+}), the first column of the basis will be a column of ones named \code{(Intercept)}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See \code{\link{spline.basis}} for details on the spectral spline basis
}
\examples{
# load data
data(exam)

# header of data
head(exam)

# make basis matrix 
x <- spline.model.matrix(Exam.score ~ ., data = exam)

# check dimension (= 3435 by 178)
dim(x)

# check term labels
attr(x, "term.labels")

# check which columns of x belong to which terms
attr(x, "assign")             # note: 0 = (Intercept)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{ regression }% use one of  RShowDoc("KEYWORDS")
\keyword{ smooth }% __ONLY ONE__ keyword per line
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
