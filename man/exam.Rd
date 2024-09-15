\name{exam}
\alias{exam}
\docType{data}
\title{
Cross-Classified Examination Data
}
\description{
Scores on secondary school leaving examinations (response) and verbal reasoning scores in primary school (fixed effect) for 3435 students in Fife, Scotland. The students are cross-classified in 148 primary schools (random effect) and 19 secondary schools (random effect).
}
\usage{data("exam")}
\format{
  A data frame with 3435 observations on the following 4 variables.
  \describe{
    \item{\code{VRQ.score}}{Verbal Reasoning Quotient obtained in primary school (integer vector ranging from 70 to 140)}
    \item{\code{Exam.score}}{Leaving examination score obtained in secondary school (integer vector ranging from 1 to 10)}
    \item{\code{Primary.school}}{Primary school identifier (factor with 148 levels)}
    \item{\code{Secondary.school}}{Secondary school identifier (factor with 19 levels)}
  }
}
\details{
The VRQ scores were obtained at age 12 (right before entering secondary school), and the Exam scores were obtained at age 16 (right before leaving secondary school). The VRQ scores are constructed to have a population mean of 100 and population standard deviation of 15. The goal is to predict the leaving Exam scores from the VRQ scores while accounting for the primary and secondary school cross-classifications.
}
\source{
Data Obtainable from:
\url{https://www.bristol.ac.uk/cmm/team/hg/msm-3rd-ed/datasets.html}
}
\references{
Goldstein, H. (2011). Multilevel Statistical Models, 4th Edition. Chapter 12: Cross-classified data structures (pages 243-254). \doi{10.1002/9780470973394} 

Paterson, L. (1991). Socio-economic status and educational attainment: a multidimensional and multilevel study. Evaluation and Research in Education, 5, 97-121. \doi{10.1080/09500799109533303}
}
\examples{
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

# variance parameters
mod$VarCorr
}
\keyword{datasets}
