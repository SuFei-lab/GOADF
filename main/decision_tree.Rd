% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/decision_tree.R
\name{decision_tree}
\alias{decision_tree}
\title{Build 100 trees}
\usage{
decision_tree(data, cp, prune)
}
\arguments{
\item{data}{A data frame in which the rows are samples, the columns are features, and the first column is the sample class.}

\item{cp}{A numeric variable with a complexity parameter. The larger the value of cp, the smaller the size of the split.}

\item{prune}{A logical variable that is choosen between 'T' and 'F' to determine whether to prune.}
}
\value{
a list of result and a plot of decision tree.
}
\description{
Build 100 trees
}
\examples{
decision_tree(data,cp=0.01,prune=T)
}
