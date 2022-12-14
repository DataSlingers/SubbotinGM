% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/subbotin_GM.R
\name{subbotin_GM}
\alias{subbotin_GM}
\title{Edge selection with the Subbotin graphical model using neighborhood selection.}
\usage{
subbotin_GM(dat, lambda, nu, sym = "and", ...)
}
\arguments{
\item{dat}{n by p data matrix.}

\item{lambda}{value of regularization parameter.}

\item{nu}{power of Subbotin graphical model.}

\item{sym}{"and" vs. "or" rule for neighborhood selection.}

\item{...}{variable arguments for neighborhood regression function.}
}
\value{
A list with the following elements:
\itemize{
\item edges: p by p adjacency matrix.
\item param: p by p asymetric parameter matrix from neighborhood regressions
}
}
\description{
Edge selection with the Subbotin graphical model using neighborhood selection.
}
