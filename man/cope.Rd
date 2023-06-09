% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cope.R
\name{cope}
\alias{cope}
\title{Covariance parameter estimation for a geostatistical linear model}
\usage{
cope(formula, data, x_column, y_column, spatial_model = "exponential",
  random_formula = NULL, use.anistropy = FALSE, thetaini = c(2, 1,
  0.02), estMeth = "REML", par = FALSE, subSampCol = NULL,
  optMeth = "Nelder-Mead", profile = TRUE)
}
\arguments{
\item{formula}{an R linear model formula}

\item{data}{an data object with spatial coordinates.  If a plain data.frame, then xcoordscol and ycoordscol must be specified.  If an sp object, xcoordscol and ycoordscol should be NULL (the default).}

\item{x_column}{name, in quotes, of the column containing the x-coordinate (if it was not possible to obtain the coordinates from the data class). Default is NULL.}

\item{y_column}{name, in quotes, of the column containing the x-coordinate (if it was not possible to obtain the coordinates from the data class). Default is NULL.}

\item{spatial_model}{spatial autocorrelation models 
for random errors.  The list of spatial autocorrelation 
models is "exponential","expRadon2","expRadon4","gaussian","stable",
"rationalQuad","cauchyGrav","cauchyMag","cauchy","circular","spherical",
"cubic","penta","cardinalSine","besselK","besselJ"  Default is "exponential".}

\item{random_formula}{a list of variance components, 
Any names in the list not given above will be searched among the columns in the data set and
used as a factor variable for levels of a traditional random effect.}

\item{estMeth}{estimation method.  Default is "REML" for restricted maximum likelihood.  Other options are "ML" for maximum likelihood}

\item{useAnisotropy}{include anistropy in parameter estimation?  Default is "FALSE"}

\item{subsample_col}{A column of factors indicating grouping for use in subsampling for large data sets.  Default is "NULL," which creates a single grouping.}
}
\value{
a list of class "slm_cope".  The functions "summary" and "print" are used to obtain and print a summary. "anova" returns just the analysis of variance table...
}
\description{
Covariance parameter estimation for a geostatistical linear model
}
\author{
Jay Ver Hoef
}
