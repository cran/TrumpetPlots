% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/powerCurves.R
\name{powerCurves}
\alias{powerCurves}
\title{Power Curves for Trumpet Plots}
\usage{
powerCurves(threshold = 0.8, N = 4e+05, alpha = 5e-08, Nfreq = 500)
}
\arguments{
\item{threshold}{user-specified power level}

\item{N}{sample size}

\item{alpha}{significance threshold}

\item{Nfreq}{Number of allele frequency data points generated to calculate the power curves}
}
\value{
A data frame with the power estimated for each allele frequency and effect size, given a: Statistical power threshold, significance threshold (alpha value), and sample size
}
\description{
This function generates curves indicating statistical power in Trumpet plots
}
\examples{
powerCurves(threshold = 0.8, N=400000, alpha = 5e-8)
}
