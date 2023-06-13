
#' Power Curves for Trumpet Plots
#'
#' This function generates curves indicating statistical power in Trumpet plots
#'
#' @param threshold user-specified power level
#' @param N sample size
#' @param alpha significance threshold
#' @param Nfreq Number of allele frequency data points generated to calculate the power curves
#'
#' @return A data frame with the power estimated for each allele frequency and effect size, given a: Statistical power threshold, significance threshold (alpha value), and sample size
#' @export
#'
#' @import data.table
#' @import stats 
#' @import magrittr
#'
#' @examples powerCurves(threshold = 0.8, N=400000, alpha = 5e-8)
powerCurves = function(threshold = 0.8, N=400000, alpha = 5e-8, Nfreq=500){

  # Function to generate frequencies with log scale
  lseq <- function(from=1, to=100000, length.out=6) {
    # logarithmic spaced sequence
    # blatantly stolen from library("emdbook"), because need only this
    exp(seq(log(from), log(to), length.out = length.out))
  }

  # --- Add curve indicating the effect sizes required to achieve power at a given AF power wanted
  pw.thresh = threshold
  # significance threshold (aka alpha)
  p.threshold = alpha
  # calculate the chi-square value corresponding to significance threshold defined in p.threshold
  q = qchisq(p.threshold, df = 1, lower.tail = F)

  # Sequence of frequencies from min to max AF (MAF=0.5)
  f = lseq(from=1e-5, to=0.5, length.out=Nfreq)
  # Sequence of effect sizes from min beta to max beta
  b = seq(0, 10, length = Nfreq)
  # create null variable b.for.f
  b.for.f = rep(NA,length(f))

  # Calculate the standard deviation of the error term after removing SNP effect
  # In some cases there is warning when NAs are produced. Using supressWarnings to silence this
  sigma = suppressWarnings(sqrt(1 - 2*f*(1-f)*b^2))

  for(i in 1:length(f)){
    # Calculate power at this allele frequency, across a range of effect sizes (b)
    pwr_con = pchisq(q, df = 1, ncp=(N*2*f[i]*(1-f[i])*b^2)/sigma^2, lower.tail=F)
    # Calculate what is the minimum b needed to read pw.thres
    b.for.f[i] = b[min( which(pwr_con > pw.thresh))]
  }
  power_dt_maf = data.table::data.table(freq=f, pos.b.for.f = b.for.f, neg.b.for.f = -b.for.f) %>%
    .[, trait := "test"] %>%
    .[, powerline := threshold]

  return(power_dt_maf)
}

utils::globalVariables(c("toy_data","imap","map","powerline","trait","."))