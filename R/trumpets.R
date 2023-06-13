#' Trumpets
#'
#' This function generates trumpet plots
#'
#' @param dataset Input text file with genetic association results. Columns required are rsID, freq, A1_beta, Analysis and Gene.
#' @param rsID (required) Single Nucleotide Polymorphism (SNP) name.
#' @param freq (required) allele frequency of effect SNP.
#' @param A1_beta (required) risk allele effect size.
#' @param Analysis (optional) adds colour to the type of analysis (e.g. GWAS, Sequencing).
#' @param Gene (optional) Candidate gene name (can be empty).
#' @param calculate_power (TRUE/FALSE) Calculate power curves. Choose TRUE to add power curves for a given threshold, alpha, sample size N and number of allele frequencies. Choose FALSE if you already ran powerCurves() outside or do not want to show power curves.
#' @param show_power_curves (TRUE/FALSE) Show power curves in plot
#' @param exist_datapwr Existing dataframe containing columns: freq, pos.b.for.f, neg.b.for.f, powerline.
#' @param threshold Required if power == TRUE. Can be a single number or a vector of statistical power thresholds.
#' @param N (Required if calculate_power == TRUE). Sample size used to test the association.
#' @param alpha (Required if calculate_power == TRUE).
#' @param Nfreq (Required if calculate_power == TRUE). Number of allele frequency data points generated to calculate the power curves. We recommend Nfreq>1000 for power curves with high resolution. Note that this will slow down the rendering of the plot.
#' @param power_color_palette A vector of colours for the power curves. Number of colors should match number of thresholds supplied.
#' @param analysis_color_palette A vector of colours for the analysis types.
#'
#' @return Creates a Trumpet plot with variant allele frequency (X axis, log10 scale) and effect size information (Y axis).
#' @export
#'
#' @importFrom purrr list_rbind
#' @import ggplot2
#' @import data.table
#' @import magrittr
#' @import stats
#'
#' @examples plot_trumpets(dataset = toy_data)
plot_trumpets = function(dataset = toy_data, rsID = "rsID", freq = "freq", A1_beta = "A1_beta", Analysis = "Analysis", Gene = "Gene",
                         calculate_power = TRUE, show_power_curves = TRUE, exist_datapwr = NULL,
                         threshold = c(0.7, 0.9), N = 100000, alpha = 5e-8, Nfreq=500,
                         power_color_palette = c("purple", "deeppink"),
                         analysis_color_palette = c("#018571", "#a6611a")
                         ){

  # Rename columns
  data.table::setnames(dataset, rsID, "rsID")
  data.table::setnames(dataset, freq, "freq")
  data.table::setnames(dataset, A1_beta, "A1_beta")
  data.table::setnames(dataset, Analysis, "Analysis")
  data.table::setnames(dataset, Gene, "Gene")

  # Flip allele frequency and effect size when freq>0.5 so we only report MAF
  dataset[, A1_beta := ifelse(freq > 0.5, -A1_beta, A1_beta)]
  dataset[, freq := ifelse(freq > 0.5, 1-freq, freq)]

  # Add placeholder variant to keep minimum value in Plot X axis = 1e-5
  placeholder = data.table::data.table(rsID = "", freq=1e-5, A1_beta=0, Analysis="", Gene="")
  dataset = rbind(dataset[, c("rsID", "freq", "A1_beta", "Analysis", "Gene")], placeholder)

  # Call powerCurves function to calculate power
  if(calculate_power & show_power_curves){
    if(length(threshold) > 1) {
      #if user supplies more than 1 threshold
      message(paste0(" Drawing ", threshold, " statistical power curve. \n Association analysis with ", as.integer(N), " individuals. \n ", Nfreq, " allele frequency data points simulated. \n"), sep=" ")
      datapwr = map(threshold, ~ powerCurves(threshold = .x, N = N, alpha = alpha, Nfreq = Nfreq)) %>%
        list_rbind()
    } else {
      message(paste0(" Drawing ", threshold, " statistical power curve. \n Association analysis with ", as.integer(N), " individuals. \n ", Nfreq, " allele frequency data points simulated. \n"), sep=" ")
      datapwr = powerCurves(threshold = threshold, N = N, alpha = alpha, Nfreq = Nfreq)
    }
  } else if(!calculate_power & !show_power_curves){
    # don't calculate power
    threshold = 0
    message(" Drawing Trumpet plots without Statistical power curves. \n")
    datapwr = data.table::data.table(freq=1, pos.b.for.f = 0, neg.b.for.f = 0, powerline = 0)
  } else if(calculate_power & !show_power_curves){
    warning(" Did you mean to choose calculate_power = T and show_power_curves = F?. \n  (To see power curves you need both calculate_power = T and show_power_curves = T)")
    threshold = 0
    message(" Drawing Trumpet plots without Statistical power curves. \n")
    datapwr = data.table::data.table(freq=1, pos.b.for.f = 0, neg.b.for.f = 0, powerline = 0)
  } else if(!calculate_power & show_power_curves){
    warning(" Did you mean to choose calculate_power = F and show_power_curves = T?. \n  (To see power curves you need both calculate_power = T and show_power_curves = T)")
  } else if(!calculate_power & show_power_curves) {
    datapwr = exist_datapwr
  }

  # Plot

  pos_lines <- imap(unique(datapwr$powerline), ~ geom_line(data = datapwr[powerline == .x],
                                                           aes(x = freq, y = pos.b.for.f), lwd = 0.15,
                                                           col = power_color_palette[.y]))
  neg_lines <- imap(unique(datapwr$powerline), ~ geom_line(data = datapwr[powerline == .x],
                                                           aes(x = freq, y = neg.b.for.f), lwd = 0.15,
                                                           col = power_color_palette[.y]))

  p = dataset %>%
    .[!Analysis == ""] %>%
    ggplot2::ggplot() +

    pos_lines +
    neg_lines +

    geom_hline(yintercept = 0, col="lightgrey", lwd=0.5, lty=3) +
    # fill and colour provides information about the Data.source
    suppressWarnings(geom_point(aes(x=`freq`, y=A1_beta, fill=Analysis, color=Analysis, size=abs(A1_beta), alpha=0.8,
                   text=paste0(' Marker: ', rsID ,'\n A1_beta: ', A1_beta, '\n Freq: ',freq, '\n Gene: ',Gene, sep=" ")))) +
    scale_size(range = c(0.05, 5)) +
    scale_x_continuous(
      trans = "log10",
      breaks = c(1e-5, 0.0001, 0.001, 0.01, 0.05, 0.1, 0.2, 0.5, 1),
      labels = c(1e-5, 0.0001, 0.001, 0.01, 0.05, 0.1, 0.2, 0.5, 1)) +
    ylab("Effect size") +
    xlab("Effect allele frequency") +
    theme_classic() +
    scale_colour_manual(values = analysis_color_palette) +
    scale_fill_manual(values = analysis_color_palette) +
    theme(
      # strip are the labels (1 and 2) that we gave to split the plot in two.
      # we dont want them to appear, so we call them element_blank()
      strip.background = element_blank(),
      strip.text = element_blank(),
      # Edit size of the axis and title
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 16),
      # Adds a bit extra space for the X axis (the two panels created)
      panel.spacing = unit(1.5, "lines"),
      # Reduce size of legend text and title
      legend.title = element_text(size = 10, face="bold"),
      legend.text = element_text(size = 10),
      legend.position = "none")

  p
}
utils::globalVariables(c("toy_data","imap","map","powerline","trait","."))