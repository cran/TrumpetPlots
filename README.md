# TrumpetPlots Rpackage

## About and preprint 
TrumpetPlots is an R package to visualize the effect size of risk variants across the allele frequency spectrum.

For more information about Trumpet plots and to cite our work, visit our preprint [Trumpet plots: Visualizing The Relationship Between Allele Frequency And Effect Size In Genetic Association Studies](https://www.medrxiv.org/content/10.1101/2023.04.21.23288923v1)

## Installation

To install, please install the "devtools" R package, then run

`devtools::install_git("https://gitlab.com/JuditGG/trumpetplots.git")`

## Usage

The function trumpets only requires a text file as input, containing association results. Please indicate the column names of the input dataset using the following arguments:

- `dataset`: Input text file with genetic association results. Columns required are rsID, freq, A1_beta, Analysis and Gene.
- `rsID`: (required) Single Nucleotide Polymorphism (SNP) name.
- `freq`: (required) allele frequency of effect SNP.
- `A1_beta`: (required) risk allele effect size.
- `Analysis`: (optional) adds colour to the type of analysis (e.g. GWAS, Sequencing).
- `Gene`: (optional) Candidate gene name (can be empty).
- `calculate_power`: (TRUE/FALSE) Calculate power curves. Choose TRUE to add power curves for a given threshold, alpha, sample size N and number of allele frequencies. Choose FALSE if you already ran powerCurves() outside or do not want to show power curves.
- `show_power_curves`: (TRUE/FALSE) Show power curves in plot. Needs argument `calculate_power = TRUE`.
- `threshold`: (Required if calculate_power == TRUE). Can be a single number or a vector of statistical power thresholds.
- `N`: (Required if calculate_power == TRUE). Sample size used to test the association.
- `alpha`: (Required if calculate_power == TRUE).  
- `Nfreq`: (Required if calculate_power == TRUE). Number of allele frequency data points generated to calculate the power curves. We recommend Nfreq>1000 for power curves with high resolution. Note that this will slow down the rendering of the plot.
- `power_color_palette`: A vector of colours for the power curves. Number of colors should match number of thresholds supplied.
- `analysis_color_palette`: A vector of colours for the argument `Analysis`. 
- `trait`: (optional) Name of the trait.

## Use 'toy_data' to run an example

A 'toy_data' dataset is provided to test the function. This toy_data contains 8000 genetic associations and seven columns to showcase function 'trumpets'. Data was obtained from GWAS and exome sequencing association analyses performend in the UK Biobank.

1. Load the toy dataset: `data(toy_data)`
2. Run the trumpets function: `plot_trumpets(dataset = toy_data)`


## R Shiny application 

For users with no knowledge of R, and to visualize Trumpet plots for >100 traits in the UK Biobank, visit the R Shiny application "Shiny Trumpets", available at: [https://juditgg.shinyapps.io/shinytrumpets/](https://juditgg.shinyapps.io/shinytrumpets/) 

## Report issues, ask questions, contribute and improve our code!
- For questions/issues related to the R package, visit [https://gitlab.com/JuditGG/trumpetplots/-/issues](https://gitlab.com/JuditGG/trumpetplots/-/issues)

- For questions/issues related to the R Shiny app, visit [https://gitlab.com/JuditGG/freq_or_plots/-/issues](https://gitlab.com/JuditGG/freq_or_plots/-/issues)