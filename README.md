# MoreyDavisStober_pcurveASA

This repository contains the data, code and text for Morey &amp; Davis-Stober's "On the poor statistical properties of the P-curve meta-analytic procedure" (2025, Journal of the American Statistical Association).

Morey, R. D., & Davis-Stober, C. P. (2025). On the poor statistical properties of the P-curve meta-analytic procedure. *Journal of the American Statistical Association, 1–19*. https://doi.org/10.1080/01621459.2025.2544397

See also the **live P-curve app** at https://richarddmorey.github.io/pcurveAppTest (whose code can be found at https://github.com/richarddmorey/pcurveAppTest).

The compiled paper, with supplement, can be downloaded from `/text/asa_article/` [[link directly to PDF](https://raw.githubusercontent.com/richarddmorey/MoreyDavisStober_pcurveASA/refs/heads/main/text/asa_article/Morey_Davis-Stober_2025_JASA_with_supplement.pdf)].

## Contents of important folders

| Folder                | Contents                   |
|:--------|:--------------------------------------------|
| `/data/`                | Data for meta-analysis of papers citing Simonsohn et al. (2015) |
| `/R_utility/`           | Primary utility functions for computing P-curve tests, power, etc |
| `/R_utility/Simonsohn/` | Archive of official P-curve app code, as used by some utility functions |
| `/text/asa_article/`    | The template files for the paper (main and supplement) |
| `/text/bib/`            | The bibtex file for references |
| `/text/sections/`       | The text of the paper, organized by section |

## Reproducing the manuscript and supplement

We use the R package `renv` for reproducibility.

In order to reproduce the manuscript, follow these steps.

1. Install R and Rstudio
2. Install the `renv` package within R
3. Download or check out the repository (and decompress if necessary)
4. Open the `pcurve` Rstudio project file within R studio
5. Activate `renv` in the project with the command `renv::activate()`
6. Install all required packages with the command `renv::restore()`
7. Compile the manuscript file `text/asa_article/asa_paper.Rmd` to PDF
8. Compile the supplement file `text/asa_article/supplement.Rmd` to PDF

