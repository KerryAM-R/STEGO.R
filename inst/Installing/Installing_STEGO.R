Sys.setenv("CXX" = "g++")

# installing STEGO.R
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# if (!require("devtools", quietly = TRUE))
#   install.packages("devtools")

if (!require("usethis", quietly = TRUE))
  install.packages("usethis")

devtools::install_github("KerryAM-R/STEGO.R")

y # for Seurat source
y # biomanager packages for the M1 or greater chip

# import all of the fonts to your system. Will only need to be done the first time.
require(extrafont)
font_import()

# when it asks if you would like to install the different fonts
y


if (!require("startup", quietly = TRUE))
  install.packages("startup")

# installing STEGO
install.packages("usethis")
usethis::edit_r_environ()
# add to the .Renviron
# R_MAX_VSIZE=200Gb
.rs.restartR()

devtools::install_github("KerryAM-R/STEGO.R")
1
n
.rs.restartR()

# Run the application
# STEGO.R::Load_required_packages()
STEGO.R::runSTEGO()

# dir.create("/Users/kerrymullan/Library/R/arm64/4.4/library", recursive = TRUE)







# remove packages for testing purposes ------
installed.packages.full <- as.data.frame(installed.packages())

!"extrafont" %in% rownames(installed.packages.full)

installed_packages <- rownames(installed.packages.full)

base_packages <- rownames(installed.packages.full[installed.packages.full$Priority %in% "base",])
base_packages

installed_packages[installed_packages %in% c(base_packages,"fonts","extrafont","dplyr","plyr")]

non_base_packages <- installed_packages[!installed_packages %in% c(base_packages,"startup")]
# non_base_packages %in% "fonts"

remove.packages(non_base_packages)

# usethis::use_vignette("STEGO")
