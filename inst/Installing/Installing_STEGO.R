Sys.setenv("CXX" = "g++")

# installing STEGO.R
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!require("devtools", quietly = TRUE))
  install.packages("devtools")

if (!require("usethis", quietly = TRUE))
install.packages("usethis")

if (!require("startup", quietly = TRUE))
  install.packages("startup")

# if (!require("RcppArmadillo", quietly = TRUE))
#   install.packages("RcppArmadillo")

devtools::install_github("KerryAM-R/STEGO.R", ref = "beta-version")

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





installed.packages.full <- as.data.frame(installed.packages())
installed_packages <- rownames(installed.packages.full)

base_packages <- rownames(installed.packages.full[installed.packages.full$Priority %in% "base",])
base_packages

installed_packages[installed_packages %in% c(base_packages,"fonts","extrafont")]

non_base_packages <- installed_packages[!installed_packages %in% c(base_packages,"fonts","extrafont")]
non_base_packages %in% "fonts"

remove.packages(non_base_packages)
