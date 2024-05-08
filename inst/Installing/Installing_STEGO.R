# installing STEGO.R
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!require("devtools", quietly = TRUE))
  install.packages("devtools")

if (!require("usethis", quietly = TRUE))
install.packages("usethis")

devtools::install_github("KerryAM-R/STEGO.R")

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
remotes::install_github("KerryAM-R/STEGO.R", ref = "beta-version")




.libPaths(.libPaths()[!grepl("/Users/kerrymullan/Library/R/arm64/4.4/library", .libPaths())])
.libPaths("/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library")
.libPaths("/Users/kerrymullan/Library/R/arm64/4.4/library")

installed.packages.full <- as.data.frame(installed.packages())
installed_packages <- rownames(installed.packages.full)

base_packages <- rownames(installed.packages.full[installed.packages.full$Priority %in% "base",])
base_packages


non_base_packages <- installed_packages[!installed_packages %in% base_packages]
non_base_packages

remove.packages(non_base_packages)
