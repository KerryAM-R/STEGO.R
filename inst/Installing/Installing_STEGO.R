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
remotes::install_github("KerryAM-R/STEGO.R", ref = "beta-version", force = T)
