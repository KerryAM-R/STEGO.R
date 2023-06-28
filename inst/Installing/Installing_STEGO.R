# installing STEGO
install.packages("usethis")
usethis::edit_r_environ()
# add to the .Renviron
R_MAX_VSIZE=200Gb
GITHUB_PAT=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
.rs.restartR()

remotes::install_github("KerryAM-R/STEGO.R")
1
n
.rs.restartR()

# Run the application
STEGO.R::Load_required_packages()
STEGO.R::runSTEGO()
