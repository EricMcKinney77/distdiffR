# R package development workflow

################### Initial setup steps.

# (1) Initialize a GitHub Repo
# (2) Create a version controlled project for the package
# (3) Run the following code:

library(devtools)
packageVersion("devtools") # At least 2.3.2.9000

create_package("~/Documents/GitHub/distdiffR")

# (4) Make an initial commit, pull, and push to the remote master



################### Additions and Maintenance, Start here.

# Reload devtools
library(devtools)

# Depending on the additions / patches, bump the version number in the
# DESCRIPTION file (see https://yihui.org/en/2013/06/r-package-versioning/).
# Use x.y.z version numbers (major.minor.patch, e.g. 0.1.12).
# Run devtools::document() if file names change before rebuilding.
# May need to remove the old version of the package from your local machine and
# restart Rstudio before rebuilding.

# (5) Write functions, objects, etc...

fbind <- function(a, b) {
  factor(c(as.character(a), as.character(b)))
}

# (6) Create a .R file for the function
use_r("fbind")

# (7) Copy the function over to fbind.R

# (8) Click Build -> Install and Restart to make sure the package development is working well
# If you change the name of a .cpp file (and function) run devtools::document() first to update namespace

# (9) Commit, pull, push

# (10) Click Build -> Check to make sure all pieces of the package still work well with each other

# (11) Add documentation with Code > Insert Roxygen Skeleton

# (12) Use document() to publish the documentation to the man folder using Roxygen
devtools::document()

# Use build_manual() after updating the Rd document files to produce a pdf
# manual in the directory just above the package directory
devtools::build_manual()

# (13) Check and commit again



# (14) usethis::use_rcpp() for setting up the package for Rcpp

# Rcpp: Create a C++ script, write the code and include documentation. Then Ctrl + Shift + D.
# Install and Restart



# (15) usethis::use_vignette("package_name") to create a vignette folder and template



# (16) Use the Build menu (Build -> More -> Build source package) to create a .tar.gz file
# This is helpful for the CHPC install
# Otherwise, use:
library(devtools)
install_github("https://github.com/EricMcKinney77/distdiffR")
