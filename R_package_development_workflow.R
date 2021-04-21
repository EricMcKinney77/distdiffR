# R package development workflow

# (1) Initialize a GitHub Repo
# (2) Create a version controlled project for the package
# (3) Run the following code:

library(devtools)
packageVersion("devtools") # At least 2.3.2.9000
library(tidyverse)
library(fs)

create_package("~/Documents/GitHub/foofactors")

# (4) Make an initial commit, pull, and push to the remote master

library(devtools)
library(tidyverse)
library(fs)


# (5) Write functions, objects, etc...

fbind <- function(a, b) {
  factor(c(as.character(a), as.character(b)))
}

# (6) Create a .R file for the function
use_r("fbind")

# (7) Copy the function over to fbind.R

# (8) Use load_all() to make sure the package development is working well
load_all()

# (9) Commit, pull, push

# (10) Use check() to make sure all pieces of the package still work well with each other
check()

# (11) Add documentation with Code > Insert Roxygen Skeleton

# (12) Use document() to publish the documentation to the man folder using Roxygen

# (13) check() and commit() again

# (14) usethis::use_rcpp() for setting up the package for Rcpp
