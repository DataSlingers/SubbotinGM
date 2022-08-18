# SubbotinGM

This package implements the estimation algorithm for the Subbotin graphical model, as described in Chang and Allen (2022+). 

# Installation
The package can be installed by running the following commands in R:
```
  # install.packages(devtools)
  library(devtools)
  devtools:install_github("DataSlingers/SubbotinGM")
```

# Usage
This package comes with two functions:
 - `subbotin_GM` serves as the main tool for graph estimation. 
 - `extreme_lasso` performs the individual neighborhood selection regressions.
