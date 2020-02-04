# MultiPhen
:exclamation: This is forked from the CRAN R package repository.  MultiPhen — A Package to Test for Pleiotropic Effects  

To run it you need to clone the repository and then run following in R (where path_to_file is the cloned directory)
The following was tested in R-3.6.1
```
path_to_file = "C:/Users/LCOIN/github/MultiPhen"
INSTALL = TRUE
if(INSTALL){
install.packages("abind")
install.packages("epitools")
install.packages("meta")
install.packages("HardyWeinberg")
install.packages("RColorBrewer")
install.packages("gplots")
install.packages(path_to_file, repos = NULL, type="source")
}
```
