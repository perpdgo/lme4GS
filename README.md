## lme4GS: lme4 for Genomic Selection

Fits several models used routinely in Genomic Selection.

#### Installation 

**From GitHub (development version, unstable)**.

```R
#Install devtools package
install.packages(pkg='devtools',repos='https://cran.r-project.org')

#Install getPass package
install.packages(pkg='getPass',repos='https://cran.r-project.org')

#Load the devtools library
library(devtools)

#Load the getPass package
library(getPass)

#Install lme4GS from github
#Here the user can be "perpdgo" or "gldc"
install_git('https://github.com/perpdgo/lme4GS',credentials=git2r::cred_user_pass("perpdgo",getPass::getPass()))

```
