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

#### Fitting several models

**Prediction based on markers and pedigree for wheat**

```R

 library(BGLR)
 library(lme4GS)

 #Example 1, wheat 
 data(wheat)
 X<-wheat.X
 Z<-scale(X,center=TRUE,scale=TRUE)
 G<-tcrossprod(Z)/ncol(Z)
 A<-wheat.A
 rownames(G)<-colnames(G)<-rownames(A)
 y<-wheat.Y[,1]

 #id a vector with ids
 random<-list(mrk=list(K=G,id=rownames(G)),
              ped=list(K=A,id=rownames(A)))

 out<-lmer_uvcov(y,fixed="1",random=random)

 summary(out)

 plot(y,predict(out))

```

**Prediction based on markers for mice**

```R

#Example 2, mice

  data(mice)

  y<-mice.pheno$Obesity.BMI
  y<-scale(y,center=TRUE,scale=TRUE)
  y<-as.vector(y)

  X<-mice.X
  Z<-scale(X,center=TRUE,scale=TRUE)
  G<-tcrossprod(Z)/ncol(Z)

  gender<-mice.pheno$GENDER
  litter<-mice.pheno$Litter
  cage<-mice.pheno$cage

  fixed<-"1+gender+litter"

  #K=NULL means that we are using the identity matrix
  random<-list(cage=list(K=NULL,id=cage),
	           mrk=list(K=G,id=rownames(G)))

  out<-lmer_uvcov(y,fixed=fixed,random=random)

  plot(y,predict(out))
```
