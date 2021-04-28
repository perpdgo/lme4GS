## lme4GS: lme4 for Genomic Selection

Fits several models used routinely in Genomic Selection.

#### 1. Installation 

##### Option 1. Installing from source.

You need to install first some extra tools, [Xcode](https://developer.apple.com/xcode/downloads/) for Mac users, 
and [Rtools](https://cran.r-project.org/bin/windows/Rtools/) for Windows users. After that, open R and type the 
following commands to install the package:

```R

#Install lme4 package
install.packages("lme4")

#Install devtools package
install.packages("devtools")

#Load the devtools  package wich includes the install_git function
library(devtools)

#Install lme4GS from github
install_git('https://github.com/perpdgo/lme4GS/',subdir='pkg_src/lme4GS')

```

##### Option 2. Installing from a binary file.

Downlod the file "lme4GS_0.1.zip" in case you are using windows or "lme4GS_0.1.tgz" if you use macOS from 
the pkg_bin directory in this repository. Then, open Rstudio and execute the following commands:

```R

#Install lme4 package
install.packages("lme4")

```

After that:

1) Go to the Tools menu -> Install.packages.
2) Select "Package archive" in "the Install from" menu.
3) Select the file "lme4GS_0.1.zip" or "lme4GS_0.1.tgz", depending on your operating system.
4) Click in the "Install" button.

<img src="https://github.com/perpdgo/lme4GS/blob/master/screenshoots/install_macOS.png" width="500">

Figure 1. Install packages diaglog in macOS.

<img src="https://github.com/perpdgo/lme4GS/blob/master/screenshoots/install_win.png" width="500">

Figure 2. Install packages dialog in Windows.

#### 2. Fitting several models

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
 y<-as.vector(wheat.Y)

 data=data.frame(y=y,m_id=rep(rownames(G),4),a_id=rep(rownames(A),4))

 out<-lmerUvcov(y~(1|m_id)+(1|a_id),data=data,Uvcov=list(m_id=list(K=G),a_id=list(K=A)))

 summary(out)

 plot(y,predict(out))

```

**Prediction based on markers for mice**

```R

#Example 2, mice

data(mice)
X<-mice.X
A<-mice.A
Z<-scale(X,center=TRUE,scale=TRUE)
G<-tcrossprod(Z)/ncol(Z)

mice.pheno$a_id<-rownames(A)
mice.pheno$m_id<-rownames(G)

out<-lmerUvcov(Obesity.BMI~GENDER+Litter+(1|cage)+(1|m_id)+(1|a_id),data=mice.pheno,
               Uvcov=list(a_id=list(K=A),m_id=list(K=G)))
summary(out)

plot(predict(out),mice.pheno$Obesity.BMI)

```
#### 3. Cross-validation 

**Wheat examples**

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

  #Predict 10/100 of records selected at random. 
  #The data were partitioned in 10 groups at random
  #and we predict individuals in group 2.

  fold<-2
  y_trn<-y[wheat.sets!=fold]
  y_tst<-y[wheat.sets==fold]

  A_trn=A[wheat.sets!=fold,wheat.sets!=fold]
  G_trn=G[wheat.sets!=fold,wheat.sets!=fold]

  pheno_trn=data.frame(y_trn=y_trn,m_id=rownames(A_trn),a_id=rownames(G_trn))

  #######################################################################################
  #Marker based prediction
  #######################################################################################
	
  fm1<-lmerUvcov(y_trn~1+(1|m_id),data=pheno_trn,Uvcov=list(m_id=list(K=G_trn)))

  plot(pheno_trn$y_trn,predict(fm1),xlab="Phenotype",ylab="Pred. Gen. Value")

  #BLUP for individuals in the testing set
  blup_tst<-ranefUvcovNew(fm1,Uvcov=list(m_id=list(K=G)))
  blup_tst<-blup_tst$m_id[,1]

  #Comparison
  #Check the names
  names(y_tst)<-rownames(G)[wheat.sets==fold]
  blup_tst<-blup_tst[match(names(y_tst),names(blup_tst))]

  yHat_tst<-fixef(fm1)[1]+blup_tst
  points(y_tst,yHat_tst,col="red",pch=19)

  #Correlation in testing set
  cor(y_tst,yHat_tst)

  #######################################################################################
  #Pedigree based prediction
  #######################################################################################

  fm2<-lmerUvcov(y_trn~1+(1|a_id),data=pheno_trn,Uvcov=list(a_id=list(K=A_trn)))

  plot(pheno_trn$y_trn,predict(fm2),xlab="Phenotype",ylab="Pred. Gen. Value")

  #BLUP for individuals in the testing set
  blup_tst<-ranefUvcovNew(fm2,Uvcov=list(a_id=list(K=A)))
  blup_tst<-blup_tst$a_id[,1]

  #Comparison
  #Check the names
  names(y_tst)<-rownames(A)[wheat.sets==fold]
  blup_tst<-blup_tst[match(names(y_tst),names(blup_tst))]

  yHat_tst<-fixef(fm2)[1]+blup_tst
  points(y_tst,yHat_tst,col="red",pch=19)

  #Correlation in testing set
  cor(y_tst,yHat_tst)

  #######################################################################################
  #Markers + Pedigree based prediction
  #######################################################################################
  
  fm3<-lmerUvcov(y_trn~1+(1|m_id)+(1|a_id),data=pheno_trn,
               Uvcov=list(m_id=list(K=G_trn),a_id=list(K=A_trn)))

  plot(pheno_trn$y_trn,predict(fm3),xlab="Phenotype",ylab="Pred. Gen. Value")

  #BLUP for individuals in the testing set
  blup_tst<-ranefUvcovNew(fm3,Uvcov=list(m_id=list(K=G),a_id=list(K=A)))

  blup_tst_m<-blup_tst$m_id[,1]
  blup_tst_a<-blup_tst$a_id[,1]

  #Comparison
  #Check the names
  names(y_tst)<-rownames(A)[wheat.sets==fold]
  blup_tst_m<-blup_tst_m[match(names(y_tst),names(blup_tst_m))]
  blup_tst_a<-blup_tst_a[match(names(y_tst),names(blup_tst_a))]

  yHat_tst<-fixef(fm3)[1] + blup_tst_m + blup_tst_a
  points(y_tst,yHat_tst,col="red",pch=19)

  #Correlation in testing set
  cor(y_tst,yHat_tst)

```
