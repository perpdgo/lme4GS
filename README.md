## lme4GS: lme4 for Genomic Selection

Fits several models used routinely in Genomic Selection.

#### 1. Installation 

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
ytrn<-y[wheat.sets!=fold]
ytst<-y[wheat.sets==fold]

#######################################################################################
#Marker based prediction
#######################################################################################

random<-list(mrk=list(K=G,id=names(ytrn)))
	
out<-lmer_uvcov(ytrn,fixed="1",random=random)

plot(ytrn,predict(out),xlab="Phenotype",ylab="Pred. Gen. Value")

#Random effect list for prediction
newrandom<-list(mrk=list(K=G,id=names(y)[wheat.sets==fold]))

blup_tst<-predict_uvcov(out,newrandom)
blup_tst<-blup_tst[[1]]
points(ytst,blup_tst,col="red",pch=19)

#Correlation in testing set
cor(ytst,blup_tst)

#######################################################################################
#Pedigree based prediction
#######################################################################################

random<-list(ped=list(K=A,id=names(ytrn)))
	
out<-lmer_uvcov(ytrn,fixed="1",random=random)

plot(ytrn,predict(out),xlab="Phenotype",ylab="Pred. Gen. Value")

#Random effect list for prediction
newrandom<-list(ped=list(K=A,id=names(y)[wheat.sets==fold]))

blup_tst<-predict_uvcov(out,newrandom)
blup_tst<-blup_tst[[1]]
points(ytst,blup_tst,col="red",pch=19)

#Correlation in testing set
cor(ytst,blup_tst)

#######################################################################################
#Markers + Pedigree based prediction
#######################################################################################
 
random<-list(mrk=list(K=G,id=names(ytrn)),
             ped=list(K=A,id=names(ytrn)))
	
out<-lmer_uvcov(ytrn,fixed="1",random=random)

plot(ytrn,predict(out),xlab="Phenotype",ylab="Pred. Gen. Value")

#Random effect list for prediction
newrandom<-list(mrk=list(K=G,id=names(y)[wheat.sets==fold]),
                ped=list(K=A,id=names(y)[wheat.sets==fold]))

blup_tst<-predict_uvcov(out,newrandom)
blup_tst<-blup_tst[[1]]+blup_tst[[2]]
points(ytst,blup_tst,col="red",pch=19)

#Correlation in testing set
cor(ytst,blup_tst)

```
