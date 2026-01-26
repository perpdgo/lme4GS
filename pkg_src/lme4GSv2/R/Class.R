#Classes for lmerUvcov and glmerUvcov

#1. Create common superclass

setClass("merUvcov",
         contains = "merMod",
         representation = list(
           relfac = "list"
         ))

#2. lmerUvcov and glmerUvcov inherits from it

setClass("lmerUvcov",
         contains = "merUvcov",
         representation = list(resp = "lmerResp"))

#glmerResp is for internal use only in lme4, is not exposed...
#is an unregistred class
setClass("glmerUvcov",
         contains = "merUvcov",
         representation = list(resp = "ANY"))
