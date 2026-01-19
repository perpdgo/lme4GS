#Classes for lmerUvcov and glmerUvcov

setClass("lmerUvcov", representation = list(resp="lmerResp",relfac = "list"),
          contains = "merMod")
setClass("glmerUvcov", representation = list(resp = "glmResp", relfac = "list"),
         contains = "merMod")

