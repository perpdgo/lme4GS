#Classes for lmerUvcov

setClass("lmerUvcov", representation = list(resp="lmerResp",relfac = "list"),
          contains = "merMod")          
