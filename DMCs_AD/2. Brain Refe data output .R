## Get DNA methylation Brain reference

library(Biobase)
library(GEOquery)

# load series and platform data from GEO

gset <- getGEO("GSE41826", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL13534", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
rawdata <- exprs(gset)
status <- pData(gset)
newstatus <- status[,c(1,36,38,41,42)]
for(i in 1:118){
    newstatus$ID[i] = as.numeric(unlist(strsplit(as.character(status$title[i]),split="-"))[1])
    newstatus$CellType[i] = unlist(strsplit(as.character(status$title[i]),split="-"))[2]
}
for(i in c(66,68,119:nrow(status))){
    newstatus$ID[i] = as.numeric(unlist(strsplit(as.character(status$title[i]),split=" "))[2])
    newstatus$CellType[i] = unlist(strsplit(as.character(status$title[i]),split=" "))[1]
}
newstatus$CellType
expG = rawdata[,newstatus$`diagnosis:ch1`=="Control"&newstatus$`CellType`=="G"]
expN = rawdata[,newstatus$`diagnosis:ch1`=="Control"&newstatus$`CellType`=="N"]
dim(expG)

RefPanel <- cbind(apply(expG,1,mean),
                  apply(expN,1,mean))

dim(RefPanel)
colnames(RefPanel) = c("Glia","Neuron")

setwd("/Users/Ziyi/Dropbox/PartialRBDeconv")
save(RefPanel, file = "Brain_RefPanel_DNAMethy.rda")
