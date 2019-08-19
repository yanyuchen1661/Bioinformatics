
##===========================================
## load data
##===========================================
setwd("/Users/yanyuchen/Desktop/Thesis/data")

## load brain reference (include glia and neuron profiles)
load("Brain_RefPanel_DNAMethy.rda")

## load real dataset
load("probe_sample_qcd_betavalues_nobmiq_nocombat.rda")

## get overlapped CpGs between real data and reference panel
ol.names = intersect(rownames(RefPanel),rownames(betas))

## 339162 overlap CpGs
length(ol.names)

Y = betas[ol.names,]
Reference = RefPanel[ol.names,]

###=========================================================
### find out which subjects are AD, which are controls
###=========================================================

## read in phenotype of data
library(readr)
id_covar <- read_tsv('ROSMAP_arrayMethylation_covariates.tsv')
id_covar$idx<-paste(id_covar$Sentrix_ID,id_covar$Sentrix_Position,sep="_")
id_phenotype<-read_csv('Savannah_ROSMAPphenotypes.csv')
id_merge<-merge(id_covar,id_phenotype,by="Sample")

Control_colname=AD_colname=vector()

for (i in 1:739)
{  
  if (id_merge$braaksc[i]==0 || id_merge$braaksc[i]==1 || id_merge$braaksc[i]==2 || id_merge$braaksc[i]==3)
    
  {Control_colname[i]<-id_merge$idx[i]}
  
  else 
    
  {AD_colname[i]<-id_merge$idx[i]}
}


## Delete NAs of AD and control 
AD_colname <- AD_colname[!is.na(AD_colname)]
Control_colname<-Control_colname[!is.na(Control_colname)]

## get the subset for AD and control
AD_data<-subset(Y, select=intersect(colnames(betas),AD_colname))
Control_data<-subset(Y, select=intersect(colnames(betas),Control_colname))

## 366 AD 368 control
dim(AD_data)
dim(Control_data)

## merge 734 subjects' phenotype and covariates information into one dataset
subject_input<-intersect(colnames(betas),c(AD_colname,Control_colname))
id_final<-id_merge[id_merge$idx %in% subject_input,]

save(AD_data, file = "AD_data.rda")
save(Control_data, file = "Control_data.rda")
save(Reference,file = "Reference.rda")