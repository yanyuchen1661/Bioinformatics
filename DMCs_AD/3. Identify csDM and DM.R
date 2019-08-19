##===========================================
## identify csDM for Glia and Neuron
##===========================================



##===========================================
## solve for proportions
##===========================================

## set your working directory
setwd("/Users/yanyuchen/Desktop/Thesis/data")
load("AD_data.rda")
load("Control_data.rda")
load("Reference.rda")

source("sim_functions.R")
source("makeDesign.R")
source("fitModel.R")
source("csTest.R")
source("DEKtissue.R")

## solve proportions using EpiDISH
library(EpiDISH)
## choose top 10000 tissue specific genes
refinx.tissue = findRefinx.cv(Reference, nmarker=1000)
epidishout_AD = epidish(AD_data[refinx.tissue,], Reference[refinx.tissue,], method = "RPC")
Prop_AD = epidishout_AD$estF  ## solve proportions for AD's

epidishout_control = epidish(Control_data[refinx.tissue,], Reference[refinx.tissue,], method = "RPC")
Prop_control = epidishout_control$estF  ## solve proportions for controls

colnames(Prop_AD) =  colnames(Prop_control) = colnames(Reference)


orig_Y = cbind(AD_data, Control_data)
Y = asin(2*orig_Y -1)
dim(Y)
Prop = rbind(Prop_AD[,1:2], Prop_control[,1:2])
design <- data.frame(AD=factor(c(rep(1,ncol(AD_data)), rep(0,ncol(Control_data)))))
colnames(Prop) <- c("Glia", "Neuron")
# generate design matrix
Design_out <- makeDesign(design, Prop)
# fit model
fitted_model <- fitModel(Design_out, orig_Y)
# do test
#### To test cell type specific differences: you should get exactly the same p values as previous results

test <- csTest(fitted_model, coef = "AD", cell_type = "Glia", contrast_matrix = NULL,sort = TRUE)
head(test$res_table)
sum(test$res_table$fdr<0.05)

test_2 <- csTest(fitted_model, coef = "AD", cell_type = "Neuron", contrast_matrix = NULL)
head(test_2$res_table)
sum(test_2$res_table$fdr<0.05)

#### jointly test AD effect in neuron and glia
test_3 <- csTest(fitted_model, coef = "AD", cell_type = "joint", contrast_matrix = NULL)
head(test_3$res_table)
sum(test_3$res_table$fdr<0.05)

Joint<-cbind(rownames(test_3$res_table[test_3$res_table$fdr<0.05,]),test_3$res_table[test_3$res_table$fdr<0.05,])
colnames(Joint)[1]=c("CpGs")
library(readr)
write_csv(Joint, "Thesis_DMC_Joint.csv")

library(minfi)
design = c(rep(1,ncol(AD_data)), rep(0,ncol(Control_data)))
result=dmpFinder(as.matrix(Y), pheno=design, type="categorical")
head(result[result$qval<0.05,])
sum(result$qval<0.05, na.rm = T)

All<-cbind(rownames(result[result$qval<0.05,]),result[result$qval<0.05,])
colnames(All)[1]=c("CpGs")

write_csv(All, "Thesis_DMC_All.csv")
