##=============================================
## overlap of DMCs between Joint and Minfi 
##=============================================
setwd("/Users/yanyuchen/Desktop/Thesis/data")

idx_joint <- read_csv("Thesis_DMC_Joint.csv") %>% 
  as.data.frame %>% 
  .[["CpGs"]]

idx_all <- read_csv("Thesis_DMC_All.csv") %>% 
  as.data.frame %>% 
  .[["CpGs"]]

load("AD_data.rda")
ol.names = rownames(AD_data)
length(ol.names)

# overlap is 1232 with Minfi
sum(idx_joint %in% idx_all)

##=============================================
## Match de jager 71 cpgs with Joint data
##=============================================

## De jager 71 cpgs
idx_dejager<-c("cg11724984","cg23968456","cg15821544","cg16733298","cg22962123","cg13076843",
               "cg25594100","cg00621289","cg19803550","cg03169557","cg05066959","cg05810363",
               "cg07012687","cg21207436", "cg21806242","cg11823178","cg12163800","cg05417607",
               "cg17474422","cg13390284","cg22904711","cg27041424","cg18556455","cg05731218",
               "cg17693222","cg12307200","cg19007269","cg15645660","cg22883290","cg14074251",
               "cg04252044","cg09448088","cg02308560","cg07180538","cg20733077","cg13639901",
               "cg20618448","cg03193328","cg07883124","cg24231804","cg12877335","cg22941668",
               "cg20546777","cg14430943","cg24676346","cg19140834","cg10920329","cg11652496",
               "cg08737189","cg25917732","cg02342148","cg16140558","cg26407544","cg13579486",
               "cg07714812","cg18659586","cg15348679","cg16459281","cg04157161","cg04126866",
               "cg18346707","cg22385702","cg27443779","cg18343862","cg06742628","cg12114584",
               "cg05322931","cg27409802","cg06753513","cg08324801","cg21644387")
idx_ad<-c("cg22883290","cg02308560")

## Identify overlap between CpGs of DMCs of Joint and dejager's
match_dejager_joint<-match(idx_dejager,idx_joint)

## Identify overlap between CpGs of DMCs of Minfi and dejager's
match_dejager_all<-match(idx_dejager,idx_all)

## Identify overlap between CpGs of our raw data and dejager's
match_allcpgs<-match(idx_dejager,ol.names)

## Output CpGs of the results of three overlap above
output_dejager<-data.frame(idx_dejager,match_dejager_joint,match_dejager_all,match_allcpgs)
colnames(output_dejager)=c("CpGs","Ranking number in DMCs of Joint","Ranking number in DMCs of Minfi","Ranking number in raw data")
write_csv(output_dejager, "Thesis_match_dejager_Joint.csv")


##=============================================
## output identified csDM information of Joint
##=============================================
setwd("/Users/yanyuchen/Desktop/Thesis/data")
idx_joint <- read_csv("Thesis_DMC_Joint.csv") %>% 
  as.data.frame %>% 
  .[["CpGs"]]

library(dplyr)
library(readr)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
data("IlluminaHumanMethylation450kanno.ilmn12.hg19")

### Macth Diffrentially methylated CpGs of Joint with its information of genes

anno_j       <- IlluminaHumanMethylation450kanno.ilmn12.hg19 %>% 
  getAnnotation %>% 
  as.data.frame %>% 
  dplyr::slice(match(idx_joint, Name))

## Writing out information of Joint to csv file
write_csv(anno_j, "Thesis_DMC_Joint_geneinfo.csv")

## output profile to conduct pathway analysis
anno_j_info<-cbind(anno_j[,1],anno_j[,2],anno_j[,2]+1)
colnames(anno_j_info)=c("chr","start","end")
write.table(anno_j_info, "Thesis_DMC_Joint_pathway.txt",sep="\t",row.names = FALSE,quote = FALSE)
