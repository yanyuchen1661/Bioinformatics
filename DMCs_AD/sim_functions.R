library(matrixStats)
library(ROCR)
library(gtools)
library(quadprog)
library(CellMix)
library(sirt)

# Ntissue=4
# pure_based_mean = logpure


# Ntissue=5
# pure_based_mean = blood.mean[,1:5]

### generate pure tissue profiles based on input paramters
# larger t_noise means larger noise on pure profiles
getOnePureTissue = function(profmean=pure_based_mean, lfc, normal=T, t_noise=1){

  N_feature = dim(profmean)[1] # number of features
  L = dim(profmean)[2] # number of pure tissues
  profvar = logpure.sd * t_noise

  if(normal){
        tissue = matrix(0,N_feature,L)
        for(i in 1:L){
            tissue[,i] = exp(rnorm(N_feature, profmean[,i], abs(profvar[,i])))
        }
  }else{
        dis.profmean = tissue = matrix(0, N_feature, L)
        for(i in 1:L){
            dis.profmean[,i] = lfc[,i] + profmean[,i]
            tissue[,i] = exp(rnorm(N_feature, dis.profmean[,i], abs(profvar[,i])))
        }
  }

  return(tissue)
}

### generate proportions based on input parameters
getProportion_old = function(N_sample, L=5, existprop=out2$Omega, prop_n_param=1000){
     prop.matrix = matrix(0, N_sample, L)
     prop.matrix.true = abs(existprop[sample(1:dim(existprop)[1],N_sample,replace=T),]) # sample proportions from ASD study

     # add noise on proportion matrix
     # larger prop_n_param means lower noise
     prop.matrix = t(apply(prop.matrix.true,1,function(x) rdirichlet(1,x*prop_n_param)))
     return(list(prop.matrix.true=prop.matrix.true, prop.matrix=prop.matrix))
}

### generate proportions based on input parameters
getProportion = function(N_sample, L=Ntissue, prop_n_param=1000){

  # get MLE for dirichlet distribution
  tmp=dirichlet.mle(out2$Omega[1:20,1:4]) # AD proportions
  alpha.case=tmp$alpha
  tmp=dirichlet.mle(out2$Omega[35:47,1:4])  # control proportions
  alpha.ctr=tmp$alpha

  prop.matrix.ctr = rdirichlet(N_sample, alpha.ctr)
  prop.matrix.case = rdirichlet(N_sample, alpha.case)

  # tmp=dirichlet.mle(out2$Omega[c(1:20,35:47),1:4])  # control proportions
  # alpha.all=tmp$alpha
  # 
  # prop.matrix.ctr = rdirichlet(N_sample, alpha.all)
  # prop.matrix.case = rdirichlet(N_sample, alpha.all)

  prop.matrix.true = rbind(prop.matrix.ctr, prop.matrix.case)

  return(list(prop.matrix.true = prop.matrix.true))
}

### get sample mixture
getSampleMix <- function(N_sample_pergroup=20, p=rep(0.05,Ntissue), n_sd=0.1,
                         profmean=pure_based_mean, prop_n_param=1000, t_noise=1){

  # # get proportion
  tmp = getProportion(N_sample=N_sample_pergroup, prop_n_param=prop_n_param)
  prop.matrix.true = as.matrix(tmp$prop.matrix.true)

  # get proportion
  # tmp = getProportion_old(N_sample=N_sample_pergroup*2, existprop=out2$Omega[,1:4], prop_n_param=prop_n_param)
  # prop.matrix = as.matrix(tmp$prop.matrix.true)

  # get tissue profile
  N_feature = dim(profmean)[1] # number of features
  L = dim(profmean)[2] # number of pure tissues

  # generate true status allowing overlaps
  trueStatus = lfc = matrix(0, N_feature, L)
  for(i in 1:L){
       trueStatus[,i] = rbinom(N_feature, 1, p[i])
       trueStatus[trueStatus[,i]==1,i] = ifelse(runif(sum(trueStatus[,i]==1),0,1)>0.5,2,1)
       lfc[trueStatus[,i]==1,i] = rnorm(sum(trueStatus[,i]==1),-1,sd=0.2)
       lfc[trueStatus[,i]==2,i] = rnorm(sum(trueStatus[,i]==2),1,sd=0.2)
  }

  Y1 = Y2 = matrix(0,N_feature,N_sample_pergroup)
  for(i in 1:N_sample_pergroup){
    # generate pure tissue profile for each person
    tissue.use = getOnePureTissue(profmean, lfc, normal=T, t_noise=t_noise) # in control
    Y1[,i] = t(tissue.use%*%prop.matrix.true[i,])
    tissue.use = getOnePureTissue(profmean, lfc, normal=F, t_noise=t_noise) # in case
    Y2[,i] = t(tissue.use%*%prop.matrix.true[i+N_sample_pergroup,])
  }

  # generate the standard deviation from the real data
  res_sd1 = -8.06 + 0.11*apply(Y1,1,mean)
  res_sd1[res_sd1<0] = 0.01
  res_sd2 = -8.06 + 0.11*apply(Y2,1,mean)
  res_sd2[res_sd2<0] = 0.01

  # add measurement error
  Y1 = Y1 + abs(rnorm(length(c(Y1)),mean=0,sd=res_sd1*n_sd))
  Y2 = Y2 + abs(rnorm(length(c(Y2)),mean=0,sd=res_sd2*n_sd))
  
  # logY1 = log(Y1)
  # logY1[logY1<0] = 0.01
  # logY2 = log(Y2)
  # logY2[logY2<0] = 0.01

  return(list(Y1 = Y1, Y2 = Y2, trueStatus = trueStatus, prop.matrix.true = prop.matrix.true))
}

### get sample mixture
getSampleMix2 <- function(N_sample_pergroup=20, p=rep(0.05,Ntissue), n_sd=1, profmean=pure_based_mean, prop_n_param=50, t_noise=1){
     
     # generate more spread-out tissue proportion
     if(Ntissue==4){
          spreadMatrix = data.frame(tissue1=rep(0.59,5),tissue2=rep(0.3,5),tissue3=rep(0.1,5),tissue4=rep(0.01,5))
     }else if(Ntissue==5){
          spreadMatrix = data.frame(tissue1=rep(0.59,5),tissue2=rep(0.25,5),tissue3=rep(0.1,5),tissue4=rep(0.05,5),tissue5=rep(0.01,5))
     }
     
     # get proportion
     tmp = getProportion_old(N_sample=N_sample_pergroup*2, existprop=spreadMatrix, prop_n_param=prop_n_param)
     prop.matrix.true = as.matrix(tmp$prop.matrix.true)
     prop.matrix = as.matrix(tmp$prop.matrix)
     
     # get tissue profile
     N_feature = dim(profmean)[1] # number of features
     L = dim(profmean)[2] # number of pure tissues
     
     # generate true status allowing overlaps
     trueStatus = lfc = matrix(0, N_feature, L)
     for(i in 1:L){
          trueStatus[,i] = rbinom(N_feature, 1, p[i])
          trueStatus[trueStatus[,i]==1,i] = ifelse(runif(sum(trueStatus[,i]==1),0,1)>0.5,2,1)
          lfc[trueStatus[,i]==1,i] = rnorm(sum(trueStatus[,i]==1),-1,sd=0.2)
          lfc[trueStatus[,i]==2,i] = rnorm(sum(trueStatus[,i]==2),1,sd=0.2)
     }
     
     Y1 = Y2 = matrix(0,N_feature,N_sample_pergroup)
     for(i in 1:N_sample_pergroup){
          # generate pure tissue profile for each person
          tissue.use = getOnePureTissue(profmean, lfc, normal=T, t_noise=t_noise) # in control
          Y1[,i] = t(tissue.use%*%prop.matrix[i,])
          tissue.use = getOnePureTissue(profmean, lfc, normal=F, t_noise=t_noise) # in case
          Y2[,i] = t(tissue.use%*%prop.matrix[i+N_sample_pergroup,])
     }
     
     # generate the standard deviation from the real data
     res_sd1 = -8.06 + 0.11*apply(Y1,1,mean)
     res_sd1[res_sd1<0] = 0.01
     res_sd2 = -8.06 + 0.11*apply(Y2,1,mean)
     res_sd2[res_sd2<0] = 0.01
     
     # add measurement error
     Y1 = Y1 + abs(rnorm(length(c(Y1)),mean=0,sd=res_sd1*n_sd))
     Y2 = Y2 + abs(rnorm(length(c(Y2)),mean=0,sd=res_sd2*n_sd))
     
     return(list(Y1 = Y1, Y2 = Y2, trueStatus = trueStatus, lfc = lfc, prop.matrix.true = prop.matrix))
}

generateReal <- function(pure,N_need){
     ## get proportions
     prop_param = c(0.15,0.5,0.25,0.1)*10
     prop.matrix = rdirichlet(N_need,prop_param)
     ## get mixed
     Y.tmp = pure%*%t(prop.matrix)
     ## add observation noise to mixed
     res_sd = -8.06 + 0.11*apply(Y.tmp,1,mean)
     res_sd[res_sd<0] = 0.01
     
     # add measurement error
     Y = Y.tmp + abs(rnorm(length(c(Y.tmp)),mean=0,sd=res_sd))
     
     return(list(simY = Y, simProp = prop.matrix))
}

# generateReal_spikein <- function(pure,N_need,spikein){
#      ## get proportions
#      prop.matrix = matrix(0,N_need,5)
#      prop_param = c(0.15,0.5,0.25,0.1)*10
#      tmp = rdirichlet(N_need,prop_param)
#      prop.matrix[,1:4] = tmp*(1-spikein)
#      prop.matrix[,5] = spikein
#      
#      ## get mixed
#      Y.tmp = pure%*%t(prop.matrix)
#      ## add observation noise to mixed
#      res_sd = -8.06 + 0.11*apply(Y.tmp,1,mean)
#      res_sd[res_sd<0] = 0.01
#      
#      # add measurement error
#      Y = Y.tmp + abs(rnorm(length(c(Y.tmp)),mean=0,sd=res_sd))
#      
#      return(list(simY = Y, simProp = prop.matrix))
# }

short2long <- function(input = rnd1, col3 = 1){
     longCor <- matrix(0,nrow(input)*ncol(input),3)
     for(i in 1:ncol(input)){
          longCor[(i-1)*nrow(input)+1:nrow(input),1] = input[,i]
          longCor[(i-1)*nrow(input)+1:nrow(input),2] = i
          longCor[(i-1)*nrow(input)+1:nrow(input),3] = col3
     }
     return(longCor)
}

### ROC
makeROC <- function(pval) {
  pred = prediction(-log10(pval), trueDM)
  perf = performance(pred,"tpr","fpr")
  return(perf)
}

calcAUC <- function(pval){
     pred = prediction(-log10(pval),trueDM)
     auc.perf = performance(pred, measure = "auc")
     return(as.numeric(auc.perf@y.values))
}

ProAvgROC <- function(rocs){
     ROC = aveROC(rocs)
     ROC = ROC[order(ROC[,1]),]
     thin = rep(FALSE,nrow(ROC))
     thin[1:200]=T
     thin[seq(201,nrow(ROC),200)]=T
     nROC = ROC[thin,]
     return(nROC)
}

AvgTDR <- function(tdr){
     p = length(tdr[[1]])
     n = length(tdr)
     tmp = matrix(unlist(tdr),n,p,byrow=T)
     avg.tdr = apply(tmp,2,mean)
     return(avg.tdr)
}


### get proportion aligned using correlation coefficient
GetPropAligned <- function(input,reference,L=4){
  colnames(input)=colnames(reference)=c(1,2,3,4)
  corMat = cor(input,reference,use="pairwise.complete.obs")
  prop_cor = rep(0,L)
  tmpmat=corMat
  for(i in 1:L){
    maxind = which(tmpmat == max(tmpmat), arr.ind = TRUE)
    prop_cor[maxind[1]] = colnames(corMat)[maxind[2]]
    tmpmat[maxind[1],]=tmpmat[,maxind[2]]=rep(-1,L)
  }
  colnames(input) = prop_cor
  trans_input = input[,colnames(reference)]
  return(trans_input)
}


#### find index for marker genes based on tissue specific profiles.
## The marker genes should be tissue specific genes.

findRefinx <- function(logpure, nmarker=1000) {
    idx = NULL
    ## rank tissue specificity in each tissue
    for(i in 1:ncol(logpure)) {
        dif = abs(logpure[,i] - logpure[,-i])
        idx[[i]] = sort(rowMeans(dif), dec=TRUE, index=TRUE)$ix
    }
    ## find markers. This is ad hoc, but need to make sure there are balanced number of marker genes for each tissue.
    K = ncol(logpure)
    nmarker.tissue = nmarker/K * 1.2 ## number of markers per tissue. Consider overlaps
    allidx = NULL
    for(i in 1:ncol(logpure)) {
        allidx = c(allidx, idx[[i]][1:nmarker.tissue])
    }
    allidx = unique(allidx) ## number of markers could be a bit off from nmarker, but it's ok
    allidx
}


findRefinx_ct <- function(logpure, nmarker=1000, ct=1:ncol(logpure)) {
     idx = NULL
     ## rank tissue specificity in each tissue
     for(i in 1:ncol(logpure)) {
          dif = abs(logpure[,i] - logpure[,-i])
          idx[[i]] = sort(rowMeans(dif), dec=TRUE, index=TRUE)$ix
     }
     ## find markers. This is ad hoc, but need to make sure there are balanced number of marker genes for each tissue.
     K = ncol(logpure)
     nmarker.tissue = nmarker/K * 1.2 ## number of markers per tissue. Consider overlaps
     allidx = NULL
     for(i in ct) {
          allidx = c(allidx, idx[[i]][1:nmarker.tissue])
     }
     allidx = unique(allidx) ## number of markers could be a bit off from nmarker, but it's ok
     allidx
}

## find index for marker genes based on raw data -- largest cv

findRefinx.cv <- function(rawdata, nmarker=1000, startn=0) {
    mm = rowMeans(rawdata)
    vv = rowVars(rawdata)
    cv = sqrt(vv) / mm
    ##cv = vv
    cv[is.na(cv)] = 0
    ix = sort(cv, dec=TRUE, index=TRUE)$ix
    ix[startn+1:nmarker]
}

## find index for marker genes based on raw data -- largest var(log(data))
findRefinx.var <- function(rawdata, nmarker=1000, startn=0) {
    vv = rowVars(log(rawdata+1))
    vv[is.na(vv)] = 0
    ix = sort(vv, dec=TRUE, index=TRUE)$ix
    ix[startn+1:nmarker]
}

## smallest score means the highest ranked
makeTDR <- function(score, trueDM, steps=seq(100, 2000, by=100)) {
    score[is.na(score)] = 10000
    ## find top ranked ones
    ix = sort(score, index=TRUE)$ix
    res = rep(NA, length(steps))
    for(i in 1:length(steps) ) {
        thisix = ix[1:steps[i]]
        res[i] = mean(trueDM[thisix])
    }
    res
}

## a function to estimate proportions
getProp <- function(Y1.raw, Y2.raw, refinx, ref, 
                    reference1=out$prop.matrix.true[1:N_sample,], 
                    reference2=out$prop.matrix.true[N_sample+1:N_sample,]) {
     Y1 = Y1.raw[refinx,]
     Y2 = Y2.raw[refinx,]
     outY1 = ged(Y1,ref)
     Prop1_raw = t(coef(outY1))
     Prop1 = GetPropAligned(input=Prop1_raw,reference=reference1)
     outY2 = ged(Y2,ref)
     Prop2_raw = t(coef(outY2))
     Prop2 = GetPropAligned(input=Prop2_raw,reference=reference2)
     rbind(Prop1, Prop2)
}

## a function to estimate proportions
getProp_single <- function(Y.raw, refinx, ref, reference) {
     Y = Y.raw[refinx,]
     outY = ged(Y,ref)
     Prop_raw = t(coef(outY))
     Prop = GetPropAligned(input=Prop_raw,reference=reference)
     return(Prop)
}

getProp_single.RF <- function(Y.raw, refinx, ref, reference) {
     Y = Y.raw[refinx,]
     outY = ged(Y,K)
     Prop_raw = t(coef(outY))
     Prop = GetPropAligned(input=Prop_raw,reference=reference)
     rbind(Prop)
}

callDE <- function(prop,k) {
     Prop1 = prop[1:N_sample,]
     Prop2 = prop[N_sample+(1:N_sample),]
     Y = cbind(Y1.raw,Y2.raw)
     design = c(rep(0,N_sample),rep(1,N_sample))
     tmp = DEKTissue(K, Y, Prop=rbind(Prop1,Prop2), design=design, WhichPar=k+K)
     KDE.pval = as.numeric(tmp$t.pval)
     
     rocs1 = makeROC(KDE.pval)
     steps=seq(100, 3000, by=100)
     tdr1 = makeTDR(log(KDE.pval), trueDM, steps)
     
     list(ROC=rocs1, TDR=tdr1)
}

DEprop <- function(prop0,refinx0, niter, method="RB",
                   refinx.tissue = refinx.tissue, 
                   propTRUE = out$prop.matrix.true[1:N_sample,]){
     
     corr = matrix(0, nrow=niter+1, ncol=K)
     nOverlap = rep(0, niter+1)
     
     corr[1,] = diag(cor(prop0, propTRUE))
     nOverlap[1] = length(intersect(refinx0, refinx.tissue))
     
     proportion=prop0
     
     for(i in 1:niter) {
          ## find tissue specific genes
          idx = NULL
          for(k in 1:K) {
               cvec = rep(-1/(K-1),K)
               cvec[k] = 1
               design = rep(0,N_sample)
               tmp = DEKTissue(K, Y=Y1.raw, Prop=proportion, design=design, contrast.vec=cvec)
               idx[[k]] = sort(abs(tmp$t.stat), dec=TRUE, index=TRUE)$ix
          }
          nmarker = nMarker
          nmarker.tissue = nmarker/K * 1.2 ## number of markers per tissue. Consider overlaps
          idxMarker = NULL
          for(k in 1:K) {
               idxMarker = c(idxMarker, idx[[k]][1:nmarker.tissue])
          }
          idxMarker = unique(idxMarker)
          nOverlap[i+1] = length(intersect(idxMarker, refinx.tissue))
          
          ## redo deconvolution using new markers
          thisRef = referSig[idxMarker,]
          if(method == "RB") ## reference based
               prop = getProp(Y1.raw, Y2.raw, idxMarker, thisRef)
          else ## reference free
               prop = getProp.RF(Y1.raw, Y2.raw, idxMarker, thisRef)
          
          corr[i+1,] = diag(cor(prop, out$prop.matrix.true))
          proportion = prop[1:N_sample,]
     }
     
     return(list(corr,nOverlap,prop))
}

DEprop_single <- function(Y, prop0,refinx0, niter, method="RB",
                   refinx.tissue = refinx.tissue, 
                   propTRUE = out$prop.matrix.true[1:N_sample,]){
     
     corr = matrix(0, nrow=niter+1, ncol=K)
     nOverlap = rep(0, niter+1)
     
     corr[1,] = diag(cor(prop0, propTRUE))
     nOverlap[1] = length(intersect(refinx0, refinx.tissue))
     
     proportion=prop0
     
     for(i in 1:niter) {
          ## find tissue specific genes
          idx = NULL
          for(k in 1:K) {
               cvec = rep(-1/(K-1),K)
               cvec[k] = 1
               design = rep(0,N_sample)
               tmp = DEKTissue(K, Y=Y, Prop=proportion, design=design, contrast.vec=cvec)
               idx[[k]] = sort(abs(tmp$t.stat), dec=TRUE, index=TRUE)$ix
          }
          nmarker = nMarker
          nmarker.tissue = nmarker/K * 1.2 ## number of markers per tissue. Consider overlaps
          idxMarker = NULL
          for(k in 1:K) {
               idxMarker = c(idxMarker, idx[[k]][1:nmarker.tissue])
          }
          idxMarker = unique(idxMarker)
          nOverlap[i+1] = length(intersect(idxMarker, refinx.tissue))
          
          ## redo deconvolution using new markers
          thisRef = referSig[idxMarker,]
          if(method == "RB") ## reference based
               prop = getProp_single(Y, idxMarker, thisRef, reference=propTRUE)
          else ## reference free
               prop = getProp_single.RF(Y, idxMarker, thisRef, reference=propTRUE)
          
          corr[i+1,] = diag(cor(prop, propTRUE))
          proportion = prop
     }
     
     return(list(corr,nOverlap,prop))
}
