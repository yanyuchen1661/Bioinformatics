# Differential expression with tissue proportion
# K tissues

DEK_wrapper <- function(K, Y, Prop, design, sort=F, var.threshold=0.1, data.threshold=TRUE){
     DEK.res <- list()

     if(!all(design==0)){
          listname = rep("nn",K)
          for(k in 1:K){
               DEK.res[[k]] = DEKTissue(K, Y, Prop, design, sort=sort, var.threshold=var.threshold, data.threshold=data.threshold, WhichPar=k+K)
               if(is.null(colnames(Prop))){
                    listname[k] = paste0("Tissue",k)
               }else if(!is.null(colnames(Prop))){
                    listname[k] = colnames(Prop)[k]
               }
          }

     }

     else if(all(design==0)){
          g=1
          listname = rep("nn",K*(K-1)/2)
          for(i in 1:(K-1)){
               for(j in 2:K){
                    cc = rep(0,K)
                    cc[i] = 1
                    cc[j] = -1
                    DEK.res[[g]] = DEKTissue(K, Y, Prop, design, sort=sort, var.threshold=var.threshold, data.threshold=data.threshold, contrast.vec = cc)
                    listname[g] = paste0("Tissue",i,"vs",j)
                    g = g+1
               }
          }
     }
     names(DEK.res) = listname
     return(DEK.res)
}


DEKTissue <- function(K, Y, Prop, design, WhichPar=NULL, contrast.vec=NULL, sort=F, var.threshold=0.1, data.threshold=TRUE){

     # K is the number of tissue types
     # Y is a G*N matrix, G is the number of features, N is the number of subjects
     # Prop is a N*K matrix, tissue proportions for all subjects
     # design is a design matrix representing the status of samples (or continuous features: need to be developed)
     #              example 1 of design: (0,...,0,1,...,1), 0 for control and 1 for cases
     #              example 2 of design: (0,...,0) for all subjects belong to one group

     # WhichPar is a number chosen from 1 to 2K (which parameter will be tested?)
     # contrast.vec if in use, WhichPar will be ignored.  contrast.vec should be a 2K*1 vector to specify a contrast vector
     # consisted of 1/-1.
     # var.threshold controls how much you want to bound your variation estimation. Default value is 0.1.

     N = dim(Y)[2]
     if(dim(Prop)[1]!=N | dim(Prop)[2]!=K){
          stop("Dimension of proportion input is not correct!")
     }

     Y = t(na.omit(Y))
     G = dim(Y)[2]

     if(!all(design==0)){
          W = cbind(Prop,Prop*design)
     }else{
          W = Prop
     }
     H = solve(t(W)%*%W)%*%t(W)
     coefs = H%*%Y
     Ypred = W%*%coefs
     resi = Y-Ypred

     s2.case = colSums(resi^2) / (N - ncol(W))
     varBeta = matrix(diag(solve(t(W)%*%W)),ncol(W),1)%*%s2.case
     
     ## bound varBeta a bit
     if(data.threshold) {
          var.threshold = quantile(c(varBeta),0.1)
     }else {
          varBeta[varBeta<var.threshold] = var.threshold
     }

     res.table = data.frame(t.stat=rep(0,G), t.pval=rep(0,G), t.fdr=rep(0,G))
     rownames(res.table) = colnames(Y)

     if(!all(design==0)){
          if(is.null(contrast.vec)){
               res.table = data.frame(beta=rep(0,G), mu=rep(0,G), effect_size=rep(0,G), t.stat=rep(0,G), t.pval=rep(0,G), t.fdr=rep(0,G))
               rownames(res.table) = colnames(Y)
               res.table$beta=coefs[WhichPar,]
               res.table$mu=coefs[WhichPar-K,]
               res.table$effect_size = res.table$beta/(res.table$mu + res.table$beta/2)
               res.table$t.stat=coefs[WhichPar,]/sqrt(varBeta[WhichPar,])
               res.table$t.pval=2*pt(-abs(res.table$t.stat),df=N-ncol(W))
               res.table$t.fdr=p.adjust(res.table$t.pval, method = 'fdr')
          }else{
               res.table = data.frame(t.stat=rep(0,G), t.pval=rep(0,G), t.fdr=rep(0,G))
               rownames(res.table) = colnames(Y)
               res.table$t.stat=as.numeric(contrast.vec%*%coefs)/as.numeric(sqrt(abs(contrast.vec)%*%varBeta))
               res.table$t.pval=2*pt(-abs(res.table$t.stat),df=N-ncol(W))
               res.table$t.fdr=p.adjust(res.table$t.pval, method = 'fdr')
          }
     }else{

          if(is.null(contrast.vec)){
               res.table = data.frame(t.stat=rep(0,G), t.pval=rep(0,G), t.fdr=rep(0,G))
               rownames(res.table) = colnames(Y)
               res.table$t.stat=coefs[WhichPar,]/sqrt(varBeta[WhichPar,])
               res.table$t.pval=2*pt(-abs(res.table$t.stat),df=N-ncol(W))
               res.table$t.fdr=p.adjust(res.table$t.pval, method = 'fdr')
          }else{
               res.table = data.frame(muA=rep(0,G), muB=rep(0,G), t.stat=rep(0,G), t.pval=rep(0,G), t.fdr=rep(0,G))
               rownames(res.table) = colnames(Y)
               i=which(contrast.vec!=0)[1]
               j=which(contrast.vec!=0)[2]
               res.table$muA = coefs[i,]
               res.table$muB = coefs[j,]
               res.table$t.stat=as.numeric(contrast.vec%*%coefs)/as.numeric(sqrt(abs(contrast.vec)%*%varBeta))
               res.table$t.pval=2*pt(-abs(res.table$t.stat),df=N-ncol(W))
               res.table$t.fdr=p.adjust(res.table$t.pval, method = 'fdr')
          }
     }

     if(sort){
          return(res.table[order(res.table$t.pval),])
     }else{
          return(res.table)
     }
}
