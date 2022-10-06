#This method is from Cox, D. R., and H. S. Battey. "Large numbers of explanatory variables, a semi-descriptive analysis." Proceedings of the National Academy of Sciences 114.32 (2017): 8592-8595.
#This is useful when we have huge number of variables but few samples.

library(HTRX)
library(HCmodelSets)
n=300
assesshap=FALSE
usebinary=TRUE
nsnp=7

set.seed(123)
if(assesshap){
  ### assessing haplotypes
  hap1=as.data.frame(matrix(0,nrow=n,ncol=nsnp))
  hap2=as.data.frame(matrix(0,nrow=n,ncol=nsnp))
  p=runif(nsnp,0.1,0.9)
  for(j in 1:nsnp){
    hap1[,j]=rbinom(n,1,p[j])
    hap2[,j]=rbinom(n,1,p[j])
  }
  hap=make_htrx(hap1,hap2)
}else{
  ### assessing SNPs
  hap1=as.data.frame(matrix(0,nrow=n,ncol=3^nsnp-1))
  hap2=as.data.frame(matrix(0,nrow=n,ncol=3^nsnp-1))
  p=runif((3^nsnp-1),0.05,0.95)
  for(j in 1:(3^nsnp-1)){
    hap1[,j]=rbinom(n,1,p[j])
    hap2[,j]=rbinom(n,1,p[j])
  }
  
  hap=hap1+hap2
}


#parameters for haplotypes/SNPs
para=rep(0,ncol(hap))

#select 6 haplotypes/SNPs with actual effect
effect_hap=sample(1:ncol(hap),size=6,replace=FALSE)

effect_hap=sort(effect_hap)

# Bigger effect size leads to better results!
para[effect_hap]=runif(6,-1,1)
#para[effect_hap]=runif(6,0.5,1)
#para[effect_hap]=runif(6,-0.5,0.5)

y=as.matrix(hap) %*% matrix(para,ncol=1)+rnorm(n,0,0.1)
if(usebinary){
  y_p=exp(y)/(1+exp(y))
  y=rbinom(length(y),1,prob=y_p)
}


split=list()
#For each analysis, the remaining variables should be around 10-20
#So we need to change the decision rules for each
for(i in 1:10){
  print(i)
  if(usebinary&&assesshap){
    #vector.signif refers to the decision rules used at different stages of reduction
    #vector.signif=c(2,0.05,0.05) for more variables
    split[[i]]=Reduction.Phase(as.matrix(hap[1:210,]),y[1:210],family=binomial,
                               seed.HC=i,dmHC=3,vector.signif=c(2,0.05))
  }
  if(usebinary&&!assesshap){
    #vector.signif=c(2,0.01,0.005) for more variables
    split[[i]]=Reduction.Phase(as.matrix(hap[1:210,]),y[1:210],family=binomial,
                               seed.HC=i,dmHC=3,vector.signif=c(2,0.005))
  }
  if(!usebinary&&assesshap){
    #vector.signif=c(2,0.0001,0.0001) for more variables
    split[[i]]=Reduction.Phase(as.matrix(hap[1:210,]),y[1:210],family=gaussian,
                               seed.HC=i,dmHC=3,vector.signif=c(2,0.0001))
  }
  if(!usebinary&&!assesshap){
    #vector.signif=c(2,0.005,0.002) for more variables
    split[[i]]=Reduction.Phase(as.matrix(hap[1:210,]),y[1:210],family=gaussian,
                               seed.HC=i,dmHC=3,vector.signif=c(2,0.002))
  }
}

v <- list()
#in the paper they even repeated 1000 times!
for(i in 1:10){
  v[[i]]=sort(split[[i]]$List.Selection$`Hypercube with dim 2`$numSelected1)
}

v_all=as.data.frame(table(unlist(v)))
if(length(which(v_all[,1]==0))==1){
  v_all=v_all[-which(v_all[,1]==0),]
}

# use the common variables selected in at least "threshold" analysis 
# threshold can be 1 to 10, the paper suggests 5
threshold=5

v_remain=as.numeric(as.character(v_all[which(v_all$Freq>=threshold),1]))

# assess whether the haplotypes with actual effects are correctly selected
list(index_selected=v_remain,index_real_effect=effect_hap)

data_true=data.frame(y,hap[,effect_hap])
data_infer=data.frame(y,hap[,v_remain])
if(usebinary){
  model_true=glm(y~.,data_true,family=binomial)
  model_infer=glm(y~.,data_infer,family=binomial)
  list(R2_true=computeR2(mypredict(model_true,data_true),data_true[,"y"],usebinary=2),
       R2_infer=computeR2(mypredict(model_infer,data_infer),data_infer[,"y"],usebinary=2))
}else{
  model_true=summary(lm(y~.,data_true))
  model_infer=summary(lm(y~.,data_infer))
  list(R2_true=model_true$r.squared,R2_infer=model_infer$r.squared)
}





#model selection
if(usebinary){
  out = ModelSelection.Phase(X=as.matrix(hap[211:300,]),Y=y[211:300],family=binomial,
                             list.reduction = v_remain, signif = 0.05,modelSize=7)
}else{
  out = ModelSelection.Phase(X=as.matrix(hap[211:300,]),Y=y[211:300],family=gaussian,
                             list.reduction = v_remain, signif = 0.05,modelSize=7) 
}

#the correct model should be reported in "out" if their effect size is big
out


#assess performance
allindex = unlist(out)
total_model=0
for(i in 1:length(out$goodModels)){
  total_model=total_model+length(out$goodModels[[i]])/i
}
correct_proportion=data.frame(index=effect_hap,proportion=rep(0,6))
for(i in 1:6){
  correct_proportion[i,2]=length(allindex[allindex==effect_hap[i]])/total_model
}
correct_proportion
